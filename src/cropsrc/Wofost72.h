/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _WOFOST72_H_
#define _WOFOST72_H_

#include "../define.h"
#include "CropState.h"        // Static crop parameters
#include "CropDynamicState.h" // Dynamic crop states (DVS, LAI, Biomass etc.)
#include "MeteoState.h"       // Weather data for driving

#include "Phenology.h"                 // Phenology module
#include "Assimilation.h"              // Assimilation module
#include "AssimilateState.h"           // Assimilation output states
#include "Respiration.h"               // Respiration module
#include "RespirationState.h"          // Respiration output states
#include "Partitioning.h"              // Partitioning module
#include "PartitioningState.h"         // Partitioning output states
#include "RootDynamics.h"              // Root dynamics module
#include "RootDynamicsState.h"         // Root dynamics output states
#include "StemDynamics.h"              // Stem dynamics module
#include "StemDynamicsState.h"         // Stem dynamics output states
#include "StorageOrganDynamics.h"      // Storage organ dynamics module
#include "StorageOrganDynamicsState.h" // Storage organ dynamics output states
#include "LeafDynamics.h"              // Leaf dynamics module
#include "LeafDynamicsState.h"         // Leaf dynamics output states
#include "Evapotranspiration.h"        // Evapotranspiration module
#include "EvapotranspirationState.h"   // Evapotranspiration output states

#include <cmath>
#include <algorithm> // For fmin, fmax

class Wofost72
{

public:
    // Instance of sub-modules
    Phenology pheno;
    Assimilation assim;
    Respiration mres;
    Evapotranspiration evtra; // Note: This is our C++ simplified ET module
    Partitioning part;
    RootDynamics ro_dynamics;
    StemDynamics st_dynamics;
    StorageOrganDynamics so_dynamics;
    LeafDynamics lv_dynamics;

    // All dynamic states
    CropDynamicState s;             // General crop states (DVS, LAI, total biomass etc.)
    AssimilateState as;             // Assimilation specific outputs
    RespirationState rs;            // Respiration specific outputs
    PartitioningState ps;           // Partitioning fractions
    RootDynamicsState rds;          // Root dynamics states
    StemDynamicsState sds;          // Stem dynamics states
    StorageOrganDynamicsState sods; // Storage organ dynamics states
    LeafDynamicsState lds;          // Leaf dynamics states
    EvapotranspirationState ets;    // Evapotranspiration states

    // Main WOFOST state variables
    realArr TAGP;  // Total above-ground Production [kg ha-1]
    realArr GASST; // Total gross assimilation [kg CH2O ha-1]
    realArr MREST; // Total gross maintenance respiration [kg CH2O ha-1]
    realArr CTRAT; // Total crop transpiration accumulated [cm]
    realArr CEVST; // Total soil evaporation accumulated [cm]
    realArr HI;    // Harvest Index [-]
    // DOF, FINISH_TYPE would be flags or dates per cell, for now simplified
    intArr CROP_FINISHED_FLAG; // 0=running, 1=finished
    realArr CROP_FINISH_DVS;   // DVS at which crop finished

    // Default Constructor
    Wofost72()
    {
        // Do nothing, wait for allocate()
    }

    // Allocation method
    void allocate(int nCells)
    {
        // Allocate all dynamic state arrays
        s.allocate(nCells);
        as.allocate(nCells);
        rs.allocate(nCells);
        ps.allocate(nCells);
        rds.allocate(nCells);
        sds.allocate(nCells);
        sods.allocate(nCells);
        lds.allocate(nCells);
        ets.allocate(nCells);

        TAGP = realArr("TAGP", nCells);
        GASST = realArr("GASST", nCells);
        MREST = realArr("MREST", nCells);
        CTRAT = realArr("CTRAT", nCells);
        CEVST = realArr("CEVST", nCells);
        HI = realArr("HI", nCells);
        CROP_FINISHED_FLAG = intArr("CROP_FINISHED_FLAG", nCells);
        CROP_FINISH_DVS = realArr("CROP_FINISH_DVS", nCells);
    }

    // Main initialization function
    // p: Static Crop Parameters (CropState)
    // initial_root_zone_SM: initial average volumetric soil moisture in root zone from SERGHEI
    void initialize(const CropState &p, const MeteoState &m, const realArr &initial_root_zone_SM, int current_day_of_year)
    {

        // --- 1. Initialize Phenology ---
        pheno.initialize(s, p);

        // --- 2. Initialize Partitioning (to get initial FR, FL, FS, FO) ---
        // Needs initial DVS from phenology
        part.initialize(ps, p, s);

        // --- CRITICAL FIX: Sync partitioning factors BEFORE organ dynamics ---
        Kokkos::parallel_for("Wofost72_SyncPartitioning", s.nCells, KOKKOS_LAMBDA(const int i) {
            s.FL(i) = ps.FL(i);
            s.FR(i) = ps.FR(i);
// Note: FS and FO are not stored in CropDynamicState, only in PartitioningState

// DEBUG: Print partitioning factors
#if DEBUG_CROP_GROWTH_MODEL
            printf("[WOFOST-INIT] cell=%d, Partitioning - FR=%.4f, FL=%.4f, FS=%.4f, FO=%.4f\n",
                   i, ps.FR(i), ps.FL(i), ps.FS(i), ps.FO(i));
#endif
        });

        // --- 3. Initialize Organ Dynamics (needs initial FR, FL, FS, FO) ---
        // Now s.FL, s.FR, s.FS, s.FO are properly set

        ro_dynamics.initialize(rds, p, s);
        ro_dynamics.set_initial_biomass(rds, p, ps.FR); // Initial WRT based on TDWI and FR

        st_dynamics.initialize(sds, p, s, ps.FR, ps.FS); // Initial WST and SAI
        so_dynamics.initialize(sods, p, ps.FR, ps.FO);   // Initial WSO and PAI
        // lv_dynamics.initialize(lds, p, s); // Initial WLV, LAIEM, LASUM, LAIEXP, LAIMAX, LAI
        lv_dynamics.initialize(lds, p, s, ps.FR, ps.FL, sds.SAI, sods.PAI);
        // --- 4. Initialize Evapotranspiration (needs initial LAI) ---
        // No explicit init needed for our simplified ET module, as it computes rates based on inputs.

        // --- 5. Initialize other states and sync from sub-modules ---
        Kokkos::parallel_for("Wofost72_MainInit", s.nCells, KOKKOS_LAMBDA(const int i) {
            // Sync biomass states from sub-modules to main state
            s.WLV(i) = lds.WLV(i);
            s.DWLV(i) = lds.DWLV(i);
            s.TWLV(i) = lds.TWLV(i);
            s.LAI(i) = lds.LAI(i);
            
            s.WST(i) = sds.WST(i);
            s.DWST(i) = sds.DWST(i);
            s.TWST(i) = sds.TWST(i);
            s.SAI(i) = sds.SAI(i);
            
            s.WRT(i) = rds.WRT(i);
            s.DWRT(i) = rds.DWRT(i);
            s.TWRT(i) = rds.TWRT(i);
            s.RD(i) = rds.RD(i);
            
            s.WSO(i) = sods.WSO(i);
            s.DWSO(i) = sods.DWSO(i);
            s.TWSO(i) = sods.TWSO(i);
            s.PAI(i) = sods.PAI(i);

            // Initial total above-ground biomass
            TAGP(i) = lds.TWLV(i) + sds.TWST(i) + sods.TWSO(i);
            GASST(i) = 0.0;
            MREST(i) = 0.0;
            CTRAT(i) = 0.0;
            CEVST(i) = 0.0;
            HI(i) = 0.0;
            CROP_FINISHED_FLAG(i) = 0;
            CROP_FINISH_DVS(i) = 0.0;

            // Check initial biomass partitioning (optional, for debugging)
            real checksum = p.p.TDWI - TAGP(i) - rds.TWRT(i);
            if (std::abs(checksum) > 0.0001) {
                printf("Error in partitioning of initial biomass (TDWI) at cell %d! Checksum: %f\n", i, checksum);
            } });
    }

    // Main daily step function for WOFOST
    // current_day_of_year: from 1 to 365/366
    // current_time_s: Current simulation time in seconds (for meteo index)
    // root_zone_SM: average volumetric soil moisture from SERGHEI for the day
    void calc_rates_and_integrate(const CropState &p, const MeteoState &m,
                                  real current_time_s, const realArr &root_zone_SM, int current_day_of_year, real delt_day = 1.0)
    {

        int drv_idx = m.get_index_at_time(current_time_s);

#if DEBUG_CROP_GROWTH_MODEL

        if (s.nCells > 0)
        {
            printf("[WOFOST] Step: DOY=%d, Time=%.1f, drv_idx=%d\n", current_day_of_year, current_time_s, drv_idx);
            printf("[WOFOST] Meteo: TMIN=%.1f, TMAX=%.1f, IRRAD=%.1f, ET0=%.1f\n",
                   m.tmin(drv_idx), m.tmax(drv_idx), m.irrad(drv_idx), m.et0(drv_idx));
            printf("[WOFOST] Before calc: LAI=%.2f, WRT=%.2f, WST=%.2f, WSO=%.2f, WLV=%.2f, DVS=%.2f\n",
                   s.LAI(0), s.WRT(0), s.WST(0), s.WSO(0), s.WLV(0), s.DVS(0));
        }

#endif

        // --- 1. Phenology ---
        pheno.calc_rates(s, p, m, drv_idx, current_day_of_year);

#if DEBUG_CROP_GROWTH_MODEL
        if (s.nCells > 0)
        {
            printf("[WOFOST-PHENO] After calc_rates: STAGE=%d, DVS=%.2f, TSUM=%.2f\n",
                   (int)s.STAGE(0), s.DVS(0), s.TSUM(0));
        }
#endif

        // Check if crop has emerged. If not, skip most calculations.
        Kokkos::parallel_for("Wofost72_CheckEmergence", s.nCells, KOKKOS_LAMBDA(const int i) {
            if (s.STAGE(i) == CropStage::EMERGING) {
                s.skip_other_modules(i) = 1; // Mark this cell to skip other modules
            } else {
                s.skip_other_modules(i) = 0;
            } });
#if DEBUG_CROP_GROWTH_MODEL
        // Debug: Check if skipping
        if (s.nCells > 0)
        {
            printf("[WOFOST-EMERG] skip_other_modules=%d\n", s.skip_other_modules(0));
        }
#endif

        // --- 2. Assimilation ---
        assim.calc_rates(as, p, s, m, drv_idx, current_day_of_year);

        // --- 3. Evapotranspiration ---
        evtra.calc_rates(ets, p, s, m, drv_idx, root_zone_SM);

// Debug: Check RFTRA and PGASS
#if DEBUG_CROP_GROWTH_MODEL
        if (s.nCells > 0)
        {
            printf("[EVAP-DEBUG] RFTRA=%.6f, TRA=%.6f, TRAMX=%.6f\n",
                   ets.RFTRA(0), ets.TRA(0), ets.TRAMX(0));
            printf("[ASSIM-DEBUG] PGASS=%.6f\n", as.PGASS(0));
        }
#endif

        // --- 4. Partitioning Rates (needed before carbon loop) ---
        part.update(ps, p, s); // Update partitioning factors based on current DVS

        // --- Calculate Respiration Rates (outside Kokkos kernel) ---
        mres.calc_rates(rs, p, s, m, drv_idx);

        // --- Main Carbon Balance Loop ---
        Kokkos::parallel_for("Wofost72_CarbonLoop", s.nCells, KOKKOS_LAMBDA(const int i) {
            if (s.skip_other_modules(i) == 1)
                return; // Skip if still emerging

            // Water stress reduction on assimilation
            // r.GASS = PGASS * k.RFTRA
            s.GASS(i) = as.PGASS(i) * ets.RFTRA(i);

            // Respiration (already calculated above)
            // r.MRES  = min(r.GASS, PMRES)
            s.MRES(i) = fmin(s.GASS(i), rs.PMRES(i));

            // Net available assimilates
            // r.ASRC  = r.GASS - r.MRES
            s.ASRC(i) = s.GASS(i) - s.MRES(i);

            // DM partitioning factors (pf)
            // part.calc_rates(day, drv) -> ps.FR, ps.FL, ps.FS, ps.FO
            real FR = ps.FR(i);
            real FL = ps.FL(i);
            real FS = ps.FS(i);
            real FO = ps.FO(i);

            // Conversion factor (CVF)
            // CVF = 1./((pf.FL/p.CVL + pf.FS/p.CVS + pf.FO/p.CVO) * (1.-pf.FR) + pf.FR/p.CVR)
            real den = (FL / p.p.CVL + FS / p.p.CVS + FO / p.p.CVO) * (1.0 - FR) + FR / p.p.CVR;
            real CVF = (den != 0.0) ? (1.0 / den) : 0.0;

            // Total dry matter increase (DMI)
            // r.DMI = CVF * r.ASRC
            s.DMI(i) = CVF * s.ASRC(i);

            // Check carbon balance (simplified)
            real py_checksum_den = fmax(0.0001, s.GASS(i));
            real py_checksum = (s.GASS(i) - s.MRES(i) - (FR + (FL + FS + FO) * (1.0 - FR)) * s.DMI(i) / CVF) * 1.0 / py_checksum_den;
            if (std::abs(py_checksum) >= 0.0001)
            {
                printf("Carbon balance error at cell %d, checksum: %f\n", i, py_checksum);
            }

            // Reallocation (WOFOST72: always 0)
            s.REALLOC_LV(i) = 0.0;
            s.REALLOC_ST(i) = 0.0;
            s.REALLOC_SO(i) = 0.0;

            // Above-ground dry matter increase
            // r.ADMI = (1. - pf.FR) * r.DMI
            s.ADMI(i) = (1.0 - FR) * s.DMI(i);

// DEBUG: Print carbon balance details for first cell
#if DEBUG_CROP_GROWTH_MODEL
            if (i == 0)
            {
                printf("[CARBON-LOOP] cell=%d, PGASS=%.6f, RFTRA=%.4f, GASS=%.6f, PMRES=%.6f, MRES=%.6f, ASRC=%.6f\n",
                       i, as.PGASS(i), ets.RFTRA(i), s.GASS(i), rs.PMRES(i), s.MRES(i), s.ASRC(i));
                printf("[CARBON-LOOP] cell=%d, FR=%.4f, FL=%.4f, FS=%.4f, FO=%.4f, CVF=%.6f, DMI=%.6f, ADMI=%.6f\n",
                       i, FR, FL, FS, FO, CVF, s.DMI(i), s.ADMI(i));
            }

#endif
        });

        // --- Call Organ Dynamics calc_rates ---
        ro_dynamics.calc_rates(rds, p, s, s.DMI, ps.FR);
        st_dynamics.calc_rates(sds, p, s, s.ADMI, ps.FS, s.REALLOC_ST);
        so_dynamics.calc_rates(sods, p, s.ADMI, ps.FO, s.REALLOC_SO);
        lv_dynamics.calc_rates(lds, p, s, m, drv_idx, s.ADMI, ets);

        // --- Integrate All States ---
        pheno.integrate(s, p, delt_day);

        // Skip integration if crop emerged just now, or was already emerging
        // Note: pheno.integrate updates STAGE. If it was EMERGING and becomes VEGETATIVE,
        // we might want to start other modules next step. Or check previous state.
        // Python code: crop_stage = pheno.get_variable("STAGE") (before integration)
        // Here we rely on s.skip_other_modules set earlier based on state BEFORE integration of this step.

        part.update(ps, p, s); // Partitioning updated based on DVS (new DVS from pheno integrate?)
        // In python: pheno.integrate -> part.integrate (which updates FR based on NEW DVS). Correct.
        // BUT: Python `integrate` uses `rates` calculated with `OLD` state.
        // Here: `ro_dynamics.integrate` uses `rds` (rates) calculated above.
        // `rds` calculated using OLD FR/DVS.
        // `pheno.integrate` updates DVS.
        // `part.update` updates FR using NEW DVS?
        // Python: part.integrate updates state (FR).
        // Logic: calc_rates (all modules) -> integrate (all modules).
        // Standard Euler integration: State_new = State_old + Rate(State_old) * dt.
        // So update of FR should happen AFTER rate calculation for next step?
        // Or FR is a state variable updated in integrate. Yes.

        ro_dynamics.integrate(rds, delt_day);
        st_dynamics.integrate(sds, p, s, delt_day);
        so_dynamics.integrate(sods, p, delt_day);
        lv_dynamics.integrate(lds, s, delt_day); // lds depends on s.SAI and s.PAI

        // --- Update Main WOFOST States ---
        Kokkos::parallel_for("Wofost72_MainIntegrate", s.nCells, KOKKOS_LAMBDA(const int i) {
            if (s.skip_other_modules(i) == 1 && s.STAGE(i) == CropStage::EMERGING) {
                return;
            }

            // Sync biomass states from sub-modules to main state (CRITICAL FIX)
            s.WLV(i) = lds.WLV(i);
            s.DWLV(i) = lds.DWLV(i);
            s.TWLV(i) = lds.TWLV(i);
            s.LAI(i) = lds.LAI(i);
            
            s.WST(i) = sds.WST(i);
            s.DWST(i) = sds.DWST(i);
            s.TWST(i) = sds.TWST(i);
            s.SAI(i) = sds.SAI(i);
            
            s.WRT(i) = rds.WRT(i);
            s.DWRT(i) = rds.DWRT(i);
            s.TWRT(i) = rds.TWRT(i);
            s.RD(i) = rds.RD(i);
            
            s.WSO(i) = sods.WSO(i);
            s.DWSO(i) = sods.DWSO(i);
            s.TWSO(i) = sods.TWSO(i);
            s.PAI(i) = sods.PAI(i);

            // Integrate total (living+dead) above-ground biomass
            TAGP(i) = lds.TWLV(i) + sds.TWST(i) + sods.TWSO(i);

           

            // Total gross assimilation and maintenance respiration
            GASST(i) += s.GASS(i) * delt_day;
            MREST(i) += s.MRES(i) * delt_day;
            
            // Total crop transpiration and soil evaporation
            CTRAT(i) += ets.TRA(i) * delt_day;
            CEVST(i) += ets.EVS(i) * delt_day; 
            
            // Check for crop maturity
            if (s.STAGE(i) == CropStage::MATURE && CROP_FINISHED_FLAG(i) == 0) {
                CROP_FINISHED_FLAG(i) = 1;
                CROP_FINISH_DVS(i) = s.DVS(i);
            } });

// Debug: Check final results after integration
#if DEBUG_CROP_GROWTH_MODEL
        if (s.nCells > 0)
        {
            printf("[WOFOST-RESULT] After integrate: LAI=%.2f, WRT=%.2f, WST=%.2f, WSO=%.2f, WLV=%.2f, DVS=%.2f\n",
                   s.LAI(0), s.WRT(0), s.WST(0), s.WSO(0), s.WLV(0), s.DVS(0));
            printf("[WOFOST-RESULT] Rates: GASS=%.2f, MRES=%.2f, DMI=%.2f, ADMI=%.2f\n",
                   s.GASS(0), s.MRES(0), s.DMI(0), s.ADMI(0));
            printf("[WOFOST-RESULT] Organ rates: DWRT=%.2f, DWLV=%.2f, DWST=%.2f, DWSO=%.2f\n",
                   rds.DWRT(0), lds.DWLV(0), sds.DWST(0), sods.DWSO(0));
        }
#endif
    }

    // // Finalization (e.g., calculate Harvest Index)
    // void finalize() {
    //     Kokkos::parallel_for("Wofost72_Finalize", s.nCells, KOKKOS_LAMBDA(const int i) {
    //         if (TAGP(i) > 0.0) {
    //             HI(i) = sods.TWSO(i) / TAGP(i);
    //         } else {
    //             HI(i) = -1.0;
    //         }
    //     });
    // }
    //
};

#endif
