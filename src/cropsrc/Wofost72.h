/**
 * @file Wofost72.h
 * @brief Main container class for the WOFOST 7.2 crop growth simulation model.
 *
 * Wofost72 is the top-level orchestrator for the WOFOST (World Food Studies)
 * crop growth model. It owns instances of all sub-modules and their state
 * containers, and provides two primary operations:
 *
 *   1. initialize()  - Sets initial conditions from crop parameters, meteo
 *                       data, and computed root-zone soil moisture.
 *   2. calc_rates_and_integrate() - Executes one simulation time step:
 *       a. Phenology: compute development rate, update DVS and stage.
 *       b. Assimilation: compute potential gross CO2 assimilation (PGASS).
 *       c. Evapotranspiration: compute water stress (RFTRA) and actual
 *          transpiration/soil evaporation.
 *       d. Partitioning: compute organ biomass fractions (FR, FL, FS, FO).
 *       e. Respiration: compute potential maintenance respiration (PMRES).
 *       f. Carbon balance: apply water stress to assimilation, subtract
 *          respiration, compute dry matter increase (DMI).
 *       g. Organ dynamics: distribute DMI to roots, stems, leaves, storage.
 *       h. Integration: advance all state variables by one time step.
 *       i. Accumulate totals: TAGP, GASST, MREST, CTRAT, CEVST, HI.
 *
 * All per-cell state is stored in Kokkos Views for GPU portability. Local
 * variable captures (auto local_X = s.X) are used inside kernels to avoid
 * capturing the 'this' pointer, which is incompatible with Kokkos lambda
 * capture on the device.
 *
 * Key accumulated outputs:
 *   - TAGP  : Total Above-Ground Production [kg/ha]
 *   - GASST : Total Gross Assimilation [kg CH2O/ha]
 *   - MREST : Total Maintenance Respiration [kg CH2O/ha]
 *   - CTRAT : Cumulative Crop Transpiration [cm]
 *   - CEVST : Cumulative Soil Evaporation [cm]
 *   - HI    : Harvest Index [-] (storage organ DM / total above-ground DM)
 *
 * Reference: Boogaard, H.L., van Diepen, C.A., Rotter, R.P., Cabrera, J.C.M.C.A.,
 *   van Laar, H.H. (1998). WOFOST 7.1 User Guide. Technical Document 52,
 *   DLO Winand Staring Centre, Wageningen.
 *
 * @see CropDynamicState.h   for the main dynamic state structure.
 * @see Phenology.h          for phenological development.
 * @see Assimilation.h       for CO2 assimilation.
 * @see Respiration.h        for maintenance respiration.
 * @see Evapotranspiration.h for evapotranspiration and water stress.
 * @see Partitioning.h       for biomass partitioning.
 * @see RootDynamics.h       for root growth.
 * @see StemDynamics.h       for stem growth.
 * @see StorageOrganDynamics.h for storage organ growth.
 * @see LeafDynamics.h       for leaf growth and senescence.
 */

/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _WOFOST72_H_
#define _WOFOST72_H_

#include "../define.h"

/* --- Static parameters and dynamic state structures --- */
#include "CropState.h"        /* Static crop parameters */
#include "CropDynamicState.h" /* Dynamic crop states (DVS, LAI, etc.) */
#include "MeteoState.h"       /* Weather driving data */

/* --- Sub-module classes --- */
#include "Phenology.h"                 /* Phenological development */
#include "Assimilation.h"              /* CO2 assimilation */
#include "AssimilateState.h"           /* Assimilation output states */
#include "Respiration.h"               /* Maintenance respiration */
#include "RespirationState.h"          /* Respiration output states */
#include "Partitioning.h"              /* Biomass partitioning fractions */
#include "PartitioningState.h"         /* Partitioning output states */
#include "RootDynamics.h"              /* Root growth and depth */
#include "RootDynamicsState.h"         /* Root dynamics output states */
#include "StemDynamics.h"              /* Stem growth */
#include "StemDynamicsState.h"         /* Stem dynamics output states */
#include "StorageOrganDynamics.h"      /* Storage organ (grain/tuber) growth */
#include "StorageOrganDynamicsState.h" /* Storage organ dynamics output states */
#include "LeafDynamics.h"              /* Leaf growth and senescence */
#include "LeafDynamicsState.h"         /* Leaf dynamics output states */
#include "Evapotranspiration.h"        /* Evapotranspiration and water stress */
#include "EvapotranspirationState.h"   /* Evapotranspiration output states */

#include <cmath>
#include <algorithm>

class Wofost72
{

public:
    /* ================================================================== */
    /* Sub-module instances                                                */
    /* ================================================================== */

    Phenology pheno;                /**< Phenological development module. */
    Assimilation assim;             /**< CO2 assimilation module. */
    Respiration mres;               /**< Maintenance respiration module. */
    Evapotranspiration evtra;       /**< Evapotranspiration and water stress module. */
    Partitioning part;              /**< Biomass partitioning module. */
    RootDynamics ro_dynamics;       /**< Root growth dynamics module. */
    StemDynamics st_dynamics;       /**< Stem growth dynamics module. */
    StorageOrganDynamics so_dynamics; /**< Storage organ dynamics module. */
    LeafDynamics lv_dynamics;       /**< Leaf growth and senescence module. */

    /* ================================================================== */
    /* Dynamic state containers (per-cell Kokkos Views)                    */
    /* ================================================================== */

    CropDynamicState s;              /**< General crop states: DVS, TSUM, LAI, RD,
                                          GASS, MRES, DMI, etc. */
    AssimilateState as;              /**< Assimilation outputs: PGASS, etc. */
    RespirationState rs;             /**< Respiration outputs: PMRES, etc. */
    PartitioningState ps;            /**< Partitioning fractions: FR, FL, FS, FO. */
    RootDynamicsState rds;           /**< Root states: WRT, TWRT, RD. */
    StemDynamicsState sds;           /**< Stem states: WST, TWST, SAI. */
    StorageOrganDynamicsState sods;  /**< Storage organ states: WSO, TWSO, PAI. */
    LeafDynamicsState lds;           /**< Leaf states: WLV, TWLV, LAI, LASUM. */
    EvapotranspirationState ets;     /**< ET outputs: RFTRA, TRA, TRAMX, EVS. */
    MeteoState meteo;                /**< Local copy of meteorological state. */

    /* ================================================================== */
    /* Main WOFOST Accumulated Output Variables                            */
    /* ================================================================== */

    realArr TAGP;   /**< Total Above-Ground Production [kg/ha].
                         TAGP = TWLV + TWST + TWSO (leaves + stems + storage).
                         Accumulated over the growing season. */
    realArr GASST;  /**< Total (cumulative) gross assimilation [kg CH2O/ha].
                         Sum of daily GASS values over the season. */
    realArr MREST;  /**< Total (cumulative) maintenance respiration [kg CH2O/ha].
                         Sum of daily MRES values over the season. */
    realArr CTRAT;  /**< Cumulative crop transpiration [cm].
                         Sum of daily TRA values. */
    realArr CEVST;  /**< Cumulative soil evaporation [cm].
                         Sum of daily EVS values. */
    realArr HI;     /**< Harvest Index [-].
                         HI = TWSO / TAGP (ratio of storage organ DM
                         to total above-ground DM). Range: 0.0 -- 1.0. */

    /* ================================================================== */
    /* Crop completion control                                             */
    /* ================================================================== */

    intArr  CROP_FINISHED_FLAG; /**< Per-cell flag: 0 = crop still growing,
                                      1 = crop reached maturity. */
    realArr CROP_FINISH_DVS;    /**< DVS value at which the crop finished [-].
                                      Typically equal to DVSEND. */

    /* ================================================================== */
    /* Default Constructor                                                 */
    /* ================================================================== */

    /**
     * @brief Default constructor. Does not allocate; call allocate() separately.
     */
    Wofost72()
    {
        /* Intentionally empty: allocation is deferred to allocate(). */
    }

    /* ================================================================== */
    /* Memory Allocation                                                   */
    /* ================================================================== */

    /**
     * @brief Allocate all Kokkos View arrays for the given number of cells.
     *
     * Must be called once before initialize(). Delegates to sub-module
     * state allocate() methods and creates the main WOFOST output arrays.
     *
     * @param nCells  Number of surface grid cells (one crop simulation per cell).
     */
    void allocate(int nCells)
    {
        /* Allocate all sub-module state arrays */
        s.allocate(nCells);
        as.allocate(nCells);
        rs.allocate(nCells);
        ps.allocate(nCells);
        rds.allocate(nCells);
        sds.allocate(nCells);
        sods.allocate(nCells);
        lds.allocate(nCells);
        ets.allocate(nCells);

        /* Allocate main WOFOST accumulated output arrays */
        TAGP  = realArr("TAGP", nCells);
        GASST = realArr("GASST", nCells);
        MREST = realArr("MREST", nCells);
        CTRAT = realArr("CTRAT", nCells);
        CEVST = realArr("CEVST", nCells);
        HI    = realArr("HI", nCells);

        /* Allocate crop completion control arrays */
        CROP_FINISHED_FLAG = intArr("CROP_FINISHED_FLAG", nCells);
        CROP_FINISH_DVS    = realArr("CROP_FINISH_DVS", nCells);
    }

    /* ================================================================== */
    /* Initialization                                                      */
    /* ================================================================== */

    /**
     * @brief Initialize all WOFOST states from crop parameters and meteo data.
     *
     * Performs the following initialization sequence:
     *   1. Initialize phenology (DVS, TSUM, STAGE).
     *   2. Initialize partitioning (FR, FL, FS, FO from initial DVS).
     *   3. Initialize organ dynamics (biomass distributed from TDWI, LAI
     *      from LAIEM, rooting depth from RDI).
     *   4. Synchronize cross-module states (LAI, SAI, PAI, RD) from
     *      sub-modules into the shared CropDynamicState.
     *   5. Zero all cumulative accumulators (TAGP, GASST, MREST, etc.).
     *   6. Verify initial biomass partitioning (TDWI = TAGP + TWRT).
     *
     * @param p                      CropState with static parameters.
     * @param m                      MeteoState with weather data.
     * @param initial_root_zone_SM   Pre-computed root-zone soil moisture [cm3/cm3].
     * @param current_day_of_year    Start day-of-year (1--365/366).
     */
    void initialize(const CropState &p, const MeteoState &m, const realArr &initial_root_zone_SM, int current_day_of_year)
    {
        /* --- 1. Initialize Phenology: set DVS, TSUM, STAGE from crop params --- */
        pheno.initialize(s, p);

        /* --- 2. Initialize Partitioning: compute FR, FL, FS, FO at initial DVS --- */
        part.initialize(ps, p, s);

#if DEBUG_CROP_GROWTH_MODEL
        /*
         * Debug: print initial partitioning fractions for each cell.
         * Local variable captures are used to avoid capturing 'this'
         * pointer inside Kokkos lambda on the device.
         */
        auto local_debug_FR = ps.FR;
        auto local_debug_FL = ps.FL;
        auto local_debug_FS = ps.FS;
        auto local_debug_FO = ps.FO;
        Kokkos::parallel_for("Wofost72_SyncPartitioning", s.nCells, KOKKOS_LAMBDA(const int i) {
            printf("[WOFOST-INIT] cell=%d, Partitioning - FR=%.4f, FL=%.4f, FS=%.4f, FO=%.4f\n",
                   i, local_debug_FR(i), local_debug_FL(i), local_debug_FS(i), local_debug_FO(i));
        });
        Kokkos::fence(); /* Ensure kernel completes before proceeding */
#endif

        /* --- 3. Initialize Organ Dynamics --- */

        /*
         * Roots: initial biomass = TDWI * FR, rooting depth = RDI.
         * Stems: initial biomass = TDWI * (1-FR) * FS.
         * Storage organs: initial biomass = TDWI * (1-FR) * FO.
         * Leaves: initial LAI = LAIEM, biomass = TDWI * (1-FR) * FL.
         */
        ro_dynamics.initialize(rds, p, s);
        ro_dynamics.set_initial_biomass(rds, p, ps.FR);

        st_dynamics.initialize(sds, p, s, ps.FR, ps.FS);
        so_dynamics.initialize(sods, p, ps.FR, ps.FO);
        lv_dynamics.initialize(lds, p, s, ps.FR, ps.FL, sds.SAI, sods.PAI);

        /* ================================================================== */
        /* Synchronize cross-module states from sub-modules into shared state  */
        /*                                                                    */
        /* Because Kokkos lambda capture does not allow capturing 'this',     */
        /* all Kokkos Views used inside the parallel_for must be captured     */
        /* as local variables by value. This is a standard Kokkos pattern.    */
        /* ================================================================== */

        /* Shared area index and rooting depth variables */
        auto local_s_LAI = s.LAI;
        auto local_s_SAI = s.SAI;
        auto local_s_PAI = s.PAI;
        auto local_s_RD  = s.RD;

        /* Sub-module area index and rooting depth sources */
        auto local_lds_LAI   = lds.LAI;
        auto local_sds_SAI   = sds.SAI;
        auto local_sods_PAI  = sods.PAI;
        auto local_rds_RD    = rds.RD;

        /* Total biomass from sub-modules (for TAGP computation) */
        auto local_lds_TWLV  = lds.TWLV;
        auto local_sds_TWST  = sds.TWST;
        auto local_sods_TWSO = sods.TWSO;
        auto local_rds_TWRT  = rds.TWRT;

        /* Main WOFOST accumulated output arrays */
        auto local_TAGP              = TAGP;
        auto local_GASST             = GASST;
        auto local_MREST             = MREST;
        auto local_CTRAT             = CTRAT;
        auto local_CEVST             = CEVST;
        auto local_HI                = HI;
        auto local_CROP_FINISHED_FLAG = CROP_FINISHED_FLAG;
        auto local_CROP_FINISH_DVS   = CROP_FINISH_DVS;

        /* Extract scalar parameter for use in device kernel */
        real p_TDWI = p.p.TDWI;

        /* --- 5. Synchronize states and initialize accumulators --- */
        Kokkos::parallel_for("Wofost72_MainInit", s.nCells, KOKKOS_LAMBDA(const int i) {

            /* Copy area indices and rooting depth from sub-modules to shared state */
            local_s_LAI(i) = local_lds_LAI(i);
            local_s_SAI(i) = local_sds_SAI(i);
            local_s_PAI(i) = local_sods_PAI(i);
            local_s_RD(i)  = local_rds_RD(i);

            /* Initial total above-ground biomass [kg/ha] */
            local_TAGP(i) = local_lds_TWLV(i) + local_sds_TWST(i) + local_sods_TWSO(i);

            /* Zero all cumulative accumulators */
            local_GASST(i) = 0.0;
            local_MREST(i) = 0.0;
            local_CTRAT(i) = 0.0;
            local_CEVST(i) = 0.0;
            local_HI(i)    = 0.0;

            /* Crop is not yet finished */
            local_CROP_FINISHED_FLAG(i) = 0;
            local_CROP_FINISH_DVS(i)    = 0.0;

            /*
             * Verify initial biomass partitioning:
             * TDWI should equal TAGP + TWRT (above-ground + root biomass).
             * A non-zero checksum indicates an error in partitioning logic.
             */
            real checksum = p_TDWI - local_TAGP(i) - local_rds_TWRT(i);
            if (Kokkos::fabs(checksum) > 0.0001)
            {
                printf("Error in partitioning of initial biomass (TDWI) at cell %d! Checksum: %f\n", i, checksum);
            }
        });
        Kokkos::fence();
    }

    /* ================================================================== */
    /* Main Simulation Step                                                */
    /* ================================================================== */

    /**
     * @brief Execute one daily time step: compute rates and integrate all states.
     *
     * This is the core simulation loop for WOFOST. The execution order is:
     *
     *   Phase 1 - Rate Calculation:
     *     1. Phenology:         development rate (DVR), temperature sum (DTSUM).
     *     2. Emergence check:   skip growth modules if crop has not emerged.
     *     3. Assimilation:      potential gross assimilation (PGASS).
     *     4. Evapotranspiration: water stress factor (RFTRA), transpiration.
     *     5. Partitioning:      organ fractions (FR, FL, FS, FO) at current DVS.
     *     6. Respiration:       potential maintenance respiration (PMRES).
     *     7. Carbon balance:    GASS = PGASS * RFTRA; MRES = min(GASS, PMRES);
     *                           ASRC = GASS - MRES; DMI = CVF * ASRC.
     *
     *   Phase 2 - Organ Rate Calculation:
     *     Root, stem, storage, and leaf growth rates from DMI partitioning.
     *
     *   Phase 3 - Integration:
     *     Advance all state variables by delt_day.
     *
     *   Phase 4 - Accumulation:
     *     Update cumulative totals (TAGP, GASST, MREST, CTRAT, CEVST).
     *
     * @param p                    CropState with static parameters.
     * @param m                    MeteoState with weather driving data.
     * @param current_time_s       Current simulation time [s] (for meteo indexing).
     * @param root_zone_SM         Current root-zone soil moisture [cm3/cm3].
     * @param current_day_of_year  Current day-of-year (1--365/366).
     * @param delt_day             Time step length [d]. Default: 1.0 (daily).
     */
    void calc_rates_and_integrate(const CropState &p, const MeteoState &m,
                                  real current_time_s, const realArr &root_zone_SM, int current_day_of_year, real delt_day = 1.0)
    {
        /* Determine the meteorological record index for this time step */
        int drv_idx = m.get_index_at_time(current_time_s);

#if DEBUG_CROP_GROWTH_MODEL
        if (s.nCells > 0)
        {
            printf("[WOFOST] Step: DOY=%d, Time=%.1f, drv_idx=%d\n", current_day_of_year, current_time_s, drv_idx);
            printf("[WOFOST] Meteo: TMIN=%.1f, TMAX=%.1f, IRRAD=%.1f, ET0=%.1f\n",
                   m.tmin(drv_idx), m.tmax(drv_idx), m.irrad(drv_idx), m.et0(drv_idx));
            printf("[WOFOST] Before calc: LAI=%.2f, DVS=%.2f\n", s.LAI(0), s.DVS(0));
        }
#endif

        /* ============================================================== */
        /* Phase 1: Rate Calculation                                      */
        /* ============================================================== */

        /* --- 1. Phenology: compute development rate and update DVS/TSUM --- */
        pheno.calc_rates(s, p, m, drv_idx, current_day_of_year);

#if DEBUG_CROP_GROWTH_MODEL
        if (s.nCells > 0)
        {
            printf("[WOFOST-PHENO] After calc_rates: STAGE=%d, DVS=%.2f, TSUM=%.2f\n",
                   (int)s.STAGE(0), s.DVS(0), s.TSUM(0));
        }
#endif

        /* --- 1b. Emergence check --- */
        /*
         * If the crop is still in the EMERGING stage (pre-emergence), set a
         * per-cell flag to skip all growth-related calculations. Only
         * phenology (temperature accumulation for emergence) is active.
         */
        auto local_STAGE = s.STAGE;
        auto local_skip  = s.skip_other_modules;

        Kokkos::parallel_for("Wofost72_CheckEmergence", s.nCells, KOKKOS_LAMBDA(const int i) {
            if (local_STAGE(i) == CropStage::EMERGING) {
                local_skip(i) = 1;  /* Skip growth modules for this cell */
            } else {
                local_skip(i) = 0;  /* Normal operation */
            }
        });
        Kokkos::fence();

#if DEBUG_CROP_GROWTH_MODEL
        if (s.nCells > 0)
        {
            printf("[WOFOST-EMERG] skip_other_modules=%d\n", s.skip_other_modules(0));
        }
#endif

        /* --- 2. Assimilation: potential gross CO2 assimilation [kg CH2O/ha/h] --- */
        assim.calc_rates(as, p, s, m, drv_idx, current_day_of_year);

        /* --- 3. Evapotranspiration: water stress factor and actual fluxes --- */
        evtra.calc_rates(ets, p, s, m, drv_idx, root_zone_SM);

#if DEBUG_CROP_GROWTH_MODEL
        if (s.nCells > 0)
        {
            printf("[EVAP-DEBUG] RFTRA=%.6f, TRA=%.6f, TRAMX=%.6f\n",
                   ets.RFTRA(0), ets.TRA(0), ets.TRAMX(0));
            printf("[ASSIM-DEBUG] PGASS=%.6f\n", as.PGASS(0));
        }
#endif

        /* --- 4. Partitioning: organ biomass fractions at current DVS --- */
        part.update(ps, p, s);

        /* --- 5. Respiration: potential maintenance respiration --- */
        mres.calc_rates(rs, p, s, m, drv_idx, lds, sds, rds, sods);

        /* ============================================================== */
        /* Main Carbon Balance Loop                                       */
        /*                                                                */
        /* Computes actual assimilation (GASS), maintenance respiration  */
        /* (MRES), net available assimilates (ASRC), dry matter increase */
        /* (DMI), and above-ground DMI (ADMI).                            */
        /*                                                                */
        /* Local variable captures are required for Kokkos lambda device */
        /* execution (cannot capture 'this' pointer).                     */
        /* ============================================================== */

        /* Capture all Views needed inside the kernel */
        auto local_GASS      = s.GASS;
        auto local_PGASS     = as.PGASS;
        auto local_RFTRA     = ets.RFTRA;
        auto local_MRES      = s.MRES;
        auto local_PMRES     = rs.PMRES;
        auto local_ASRC      = s.ASRC;
        auto local_FR        = ps.FR;
        auto local_FL        = ps.FL;
        auto local_FS        = ps.FS;
        auto local_FO        = ps.FO;
        auto local_DMI       = s.DMI;
        auto local_ADMI      = s.ADMI;
        auto local_REALLOC_LV = s.REALLOC_LV;
        auto local_REALLOC_ST = s.REALLOC_ST;
        auto local_REALLOC_SO = s.REALLOC_SO;

        /* Capture scalar conversion efficiencies for device use */
        real p_CVL = p.p.CVL;
        real p_CVS = p.p.CVS;
        real p_CVO = p.p.CVO;
        real p_CVR = p.p.CVR;

        Kokkos::parallel_for("Wofost72_CarbonLoop", s.nCells, KOKKOS_LAMBDA(const int i) {
            /* Skip non-emerged cells */
            if (local_skip(i) == 1) return;

            /* --- Actual gross assimilation (water-stress corrected) --- */
            local_GASS(i) = local_PGASS(i) * local_RFTRA(i);

            /* --- Actual maintenance respiration (cannot exceed GASS) --- */
            local_MRES(i) = fmin(local_GASS(i), local_PMRES(i));

            /* --- Net available assimilates --- */
            local_ASRC(i) = local_GASS(i) - local_MRES(i);

            /* --- Weighted conversion efficiency (CVF) --- */
            /*
             * CVF converts CH2O assimilates to dry matter, accounting for
             * different conversion costs per organ:
             *   CVF = 1 / ( FL/CVL*(1-FR) + FS/CVS*(1-FR) + FO/CVO*(1-FR) + FR/CVR )
             */
            real FR = local_FR(i);
            real FL = local_FL(i);
            real FS = local_FS(i);
            real FO = local_FO(i);

            real den = (FL / p_CVL + FS / p_CVS + FO / p_CVO) * (1.0 - FR) + FR / p_CVR;
            real CVF = (den != 0.0) ? (1.0 / den) : 0.0;

            /* --- Total dry matter increase [kg DM/ha/d] --- */
            local_DMI(i) = CVF * local_ASRC(i);

            /*
             * Carbon balance verification checksum:
             * GASS - MRES - (FR + (FL+FS+FO)*(1-FR)) * DMI/CVF should equal 0.
             * Normalised by GASS to give a relative error.
             */
            real py_checksum_den = fmax(0.0001, local_GASS(i));
            real py_checksum = (local_GASS(i) - local_MRES(i) - (FR + (FL + FS + FO) * (1.0 - FR)) * local_DMI(i) / CVF) * 1.0 / py_checksum_den;
            if (Kokkos::fabs(py_checksum) >= 0.0001)
            {
                printf("Carbon balance error at cell %d, checksum: %f\n", i, py_checksum);
            }

            /* --- Reallocation fluxes (currently inactive, set to 0) --- */
            local_REALLOC_LV(i) = 0.0;
            local_REALLOC_ST(i) = 0.0;
            local_REALLOC_SO(i) = 0.0;

            /* --- Above-ground dry matter increase [kg DM/ha/d] --- */
            local_ADMI(i) = (1.0 - FR) * local_DMI(i);

#if DEBUG_CROP_GROWTH_MODEL
            if (i == 0)
            {
                printf("[CARBON-LOOP] cell=%d, PGASS=%.6f, RFTRA=%.4f, GASS=%.6f, PMRES=%.6f, MRES=%.6f, ASRC=%.6f\n",
                       i, local_PGASS(i), local_RFTRA(i), local_GASS(i), local_PMRES(i), local_MRES(i), local_ASRC(i));
                printf("[CARBON-LOOP] cell=%d, FR=%.4f, FL=%.4f, FS=%.4f, FO=%.4f, CVF=%.6f, DMI=%.6f, ADMI=%.6f\n",
                       i, FR, FL, FS, FO, CVF, local_DMI(i), local_ADMI(i));
            }
#endif
        });
        Kokkos::fence();

        /* ============================================================== */
        /* Phase 2: Organ Rate Calculation                                */
        /*                                                                */
        /* Each organ module computes its growth rate from the partitioned */
        /* DMI and applies senescence/mortality losses.                    */
        /* ============================================================== */

        ro_dynamics.calc_rates(rds, p, s, s.DMI, ps.FR);        /* Roots: growth from FR * DMI */
        st_dynamics.calc_rates(sds, p, s, s.ADMI, ps.FS, s.REALLOC_ST);  /* Stems: growth from FS * ADMI */
        so_dynamics.calc_rates(sods, p, s.ADMI, ps.FO, s.REALLOC_SO);    /* Storage: growth from FO * ADMI */
        lv_dynamics.calc_rates(lds, p, s, m, drv_idx, s.ADMI, ets, ps);  /* Leaves: growth from FL * ADMI + senescence */

        /* ============================================================== */
        /* Phase 3: Integration                                           */
        /*                                                                */
        /* Advance all state variables by the time step (delt_day).       */
        /* ============================================================== */

        /* Phenology: advance DVS, TSUM, STAGE */
        pheno.integrate(s, p, delt_day);

        /* Partitioning: re-compute fractions at new DVS */
        part.update(ps, p, s);

        /* Organ states: integrate biomass, LAI, rooting depth */
        ro_dynamics.integrate(rds, delt_day);
        st_dynamics.integrate(sds, p, s, delt_day);
        so_dynamics.integrate(sods, p, delt_day);
        lv_dynamics.integrate(lds, s, delt_day);

        /* ============================================================== */
        /* Phase 4: Accumulation and State Synchronization                */
        /*                                                                */
        /* After integration, synchronize the shared cross-module states  */
        /* (LAI, SAI, PAI, RD) from organ modules, update cumulative     */
        /* totals (TAGP, GASST, MREST, CTRAT, CEVST), and detect crop    */
        /* maturity.                                                      */
        /* ============================================================== */

        /* Capture Views for the final synchronization kernel */
        auto local_s_LAI     = s.LAI;
        auto local_s_SAI     = s.SAI;
        auto local_s_PAI     = s.PAI;
        auto local_s_RD      = s.RD;

        auto local_lds_LAI   = lds.LAI;
        auto local_sds_SAI   = sds.SAI;
        auto local_sods_PAI  = sods.PAI;
        auto local_rds_RD    = rds.RD;

        auto local_TAGP      = TAGP;
        auto local_lds_TWLV  = lds.TWLV;
        auto local_sds_TWST  = sds.TWST;
        auto local_sods_TWSO = sods.TWSO;

        auto local_GASST     = GASST;
        auto local_s_GASS    = s.GASS;
        auto local_MREST     = MREST;
        auto local_s_MRES    = s.MRES;

        auto local_CTRAT     = CTRAT;
        auto local_ets_TRA   = ets.TRA;
        auto local_CEVST     = CEVST;
        auto local_ets_EVS   = ets.EVS;

        auto local_CROP_FINISHED_FLAG = CROP_FINISHED_FLAG;
        auto local_CROP_FINISH_DVS    = CROP_FINISH_DVS;
        auto local_s_DVS              = s.DVS;

        real local_delt_day = delt_day;

        Kokkos::parallel_for("Wofost72_MainIntegrate", s.nCells, KOKKOS_LAMBDA(const int i) {
            /* Skip non-emerged cells (they only accumulate phenology for emergence) */
            if (local_skip(i) == 1 && local_STAGE(i) == CropStage::EMERGING) {
                return;
            }

            /* --- Synchronize area indices and rooting depth from sub-modules --- */
            local_s_LAI(i) = local_lds_LAI(i);
            local_s_SAI(i) = local_sds_SAI(i);
            local_s_PAI(i) = local_sods_PAI(i);
            local_s_RD(i)  = local_rds_RD(i);

            /* --- Update total above-ground production [kg/ha] --- */
            local_TAGP(i) = local_lds_TWLV(i) + local_sds_TWST(i) + local_sods_TWSO(i);

            /* --- Accumulate daily carbon fluxes into seasonal totals --- */
            local_GASST(i) += local_s_GASS(i) * local_delt_day;  /* [kg CH2O/ha] */
            local_MREST(i) += local_s_MRES(i) * local_delt_day;  /* [kg CH2O/ha] */

            /* --- Accumulate daily water fluxes into seasonal totals --- */
            local_CTRAT(i) += local_ets_TRA(i) * local_delt_day;  /* [cm] */
            local_CEVST(i) += local_ets_EVS(i) * local_delt_day;  /* [cm] */

            /* --- Detect crop maturity --- */
            if (local_STAGE(i) == CropStage::MATURE && local_CROP_FINISHED_FLAG(i) == 0) {
                local_CROP_FINISHED_FLAG(i) = 1;
                local_CROP_FINISH_DVS(i) = local_s_DVS(i);
            }
        });
        Kokkos::fence();

#if DEBUG_CROP_GROWTH_MODEL
        if (s.nCells > 0)
        {
            printf("[WOFOST-RESULT] After integrate: LAI=%.2f, DVS=%.2f\n",
                   s.LAI(0), s.DVS(0));
            printf("[WOFOST-RESULT] Rates: GASS=%.2f, MRES=%.2f, DMI=%.2f, ADMI=%.2f\n",
                   s.GASS(0), s.MRES(0), s.DMI(0), s.ADMI(0));
        }
#endif
    }
};

#endif /* _WOFOST72_H_ */
