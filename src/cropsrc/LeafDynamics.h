/* -*- mode: c++; c-default-style: "linux" -*- */
/*******************************************************************************
 * @file    LeafDynamics.h
 * @brief   WOFOST leaf growth and senescence dynamics (compute kernel)
 *
 * This header implements the three core simulation phases of leaf dynamics
 * in the WOFOST (World Food Studies) crop growth model:
 *
 *   1. initialize()  - Set initial leaf biomass and leaf area index at
 *                       crop emergence from total dry weight (TDWI) and
 *                       partitioning fractions (FR, FL).
 *   2. calc_rates()  - Compute daily rate variables: leaf growth rate,
 *                       leaf death rates (water stress, self-shading,
 *                       frost, ageing), physiological ageing, specific
 *                       leaf area (SLA), and leaf area expansion rates.
 *   3. integrate()   - Advance leaf state by one time step (delt [d]):
 *                       remove dead biomass from oldest leaf age classes,
 *                       shift classes, add new growth, update aggregate
 *                       biomass totals and leaf area index (LAI).
 *
 * Key physical quantities
 * -----------------------
 *   LAI   [-]         Leaf area index (leaf area per unit ground area)
 *   WLV   [kg/ha]     Living leaf dry biomass
 *   DWLV  [kg/ha]     Dead  leaf dry biomass
 *   SLA   [ha/kg]     Specific leaf area (leaf area per unit leaf mass)
 *   GRLV  [kg/(ha d)] Leaf biomass growth rate
 *   DRLV  [kg/(ha d)] Total leaf death rate (max of stress, shading, ageing)
 *
 * Reference: WOFOST 7.1 User Guide, Chapter 4 -- Leaf Dynamics
 ******************************************************************************/

#ifndef _LEAF_DYNAMICS_H_
#define _LEAF_DYNAMICS_H_

#include "../define.h"
#include "CropState.h"               // Static crop parameters (TDWI, SPAN, TBASE, RGRLAI, LAIEM, PERDL)
#include "CropDynamicState.h"        // Dynamic state variables (DVS, ADMI, SAI, PAI)
#include "LeafDynamicsState.h"       // Leaf-specific state arrays (WLV, DWLV, LAI, LV, SLA, LVAGE)
#include "EvapotranspirationState.h" // Transpiration reduction factor RFTRA [-]
#include "MeteoState.h"              // Driving meteorological variables (TMIN, TMAX)
#include "PartitioningState.h"       // Biomass partitioning fractions (FR, FL, FS, FO)
#include "Afgen.h"                   // AFGEN table interpolation utility
#include <cmath>
#include <algorithm>
#include <cstdio>

class LeafDynamics
{

public:
    /*=====================================================================*/
    /*  Initialization                                                     */
    /*=====================================================================*/

    /**
     * @brief Initialize leaf state at crop emergence.
     *
     * Computes the initial living leaf biomass WLV from the total initial
     * crop dry weight (TDWI), root fraction (FR), and leaf fraction (FL).
     * If the resulting LAI falls below the prescribed emergence LAI (LAIEM),
     * WLV is back-calculated so that LAI = LAIEM (WOFOST standard practice).
     *
     * A single leaf age class is created with biomass WLV, SLA from the
     * SLATB table at the current DVS, and age zero.
     *
     * @param[out] lds   Leaf dynamics state arrays (per grid cell)
     * @param[in]  p     Static crop parameters (TDWI, LAIEM, tables)
     * @param[in]  s     Dynamic crop state (DVS)
     * @param[in]  FR    Root partitioning fraction [-] (from PartitioningState)
     * @param[in]  FL    Leaf partitioning fraction [-] (from PartitioningState)
     * @param[in]  SAI   Stem area index [-]
     * @param[in]  PAI   Pod area index [-]
     */
    void initialize(LeafDynamicsState &lds, const CropState &p, const CropDynamicState &s,
                    const realArr &FR, const realArr &FL,
                    const realArr &SAI, const realArr &PAI)
    {

        auto SLATB = p.tables.at("SLATB");   // Specific leaf area table [ha/kg] vs DVS
        auto KDIFTB = p.tables.at("KDIFTB"); // Extinction coefficient table [-] vs DVS

        Kokkos::parallel_for("LeafDynamics_Init", lds.nCells, KOKKOS_LAMBDA(const int i) {
            real TDWI = p.p.TDWI;             // Total initial dry weight [kg/ha]
            real fr = FR(i);                   // Root fraction [-]
            real fl = FL(i);                   // Leaf fraction [-]
            real DVS = s.DVS(i);               // Development stage [-]

#if DEBUG_CROP_GROWTH_MODEL
            if (i == 0)
            {
                printf("[LEAF-INIT-KERNEL] TDWI=%.4f, FL=%.4f, FR=%.4f, DVS=%.4f\n", TDWI, fl, fr, DVS);
            }
#endif
            // --- Initial leaf biomass ---
            // Above-ground fraction = (1 - FR); leaf fraction of above-ground = FL
            real wlv = (TDWI * (1.0 - fr)) * fl;

            // --- WOFOST correction: enforce minimum LAI at emergence ---
            // If the computed LAI from initial biomass and SLA is below the
            // prescribed emergence LAI (LAIEM, typically ~0.07 [-]), the leaf
            // biomass is adjusted upward so that LAI = LAIEM.
            real sla_val = Afgen::lookup(SLATB, DVS);  // SLA at current DVS [ha/kg]
            real calc_lai = wlv * sla_val;              // Computed LAI [-]
            real LAIEM = p.p.LAIEM;                     // Emergence LAI [-], e.g. 0.07

            if (calc_lai < LAIEM)
            {
                // Back-calculate WLV so that LAI equals LAIEM
                wlv = LAIEM / sla_val;
                calc_lai = LAIEM;
            }

            // --- Set scalar state variables ---
            lds.WLV(i) = wlv;               // Living leaf biomass [kg/ha]
            lds.DWLV(i) = 0.0;              // Dead leaf biomass [kg/ha] (no dead leaves at emergence)
            lds.TWLV(i) = wlv;              // Total leaf biomass [kg/ha] (= WLV + DWLV)

            // --- Initialize first leaf age class ---
            lds.num_leaf_classes(i) = 1;
            lds.LV(i, 0)    = wlv;          // Biomass of class 0 [kg/ha]
            lds.SLA(i, 0)   = sla_val;      // Specific leaf area of class 0 [ha/kg]
            lds.LVAGE(i, 0) = 0.0;          // Physiological age of class 0 [d]

            // --- LAI tracking variables ---
            lds.LAIEM(i)  = calc_lai;        // LAI at emergence [-]
            lds.LASUM(i)  = calc_lai;        // Total leaf area from leaf classes [-]
            lds.LAIEXP(i) = calc_lai;        // LAI under exponential growth regime [-]
            lds.LAIMAX(i) = calc_lai;        // Maximum LAI reached so far [-]

            // --- Total green area index = leaf + stem + pod ---
            real sai = SAI(i);               // Stem area index [-]
            real pai = PAI(i);               // Pod area index [-]
            lds.LAI(i) = lds.LASUM(i) + sai + pai;
#if DEBUG_CROP_GROWTH_MODEL
            printf("[LEAF-INIT-KERNEL] Final Init LAI=%.4f (LASUM=%.4f, SAI=%.4f, PAI=%.4f)\n",
                   lds.LAI(i), lds.LASUM(i), sai, pai);
#endif
        });
    }

    /*=====================================================================*/
    /*  Rate Calculation                                                   */
    /*=====================================================================*/

    /**
     * @brief Compute daily rate variables for leaf dynamics.
     *
     * This function evaluates all source and sink terms for leaf biomass
     * and leaf area over one simulation day:
     *
     *  1. Leaf growth rate GRLV = ADMI * FL   [kg/(ha d)]
     *  2. Death by water stress  DSLV1         [kg/(ha d)]
     *  3. Death by self-shading  DSLV2         [kg/(ha d)]
     *  4. Death by frost         DSLV3         [kg/(ha d)]  (placeholder, set to 0)
     *  5. Death by ageing        DALV          [kg/(ha d)]
     *  6. Physiological ageing   FYSAGE        [-]
     *  7. Specific leaf area     SLAT          [ha/kg]
     *  8. Leaf area expansion    GLAIEX, GLASOL [-/d]
     *
     * @param[in,out] lds   Leaf dynamics state (rate arrays are written)
     * @param[in]     p     Static crop parameters
     * @param[in]     s     Dynamic crop state (DVS)
     * @param[in]     m     Meteorological driving variables (TMIN, TMAX)
     * @param[in]     drv_idx  Index into the meteorological time-series arrays
     * @param[in]     ADMI  Above-ground dry matter increase [kg/(ha d)]
     * @param[in]     ets   Evapotranspiration state (RFTRA [-])
     * @param[in]     ps    Partitioning state (FL [-])
     */
    void calc_rates(LeafDynamicsState &lds, const CropState &p, const CropDynamicState &s,
                    const MeteoState &m, int drv_idx,
                    const realArr &ADMI, const EvapotranspirationState &ets,
                    const PartitioningState &ps)
    {

        auto KDIFTB = p.tables.at("KDIFTB"); // Light extinction coefficient [-] vs DVS
        auto SLATB = p.tables.at("SLATB");   // Specific leaf area [ha/kg] vs DVS

        Kokkos::parallel_for("LeafDynamics_Rates", lds.nCells, KOKKOS_LAMBDA(const int i) {
            real DVS = s.DVS(i);               // Development stage [-] (0 = emergence, 1 = flowering, 2 = maturity)
            real TMIN = m.tmin(drv_idx);       // Daily minimum temperature [deg C]
            real TMAX = m.tmax(drv_idx);       // Daily maximum temperature [deg C]
            real TEMP = (TMIN + TMAX) / 2.0;   // Daily mean temperature [deg C]

            /*-------------------------------------------------------------*/
            /*  1. Leaf biomass growth rate                                */
            /*      GRLV = ADMI * FL                                       */
            /*      ADMI = above-ground dry matter increase [kg/(ha d)]    */
            /*      FL   = leaf partitioning fraction [-]                  */
            /*-------------------------------------------------------------*/
            lds.GRLV(i) = ADMI(i) * ps.FL(i);

            /*-------------------------------------------------------------*/
            /*  2. Leaf death rate due to water stress (DSLV1)             */
            /*      DSLV1 = WLV * (1 - RFTRA) * PERDL                     */
            /*      RFTRA = transpiration reduction factor [-] (0..1)      */
            /*      PERDL = maximum relative death rate of leaves          */
            /*               due to water stress [d^-1] (typically 0.03)   */
            /*-------------------------------------------------------------*/
            lds.DSLV1(i) = lds.WLV(i) * (1.0 - ets.RFTRA(i)) * p.p.PERDL;

            /*-------------------------------------------------------------*/
            /*  3. Leaf death rate due to self-shading (DSLV2)             */
            /*      LAICR = critical LAI above which self-shading occurs   */
            /*              = 3.2 / KDIF   [-]                             */
            /*      term = 0.03 * (LAI - LAICR) / LAICR  (clamped 0..0.03)*/
            /*      DSLV2 = WLV * term                                     */
            /*-------------------------------------------------------------*/
            real kdif = Afgen::lookup(KDIFTB, DVS);   // Light extinction coefficient [-]
            real LAICR = 3.2 / kdif;                   // Critical LAI for self-shading [-]
            real term = 0.03 * (lds.LAI(i) - LAICR) / LAICR;
            if (term < 0.0)
                term = 0.0;
            if (term > 0.03)
                term = 0.03;
            lds.DSLV2(i) = lds.WLV(i) * term;

            /*-------------------------------------------------------------*/
            /*  4. Leaf death rate due to frost (DSLV3)                    */
            /*      Currently a placeholder set to 0 [kg/(ha d)]           */
            /*      In a full WOFOST implementation, this would activate   */
            /*      when TMIN drops below a crop-specific frost threshold.  */
            /*-------------------------------------------------------------*/
            lds.DSLV3(i) = 0.0;

            /*-------------------------------------------------------------*/
            /*  5. Maximum stress-related death rate (DSLV)                */
            /*      The dominant stress is the maximum of the three        */
            /*      stress-related death components:                       */
            /*        DSLV = max(DSLV1, DSLV2, DSLV3)                     */
            /*-------------------------------------------------------------*/
            real dslv = fmax(lds.DSLV1(i), fmax(lds.DSLV2(i), lds.DSLV3(i)));
            lds.DSLV(i) = dslv;

            /*-------------------------------------------------------------*/
            /*  6. Leaf death rate due to ageing (DALV)                    */
            /*      Leaves older than SPAN [d] are considered senescent    */
            /*      and their biomass is added to the ageing death rate.   */
            /*      SPAN is typically 28--40 days depending on crop.        */
            /*-------------------------------------------------------------*/
            real dalv_sum = 0.0;
            real SPAN = p.p.SPAN;                       // Life span of leaves [d]
            int n_classes = lds.num_leaf_classes(i);

            for (int c = 0; c < n_classes; ++c)
            {
                if (lds.LVAGE(i, c) > SPAN)
                {
                    dalv_sum += lds.LV(i, c);
                }
            }
            lds.DALV(i) = dalv_sum;

            /*-------------------------------------------------------------*/
            /*  7. Total leaf death rate (DRLV)                            */
            /*      DRLV = max(stress death, ageing death)                 */
            /*      Only the larger of the two mechanisms is applied to    */
            /*      avoid double-counting leaf loss.                       */
            /*-------------------------------------------------------------*/
            lds.DRLV(i) = fmax(dslv, dalv_sum);

            /*-------------------------------------------------------------*/
            /*  8. Physiological ageing factor (FYSAGE)                    */
            /*      FYSAGE = (TEMP - TBASE) / (35 - TBASE)                 */
            /*      Clamped to [0, +inf).                                   */
            /*      TBASE = base temperature for development [deg C]       */
            /*      (typically 0--10 deg C depending on crop)               */
            /*      This factor is accumulated in LVAGE during integrate.  */
            /*-------------------------------------------------------------*/
            real TBASE = p.p.TBASE;
            real fysage = (TEMP - TBASE) / (35.0 - TBASE);
            if (fysage < 0.0)
                fysage = 0.0;
            lds.FYSAGE(i) = fysage;

            /*-------------------------------------------------------------*/
            /*  9. Specific leaf area at current DVS (SLAT)                */
            /*      Interpolated from the SLATB table.  Units: [ha/kg]     */
            /*      Typical range: 0.0015--0.0030 ha/kg                    */
            /*      (i.e. 15--30 m^2/kg)                                    */
            /*-------------------------------------------------------------*/
            lds.SLAT(i) = Afgen::lookup(SLATB, DVS);

            /*-------------------------------------------------------------*/
            /*  10. Leaf area expansion rate                                */
            /*      Two competing expansion pathways are computed and the  */
            /*      minimum determines the actual expansion:               */
            /*                                                               */
            /*      a) Exponential growth pathway (GLAIEX):                 */
            /*         GLAIEX_pot = LAIEXP * RGRLAI * dteff                 */
            /*           LAIEXP  = LAI under exponential growth [-]         */
            /*           RGRLAI  = maximum relative increase in LAI [d^-1]  */
            /*                      (typically 0.0085 for C3 crops)          */
            /*           dteff   = effective temperature [deg C]            */
            /*                                                               */
            /*      b) Source-limited pathway (GLASOL):                     */
            /*         GLASOL = GRLV * SLAT                                 */
            /*           Growth rate [kg/(ha d)] * SLA [ha/kg]              */
            /*                                                               */
            /*      The realized SLA (SLAT) is adjusted to match the        */
            /*      actual expansion: SLAT = GLA / GRLV                     */
            /*                                                               */
            /*      This exponential regime applies only when LAIEXP < 6.0. */
            /*      Above that threshold, the crop relies solely on the     */
            /*      source-limited pathway and GLAIEX = GLASOL = 0.         */
            /*-------------------------------------------------------------*/
            if (lds.LAIEXP(i) < 6.0)
            {
                real dteff = TEMP - TBASE;    // Effective temperature for growth [deg C]
                if (dteff < 0.0)
                    dteff = 0.0;

                real glaiex_pot = lds.LAIEXP(i) * p.p.RGRLAI * dteff;  // Potential exponential LAI increase [-/d]
                lds.GLAIEX(i) = glaiex_pot;    // Store potential exponential rate for integrate()
                lds.GLASOL(i) = lds.GRLV(i) * lds.SLAT(i); // Source-limited LAI increase [-/d]

                real gla = fmin(glaiex_pot, lds.GLASOL(i)); // Actual LAI increase = min of both

                // Adjust effective SLA to match the realized growth
                if (lds.GRLV(i) > 0.0)
                {
                    lds.SLAT(i) = gla / lds.GRLV(i);
                }

#if DEBUG_CROP_GROWTH_MODEL
                if (i == 0)
                {
                    printf("[CPP-LEAF-RATES-EXPANSION] cell=%d, LAIEXP=%.6f, GLAIEX=%.6f, GLASOL=%.6f, GLA=%.6f, SLAT=%.6f\n",
                           i, lds.LAIEXP(i), lds.GLAIEX(i), lds.GLASOL(i), gla, lds.SLAT(i));
                }
#endif
            }
            else
            {
                // Exponential growth phase has ended; LAIEXP >= 6.0 [-]
                lds.GLAIEX(i) = 0.0;
                lds.GLASOL(i) = 0.0;
            }

#if DEBUG_CROP_GROWTH_MODEL
            if (i == 0)
            {
                printf("[CPP-LEAF-RATES-DETAIL] cell=%d, FL=%.6f, ADMI=%.6f, GRLV=%.6f, WLV=%.6f\n",
                       i, ps.FL(i), ADMI(i), lds.GRLV(i), lds.WLV(i));
                printf("[CPP-LEAF-RATES-DETAIL] cell=%d, DVS=%.4f, RFTRA=%.6f\n",
                       i, s.DVS(i), ets.RFTRA(i));
                printf("[CPP-LEAF-RATES-DEATH] cell=%d, DSLV1=%.6f, DSLV2=%.6f, DALV=%.6f, DRLV=%.6f\n",
                       i, lds.DSLV1(i), lds.DSLV2(i), lds.DALV(i), lds.DRLV(i));
                printf("[CPP-LEAF-RATES-AGE] cell=%d, TEMP=%.2f, TBASE=%.2f, FYSAGE=%.6f\n",
                       i, TEMP, TBASE, lds.FYSAGE(i));
            }
#endif
        });
    }

    /*=====================================================================*/
    /*  State Integration                                                  */
    /*=====================================================================*/

    /**
     * @brief Advance leaf state variables by one time step.
     *
     * Integration steps performed in order:
     *   1. Remove dead leaf biomass from oldest age classes first (FIFO).
     *   2. Advance physiological age (LVAGE) of all surviving classes.
     *   3. Shift existing classes to make room for new growth at class 0.
     *   4. Create a new class 0 with the current growth (GRLV * delt).
     *   5. Recompute aggregate totals (WLV, LASUM, DWLV, TWLV, LAI).
     *   6. Update LAIMAX if the new LAI exceeds the previous maximum.
     *   7. Advance exponential LAI tracker (LAIEXP).
     *
     * @param[in,out] lds   Leaf dynamics state arrays
     * @param[in]     s     Dynamic crop state (SAI, PAI for total LAI)
     * @param[in]     delt  Time step size [d] (default 1.0)
     */
    void integrate(LeafDynamicsState &lds, const CropDynamicState &s, real delt = 1.0)
    {

        Kokkos::parallel_for("LeafDynamics_Integrate", lds.nCells, KOKKOS_LAMBDA(const int i) {

#if DEBUG_CROP_GROWTH_MODEL
            if (i == 0)
            {
                printf("[CPP-LEAF-INTEGRATE-START] cell=%d, WLV=%.6f, LASUM=%.6f, LAI=%.6f, GRLV=%.6f, DRLV=%.6f\n",
                       i, lds.WLV(i), lds.LASUM(i), lds.LAI(i), lds.GRLV(i), lds.DRLV(i));
            }
#endif

            real tDRLV = lds.DRLV(i) * delt;  // Total death to apply this time step [kg/ha]
            int n_classes = lds.num_leaf_classes(i);

#if DEBUG_CROP_GROWTH_MODEL
            if (i == 0)
            {
                printf("[CPP-LEAF-INTEGRATE-DEATH] cell=%d, tDRLV=%.6f, n_classes=%d\n",
                       i, tDRLV, n_classes);
            }
#endif

            /*-------------------------------------------------------------*/
            /*  Step 1: Remove dead biomass from oldest classes first       */
            /*      Leaf classes are ordered from newest (0) to oldest.    */
            /*      Dead biomass is deducted starting from the oldest       */
            /*      class and working backwards. Entire classes may be      */
            /*      removed if their biomass is fully consumed.             */
            /*-------------------------------------------------------------*/
            while (n_classes > 0 && tDRLV > 0.0)
            {
                int oldest_idx = n_classes - 1;
                real lv_weight = lds.LV(i, oldest_idx);

                if (tDRLV >= lv_weight)
                {
                    // Entire oldest class dies; remove it
                    tDRLV -= lv_weight;
                    n_classes--;
                }
                else
                {
                    // Partial removal of the oldest class
                    lds.LV(i, oldest_idx) -= tDRLV;
                    tDRLV = 0.0;
                }
            }
            lds.num_leaf_classes(i) = n_classes;

            /*-------------------------------------------------------------*/
            /*  Step 2: Advance physiological age of all surviving classes  */
            /*      LVAGE[c] += FYSAGE * delt                              */
            /*      FYSAGE is the daily physiological ageing increment      */
            /*      computed in calc_rates().                               */
            /*-------------------------------------------------------------*/
            real fysage = lds.FYSAGE(i) * delt;
            for (int c = 0; c < n_classes; ++c)
            {
                lds.LVAGE(i, c) += fysage;
            }

            /*-------------------------------------------------------------*/
            /*  Step 3-4: Shift classes and create new growth class         */
            /*      All existing classes are shifted one position to make   */
            /*      room for the newest growth at index 0.                  */
            /*      Class 0 is initialized with:                            */
            /*        LV[0]    = GRLV * delt   (new leaf biomass [kg/ha])  */
            /*        SLA[0]   = SLAT          (specific leaf area [ha/kg])*/
            /*        LVAGE[0] = 0              (new leaf, age = 0 [d])    */
            /*-------------------------------------------------------------*/
            if (n_classes < LeafDynamicsState::MAX_LEAF_CLASSES)
            {
                for (int c = n_classes; c > 0; --c)
                {
                    lds.LV(i, c)    = lds.LV(i, c - 1);
                    lds.SLA(i, c)   = lds.SLA(i, c - 1);
                    lds.LVAGE(i, c) = lds.LVAGE(i, c - 1);
                }
                lds.LV(i, 0)    = lds.GRLV(i) * delt;   // New leaf biomass [kg/ha]
                lds.SLA(i, 0)   = lds.SLAT(i);           // SLA of new leaves [ha/kg]
                lds.LVAGE(i, 0) = 0.0;                   // Age of new leaves [d]
                lds.num_leaf_classes(i)++;
            }

            /*-------------------------------------------------------------*/
            /*  Step 5: Recompute aggregate totals                          */
            /*      WLV   = sum of LV[c] over all classes [kg/ha]          */
            /*      LASUM = sum of LV[c] * SLA[c] over all classes [-]     */
            /*              (leaf area from biomass × specific leaf area)   */
            /*      DWLV  += DRLV * delt (accumulate dead leaf biomass)    */
            /*      TWLV  = WLV + DWLV (total leaf biomass [kg/ha])        */
            /*      LAI   = LASUM + SAI + PAI (total green area index [-]) */
            /*-------------------------------------------------------------*/
            real sum_lv = 0.0;
            real sum_lasum = 0.0;
            int current_n = lds.num_leaf_classes(i);

            for (int c = 0; c < current_n; ++c)
            {
                sum_lv += lds.LV(i, c);
                sum_lasum += lds.LV(i, c) * lds.SLA(i, c);
            }

            lds.WLV(i) = sum_lv;
            lds.LASUM(i) = sum_lasum;

            lds.DWLV(i) += lds.DRLV(i) * delt;
            lds.TWLV(i) = lds.WLV(i) + lds.DWLV(i);

            real sai = s.SAI(i);                   // Stem area index [-]
            real pai = s.PAI(i);                   // Pod area index [-]
            real old_lai = lds.LAI(i);
            lds.LAI(i) = lds.LASUM(i) + sai + pai;

#if DEBUG_CROP_GROWTH_MODEL
            printf("Debug Integrate: cell=%d, LASUM=%.6f, SAI=%.6f, PAI=%.6f, LAI=%.6f\n",
                   i, lds.LASUM(i), sai, pai, lds.LAI(i));
#endif

            /*-------------------------------------------------------------*/
            /*  Step 6: Update maximum LAI tracker                          */
            /*      LAIMAX records the highest LAI reached during the       */
            /*      growing season, used for diagnostics and output.        */
            /*-------------------------------------------------------------*/
            if (lds.LAI(i) > lds.LAIMAX(i))
            {
                lds.LAIMAX(i) = lds.LAI(i);
            }

            /*-------------------------------------------------------------*/
            /*  Step 7: Advance exponential LAI tracker                     */
            /*      LAIEXP tracks the LAI under the exponential growth      */
            /*      regime. Once LAIEXP >= 6.0, the exponential phase ends  */
            /*      and source-limited growth takes over.                   */
            /*-------------------------------------------------------------*/
            lds.LAIEXP(i) += lds.GLAIEX(i) * delt;

#if DEBUG_CROP_GROWTH_MODEL
            if (i == 0)
            {
                printf("[CPP-LEAF-INTEGRATE-END] cell=%d, WLV=%.6f, LASUM=%.6f, LAI=%.6f (old=%.6f), LAIEXP=%.6f\n",
                       i, lds.WLV(i), lds.LASUM(i), lds.LAI(i), old_lai, lds.LAIEXP(i));
                printf("[CPP-LEAF-INTEGRATE-END] cell=%d, SAI=%.6f, PAI=%.6f, LAIMAX=%.6f\n",
                       i, sai, pai, lds.LAIMAX(i));
            }
#endif
        });
    }
};

#endif
