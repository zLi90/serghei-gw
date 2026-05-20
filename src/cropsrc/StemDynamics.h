/* -*- mode: c++; c-default-style: "linux" -*- */
/*******************************************************************************
 * @file    StemDynamics.h
 * @brief   WOFOST stem growth and senescence dynamics (compute kernel)
 *
 * Implements the three simulation phases for stem biomass in the WOFOST
 * crop growth model:
 *
 *   1. initialize()  - Set initial stem biomass from total dry weight
 *                       (TDWI), root fraction (FR), and stem fraction (FS).
 *   2. calc_rates()  - Compute daily stem growth rate (GRST), stem death
 *                       rate (DRST), net change (GWST), and the effect of
 *                       biomass reallocation from stems to storage organs.
 *   3. integrate()   - Advance stem state by one time step: update living
 *                       and dead biomass, total biomass, and stem area
 *                       index (SAI).
 *
 * Key physical quantities
 * -----------------------
 *   WST   [kg/ha]     Living stem dry biomass
 *   DWST  [kg/ha]     Dead stem dry biomass
 *   TWST  [kg/ha]     Total stem dry biomass (WST + DWST)
 *   SAI   [-]         Stem area index (stem area per unit ground area)
 *   GRST  [kg/(ha d)] Stem biomass growth rate
 *   DRST  [kg/(ha d)] Stem biomass death rate
 *   GWST  [kg/(ha d)] Net stem biomass change rate
 *
 * Reference: WOFOST 7.1 User Guide, Chapter 5 -- Stem Dynamics
 ******************************************************************************/

#ifndef _STEM_DYNAMICS_H_
#define _STEM_DYNAMICS_H_

#include "../define.h"
#include "CropState.h"        // Static crop parameters (TDWI, tables: SSATB, RDRSTB)
#include "CropDynamicState.h" // Dynamic state variables (DVS, ADMI)
#include "StemDynamicsState.h" // Stem-specific state arrays (WST, DWST, TWST, SAI)
#include "Afgen.h"            // AFGEN table interpolation utility
#include <cmath>

class StemDynamics
{

public:
    /*=====================================================================*/
    /*  Initialization                                                     */
    /*=====================================================================*/

    /**
     * @brief Initialize stem state at crop emergence.
     *
     * Computes the initial living stem biomass from total dry weight (TDWI),
     * subtracting the root fraction (FR) and applying the stem fraction (FS)
     * of the above-ground portion:
     *   WST = TDWI * (1 - FR) * FS
     *
     * Dead stem biomass starts at zero.  The initial stem area index (SAI)
     * is computed as WST multiplied by the specific stem area interpolated
     * from the SSATB table at the current development stage (DVS).
     *
     * @param[out] sds  Stem dynamics state arrays (per grid cell)
     * @param[in]  p    Static crop parameters
     * @param[in]  s    Dynamic crop state (DVS)
     * @param[in]  FR   Root partitioning fraction [-]
     * @param[in]  FS   Stem partitioning fraction [-]
     */
    void initialize(StemDynamicsState &sds, const CropState &p, const CropDynamicState &s,
                    const realArr &FR, const realArr &FS)
    {

        auto SSATB = p.tables.at("SSATB");  // Specific stem area table [ha/kg] vs DVS

        Kokkos::parallel_for("StemDynamics_Init", sds.nCells, KOKKOS_LAMBDA(const int i) {

            real TDWI = p.p.TDWI;            // Total initial dry weight [kg/ha]
            real fr = FR(i);                  // Root partitioning fraction [-]
            real fs = FS(i);                  // Stem partitioning fraction [-]

            // --- Initial stem biomass ---
            // Above-ground portion = TDWI * (1 - FR); stem fraction of that = FS
            // WST = (TDWI * (1 - FR)) * FS   [kg/ha]
            sds.WST(i) = (TDWI * (1.0 - fr)) * fs;
            sds.DWST(i) = 0.0;               // No dead stems at emergence [kg/ha]
            sds.TWST(i) = sds.WST(i) + sds.DWST(i);  // Total stem biomass [kg/ha]

            // --- Initial stem area index ---
            // SAI = WST * specific stem area (from SSATB table at current DVS)
            // Units: [kg/ha] * [ha/kg] = [-] (dimensionless area ratio)
            real DVS = s.DVS(i);              // Development stage [-]
            real ssa = Afgen::lookup(SSATB, DVS);  // Specific stem area [ha/kg]
            sds.SAI(i) = sds.WST(i) * ssa;
        });
    }

    /*=====================================================================*/
    /*  Rate Calculation                                                   */
    /*=====================================================================*/

    /**
     * @brief Compute daily rate variables for stem dynamics.
     *
     * Three rate terms are calculated:
     *
     *   1. GRST = ADMI * FS                    [kg/(ha d)]
     *      Stem growth rate from above-ground dry matter increase (ADMI)
     *      and stem partitioning fraction (FS).
     *
     *   2. DRST = WST * RDRSTB(DVS)            [kg/(ha d)]
     *      Stem death rate from the relative death rate table (RDRSTB)
     *      evaluated at the current development stage.
     *
     *   3. GWST = GRST - DRST - REALLOC_ST     [kg/(ha d)]
     *      Net change in living stem biomass, including reallocation
     *      (biomass translocated from stems to storage organs during
     *      grain filling).  REALLOC_ST is positive when stems lose
     *      biomass to storage organs.
     *
     * @param[in,out] sds          Stem dynamics state (rate arrays are written)
     * @param[in]     p            Static crop parameters (RDRSTB table)
     * @param[in]     s            Dynamic crop state (DVS)
     * @param[in]     ADMI         Above-ground dry matter increase [kg/(ha d)]
     * @param[in]     FS           Stem partitioning fraction [-]
     * @param[in]     REALLOC_ST   Biomass reallocated from stems to storage organs [kg/(ha d)]
     */
    void calc_rates(StemDynamicsState &sds, const CropState &p, const CropDynamicState &s,
                    const realArr &ADMI, const realArr &FS, const realArr &REALLOC_ST)
    {

        auto RDRRTB = p.tables.at("RDRSTB");  // Relative death rate of stems [d^-1] vs DVS

        Kokkos::parallel_for("StemDynamics_Rates", sds.nCells, KOKKOS_LAMBDA(const int i) {

            real DVS = s.DVS(i);              // Development stage [-]
            real admi = ADMI(i);              // Above-ground dry matter increase [kg/(ha d)]
            real fs = FS(i);                   // Stem partitioning fraction [-]
            real realloc_st = REALLOC_ST(i);  // Reallocation from stems [kg/(ha d)]

            /*-------------------------------------------------------------*/
            /*  Stem growth rate                                           */
            /*    GRST = ADMI * FS                                         */
            /*    The fraction of above-ground assimilates allocated to     */
            /*    stem growth.                                             */
            /*-------------------------------------------------------------*/
            sds.GRST(i) = admi * fs;

            /*-------------------------------------------------------------*/
            /*  Stem death rate                                            */
            /*    DRST = WST * RDRSTB(DVS)                                 */
            /*    The relative death rate is interpolated from the RDRSTB   */
            /*    table at the current development stage.                   */
            /*-------------------------------------------------------------*/
            real rdr = Afgen::lookup(RDRRTB, DVS);
            sds.DRST(i) = rdr * sds.WST(i);

            /*-------------------------------------------------------------*/
            /*  Net stem biomass change                                    */
            /*    GWST = GRST - DRST - REALLOC_ST                         */
            /*    REALLOC_ST represents translocation of stem reserves     */
            /*    to storage organs (e.g. grain filling from stem          */
            /*    carbohydrates). Positive REALLOC_ST reduces WST.         */
            /*-------------------------------------------------------------*/
            sds.GWST(i) = sds.GRST(i) - sds.DRST(i) - realloc_st;
        });
    }

    /*=====================================================================*/
    /*  State Integration                                                  */
    /*=====================================================================*/

    /**
     * @brief Advance stem state variables by one time step.
     *
     * Updates living stem biomass (WST), accumulates dead stem biomass (DWST),
     * recomputes total stem biomass (TWST), and updates the stem area index
     * (SAI) from the SSATB table at the current DVS.
     *
     * @param[in,out] sds   Stem dynamics state arrays
     * @param[in]     p     Static crop parameters (SSATB table)
     * @param[in]     s     Dynamic crop state (DVS)
     * @param[in]     delt  Time step size [d] (default 1.0)
     */
    void integrate(StemDynamicsState &sds, const CropState &p, const CropDynamicState &s, real delt = 1.0)
    {

        auto SSATB = p.tables.at("SSATB");  // Specific stem area table [ha/kg] vs DVS

        Kokkos::parallel_for("StemDynamics_Integrate", sds.nCells, KOKKOS_LAMBDA(const int i) {

            /*-------------------------------------------------------------*/
            /*  Update stem biomass                                        */
            /*    WST += GWST * delt   (living biomass [kg/ha])            */
            /*    DWST += DRST * delt  (dead biomass accumulation [kg/ha]) */
            /*    TWST = WST + DWST    (total stem biomass [kg/ha])        */
            /*-------------------------------------------------------------*/
            sds.WST(i) += sds.GWST(i) * delt;
            sds.DWST(i) += sds.DRST(i) * delt;
            sds.TWST(i) = sds.WST(i) + sds.DWST(i);

            /*-------------------------------------------------------------*/
            /*  Update stem area index                                     */
            /*    SAI = WST * SSATB(DVS)                                   */
            /*    The specific stem area may change with development stage. */
            /*    Units: [kg/ha] * [ha/kg] = [-]                           */
            /*-------------------------------------------------------------*/
            real DVS = s.DVS(i);
            real ssa = Afgen::lookup(SSATB, DVS);
            sds.SAI(i) = sds.WST(i) * ssa;
        });
    }
};

#endif
