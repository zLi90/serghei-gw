/* -*- mode: c++; c-default-style: "linux" -*- */
/*******************************************************************************
 * @file    StorageOrganDynamics.h
 * @brief   WOFOST storage organ (grain/tuber) growth dynamics (compute kernel)
 *
 * Implements the three simulation phases for storage organ biomass in the
 * WOFOST crop growth model.  Storage organs represent the economically
 * harvested part of the crop: grains in cereals, tubers in potatoes,
 * pods in legumes, etc.
 *
 *   1. initialize()  - Set initial storage organ biomass from total dry
 *                       weight (TDWI), root fraction (FR), and storage
 *                       organ fraction (FO).  Compute initial pod area
 *                       index (PAI).
 *   2. calc_rates()  - Compute daily storage organ growth rate (GRSO),
 *                       death rate (DRSO, always zero in WOFOST), and net
 *                       change including reallocation gains (GWSO).
 *   3. integrate()   - Advance storage organ state by one time step:
 *                       update living and dead biomass, total biomass,
 *                       and pod area index (PAI).
 *
 * Key physical quantities
 * -----------------------
 *   WSO   [kg/ha]     Living storage organ dry biomass (grain/tuber yield)
 *   DWSO  [kg/ha]     Dead storage organ biomass (always 0 in WOFOST)
 *   TWSO  [kg/ha]     Total storage organ biomass (WSO + DWSO)
 *   PAI   [-]         Pod area index (pod area per unit ground area)
 *   GRSO  [kg/(ha d)] Storage organ growth rate
 *   DRSO  [kg/(ha d)] Storage organ death rate (always 0)
 *   GWSO  [kg/(ha d)] Net storage organ biomass change rate
 *
 * Reference: WOFOST 7.1 User Guide, Chapter 7 -- Storage Organ Dynamics
 ******************************************************************************/

#ifndef _STORAGE_ORGAN_DYNAMICS_H_
#define _STORAGE_ORGAN_DYNAMICS_H_

#include "../define.h"
#include "CropState.h"                  // Static crop parameters (TDWI, SPA)
#include "CropDynamicState.h"           // Dynamic state variables (ADMI)
#include "StorageOrganDynamicsState.h"  // Storage organ state arrays (WSO, DWSO, TWSO, PAI)
#include <cmath>

class StorageOrganDynamics
{

public:
    /*=====================================================================*/
    /*  Initialization                                                     */
    /*=====================================================================*/

    /**
     * @brief Initialize storage organ state at crop emergence.
     *
     * Computes the initial living storage organ biomass from total dry
     * weight (TDWI), subtracting the root fraction (FR) and applying the
     * storage organ fraction (FO):
     *   WSO = TDWI * (1 - FR) * FO
     *
     * At emergence, WSO is typically zero or very small because the
     * storage organ fraction (FO) is near zero early in development.
     * The initial pod area index (PAI) is computed as WSO * SPA.
     *
     * @param[out] sods  Storage organ dynamics state arrays (per grid cell)
     * @param[in]  p     Static crop parameters (TDWI, SPA)
     * @param[in]  FR    Root partitioning fraction [-]
     * @param[in]  FO    Storage organ partitioning fraction [-]
     */
    void initialize(StorageOrganDynamicsState &sods, const CropState &p,
                    const realArr &FR, const realArr &FO)
    {

        Kokkos::parallel_for("StorageOrganDynamics_Init", sods.nCells, KOKKOS_LAMBDA(const int i) {

            real TDWI = p.p.TDWI;           // Total initial dry weight [kg/ha]
            real SPA = p.p.SPA;             // Specific pod area [ha/kg]
            real fr = FR(i);                 // Root partitioning fraction [-]
            real fo = FO(i);                 // Storage organ partitioning fraction [-]

            // --- Initial storage organ biomass ---
            // Above-ground portion = TDWI * (1 - FR); storage organ fraction = FO
            // WSO = (TDWI * (1 - FR)) * FO   [kg/ha]
            // Typically near zero at emergence (FO ~ 0 at DVS < 1)
            sods.WSO(i) = (TDWI * (1.0 - fr)) * fo;
            sods.DWSO(i) = 0.0;              // No dead storage organ biomass [kg/ha]
            sods.TWSO(i) = sods.WSO(i) + sods.DWSO(i);  // Total [kg/ha]

            // --- Initial pod area index ---
            // PAI = WSO * SPA   (biomass * specific pod area)
            // Units: [kg/ha] * [ha/kg] = [-] (dimensionless area ratio)
            sods.PAI(i) = sods.WSO(i) * SPA;
        });
    }

    /*=====================================================================*/
    /*  Rate Calculation                                                   */
    /*=====================================================================*/

    /**
     * @brief Compute daily rate variables for storage organ dynamics.
     *
     * Three rate terms are calculated:
     *
     *   1. GRSO = ADMI * FO                [kg/(ha d)]
     *      Storage organ growth rate from above-ground dry matter increase
     *      (ADMI) and storage organ partitioning fraction (FO).  This is
     *      the primary pathway for grain filling.
     *
     *   2. DRSO = 0.0                      [kg/(ha d)]
     *      Storage organ death rate.  Always zero in WOFOST -- storage
     *      organs do not senesce during the simulation.
     *
     *   3. GWSO = GRSO - DRSO + REALLOC_SO [kg/(ha d)]
     *      Net change in storage organ biomass, including gains from
     *      reallocation (biomass translocated from stems and leaves to
     *      storage organs during grain filling).  REALLOC_SO is positive
     *      when storage organs gain biomass from other organs.
     *
     * @param[in,out] sods        Storage organ dynamics state (rate arrays written)
     * @param[in]     p           Static crop parameters
     * @param[in]     ADMI        Above-ground dry matter increase [kg/(ha d)]
     * @param[in]     FO          Storage organ partitioning fraction [-]
     * @param[in]     REALLOC_SO  Biomass reallocated to storage organs [kg/(ha d)]
     */
    void calc_rates(StorageOrganDynamicsState &sods, const CropState &p,
                    const realArr &ADMI, const realArr &FO, const realArr &REALLOC_SO)
    {

        Kokkos::parallel_for("StorageOrganDynamics_Rates", sods.nCells, KOKKOS_LAMBDA(const int i) {

            real admi = ADMI(i);              // Above-ground dry matter increase [kg/(ha d)]
            real fo = FO(i);                   // Storage organ partitioning fraction [-]
            real realloc_so = REALLOC_SO(i);  // Reallocation gain from other organs [kg/(ha d)]

            /*-------------------------------------------------------------*/
            /*  Storage organ growth rate                                   */
            /*    GRSO = ADMI * FO                                          */
            /*    This is the main grain-filling pathway: assimilates        */
            /*    partitioned to storage organs.                            */
            /*-------------------------------------------------------------*/
            sods.GRSO(i) = admi * fo;

            /*-------------------------------------------------------------*/
            /*  Storage organ death rate                                    */
            /*    DRSO = 0.0  (always in WOFOST)                            */
            /*    Storage organs (grains, tubers) are assumed not to die    */
            /*    during the simulation period.                             */
            /*-------------------------------------------------------------*/
            sods.DRSO(i) = 0.0;

            /*-------------------------------------------------------------*/
            /*  Net storage organ biomass change                            */
            /*    GWSO = GRSO - DRSO + REALLOC_SO                          */
            /*    REALLOC_SO is ADDED (positive = gain to storage organs).  */
            /*    This represents translocation of carbohydrates from       */
            /*    stems and leaves during grain filling, a key process      */
            /*    in cereals where stem reserves contribute to yield.       */
            /*-------------------------------------------------------------*/
            sods.GWSO(i) = sods.GRSO(i) - sods.DRSO(i) + realloc_so;
        });
    }

    /*=====================================================================*/
    /*  State Integration                                                  */
    /*=====================================================================*/

    /**
     * @brief Advance storage organ state variables by one time step.
     *
     * Updates living storage organ biomass (WSO), accumulates dead biomass
     * (DWSO, always zero since DRSO = 0), recomputes total biomass (TWSO),
     * and updates the pod area index (PAI).
     *
     * @param[in,out] sods  Storage organ dynamics state arrays
     * @param[in]     p     Static crop parameters (SPA)
     * @param[in]     delt  Time step size [d] (default 1.0)
     */
    void integrate(StorageOrganDynamicsState &sods, const CropState &p, real delt = 1.0)
    {

        Kokkos::parallel_for("StorageOrganDynamics_Integrate", sods.nCells, KOKKOS_LAMBDA(const int i) {

            real SPA = p.p.SPA;               // Specific pod area [ha/kg]

            /*-------------------------------------------------------------*/
            /*  Update storage organ biomass                                */
            /*    WSO += GWSO * delt  (living biomass [kg/ha])              */
            /*    DWSO += DRSO * delt (dead biomass [kg/ha], always 0)      */
            /*    TWSO = WSO + DWSO   (total biomass [kg/ha])               */
            /*-------------------------------------------------------------*/
            sods.WSO(i) += sods.GWSO(i) * delt;
            sods.DWSO(i) += sods.DRSO(i) * delt;
            sods.TWSO(i) = sods.WSO(i) + sods.DWSO(i);

            /*-------------------------------------------------------------*/
            /*  Update pod area index                                       */
            /*    PAI = WSO * SPA                                            */
            /*    The specific pod area (SPA) converts storage organ biomass */
            /*    to green area.  Units: [kg/ha] * [ha/kg] = [-]            */
            /*-------------------------------------------------------------*/
            sods.PAI(i) = sods.WSO(i) * SPA;
        });
    }
};

#endif
