/* -*- mode: c++; c-default-style: "linux" -*- */
/*******************************************************************************
 * @file    RootDynamics.h
 * @brief   WOFOST root growth and rooting depth dynamics (compute kernel)
 *
 * Implements the three simulation phases for root biomass and rooting depth
 * in the WOFOST crop growth model:
 *
 *   1. initialize()       - Set initial root depth (RD) from crop and soil
 *                            parameters, and initial root biomass from the
 *                            root partitioning fraction (FR).
 *   2. set_initial_biomass() - Separate helper to set WRT after partitioning
 *                            has been initialised (called after Partitioning).
 *   3. calc_rates()       - Compute daily root biomass growth rate (GRRT),
 *                            root death rate (DRRT), net biomass change (GWRT),
 *                            and root depth growth rate (RR).
 *   4. integrate()        - Advance root state by one time step.
 *
 * Key physical quantities
 * -----------------------
 *   RD    [cm]        Current rooting depth (depth of the root zone)
 *   RDM   [cm]        Maximum attainable rooting depth
 *   WRT   [kg/ha]     Living root dry biomass
 *   DWRT  [kg/ha]     Dead root dry biomass
 *   TWRT  [kg/ha]     Total root dry biomass (WRT + DWRT)
 *   RR    [cm/d]      Rate of root depth increase
 *   GRRT  [kg/(ha d)] Root biomass growth rate
 *   DRRT  [kg/(ha d)] Root biomass death rate
 *   GWRT  [kg/(ha d)] Net root biomass change rate
 *
 * Reference: WOFOST 7.1 User Guide, Chapter 6 -- Root Dynamics
 ******************************************************************************/

#ifndef _ROOT_DYNAMICS_H_
#define _ROOT_DYNAMICS_H_

#include "../define.h"
#include "CropState.h"         // Static crop parameters (TDWI, RDI, RDMCR, RRI, tables: RDRRTB)
#include "CropDynamicState.h"  // Dynamic state variables (DVS)
#include "RootDynamicsState.h" // Root-specific state arrays (RDM, RD, WRT, DWRT, TWRT, rates)
#include "Afgen.h"             // AFGEN table interpolation utility
#include <cmath>
#include <algorithm> // for std::min, std::max

class RootDynamics
{

public:
    /*=====================================================================*/
    /*  Initialization                                                     */
    /*=====================================================================*/

    /**
     * @brief Initialize root state at crop emergence.
     *
     * Sets the maximum attainable rooting depth (RDM) as the minimum of
     * the crop-specific maximum root depth (RDMCR) and the soil-limited
     * maximum root depth (RDMSOL), bounded below by the initial rooting
     * depth (RDI):
     *   RDM = max(RDI, min(RDMCR, RDMSOL))
     *
     * The initial rooting depth (RD) is set to RDI.
     *
     * Root biomass (WRT) is initialised to zero here and must be set
     * separately via set_initial_biomass() after the partitioning module
     * has been initialised, because the root fraction (FR) is computed
     * by the partitioning module.
     *
     * @param[out] rs  Root dynamics state arrays (per grid cell)
     * @param[in]  p   Static crop parameters (RDI, RDMCR)
     * @param[in]  s   Dynamic crop state (DVS -- used if FR lookup needed)
     */
    void initialize(RootDynamicsState &rs, const CropState &p, const CropDynamicState &s)
    {

        auto FRTB = p.tables.at("FRTB");  // Root partitioning fraction table [-] vs DVS

        Kokkos::parallel_for("RootDynamics_Init", rs.nCells, KOKKOS_LAMBDA(const int i) {

            /*-------------------------------------------------------------*/
            /*  Initial rooting depth                                       */
            /*    RDM = max(RDI, min(RDMCR, RDMSOL))                       */
            /*    RDI   = initial rooting depth at emergence [cm]           */
            /*    RDMCR = crop-specific maximum rooting depth [cm]          */
            /*            (e.g. 80--150 cm for cereals)                     */
            /*    RDMSOL = soil-limited maximum rooting depth [cm]          */
            /*             (depends on soil profile; placeholder 150 cm)    */
            /*-------------------------------------------------------------*/
            real RDI = p.p.RDI;                // Initial rooting depth [cm]
            real RDMCR = p.p.RDMCR;            // Crop maximum rooting depth [cm]
            // TODO: Ensure RDMSOL is available as a soil parameter.
            // Currently using a placeholder default of 150 cm.
            // In a full implementation, this should come from a soil map
            // or be passed as a spatially varying array.
            real RDMSOL = 150.0;               // Soil-limited max root depth [cm] (placeholder)

            real rdmax = fmax(RDI, fmin(RDMCR, RDMSOL));

            rs.RDM(i) = rdmax;                 // Maximum attainable rooting depth [cm]
            rs.RD(i) = RDI;                    // Current rooting depth [cm]

            /*-------------------------------------------------------------*/
            /*  Initial root biomass                                        */
            /*    WRT = TDWI * FR(DVS)                                      */
            /*    Recalculate FR from the FRTB table at the initial DVS     */
            /*    since the PartitioningState may not be initialised yet.    */
            /*-------------------------------------------------------------*/
            real DVS = s.DVS(i);
            real FR = Afgen::lookup(FRTB, DVS);
            real TDWI = p.p.TDWI;              // Total initial dry weight [kg/ha]

            rs.WRT(i) = TDWI * FR;             // Living root biomass [kg/ha]
            rs.DWRT(i) = 0.0;                  // Dead root biomass [kg/ha]
            rs.TWRT(i) = rs.WRT(i);            // Total root biomass [kg/ha]
        });
    }

    /**
     * @brief Set initial root biomass after partitioning is initialised.
     *
     * This is an alternative initialisation path that uses the partitioning
     * fraction (FR) computed by the Partitioning module rather than
     * recalculating it from the FRTB table.  It should be called after
     * Partitioning::initialize() has run.
     *
     * @param[out] rs  Root dynamics state arrays
     * @param[in]  p   Static crop parameters (TDWI)
     * @param[in]  FR  Root partitioning fraction [-] (from PartitioningState)
     */
    void set_initial_biomass(RootDynamicsState &rs, const CropState &p, const realArr &FR)
    {
        Kokkos::parallel_for("RootDynamics_InitBiomass", rs.nCells, KOKKOS_LAMBDA(const int i) {
            real TDWI = p.p.TDWI;              // Total initial dry weight [kg/ha]
            rs.WRT(i) = TDWI * FR(i);          // Living root biomass [kg/ha]
            rs.TWRT(i) = rs.WRT(i);            // Total root biomass [kg/ha]
        });
    }

    /*=====================================================================*/
    /*  Rate Calculation                                                   */
    /*=====================================================================*/

    /**
     * @brief Compute daily rate variables for root dynamics.
     *
     * Calculates four rate terms:
     *
     *   1. GRRT = FR * DMI                [kg/(ha d)]
     *      Root biomass growth rate from total dry matter increase (DMI)
     *      and root partitioning fraction (FR).
     *
     *   2. DRRT = WRT * RDRRTB(DVS)       [kg/(ha d)]
     *      Root death rate from the relative death rate table (RDRRTB)
     *      evaluated at the current development stage.
     *
     *   3. GWRT = GRRT - DRRT             [kg/(ha d)]
     *      Net change in living root biomass.
     *
     *   4. RR = min(RDM - RD, RRI)        [cm/d]
     *      Root depth growth rate.  Increases by the daily increment (RRI)
     *      until the maximum depth (RDM) is reached.  Set to zero when
     *      FR = 0 (no root growth after partitioning shifts away from roots).
     *
     * @param[in,out] rs   Root dynamics state (rate arrays are written)
     * @param[in]     p    Static crop parameters (RRI, RDRRTB table)
     * @param[in]     s    Dynamic crop state (DVS)
     * @param[in]     DMI  Total dry matter increase [kg/(ha d)]
     * @param[in]     FR   Root partitioning fraction [-]
     */
    void calc_rates(RootDynamicsState &rs, const CropState &p, const CropDynamicState &s,
                    const realArr &DMI, const realArr &FR)
    {

        auto RDRRTB = p.tables.at("RDRRTB");  // Relative death rate of roots [d^-1] vs DVS

        Kokkos::parallel_for("RootDynamics_Rates", rs.nCells, KOKKOS_LAMBDA(const int i) {

            real DVS = s.DVS(i);               // Development stage [-]

            /*-------------------------------------------------------------*/
            /*  Root biomass growth rate                                    */
            /*    GRRT = FR * DMI                                           */
            /*    DMI = total dry matter increase [kg/(ha d)]               */
            /*    FR  = root partitioning fraction [-] (0..1)               */
            /*-------------------------------------------------------------*/
            rs.GRRT(i) = FR(i) * DMI(i);

            /*-------------------------------------------------------------*/
            /*  Root biomass death rate                                     */
            /*    DRRT = WRT * RDRRTB(DVS)                                  */
            /*    The relative death rate is interpolated from the RDRRTB    */
            /*    table at the current development stage.                   */
            /*-------------------------------------------------------------*/
            real rdr = Afgen::lookup(RDRRTB, DVS);
            rs.DRRT(i) = rs.WRT(i) * rdr;

            /*-------------------------------------------------------------*/
            /*  Net root biomass change                                     */
            /*    GWRT = GRRT - DRRT                                        */
            /*-------------------------------------------------------------*/
            rs.GWRT(i) = rs.GRRT(i) - rs.DRRT(i);

            /*-------------------------------------------------------------*/
            /*  Root depth growth rate                                      */
            /*    RR = min(RDM - RD, RRI)  when FR > 0                     */
            /*    RR = 0                    when FR = 0                     */
            /*    RRI = daily maximum root depth increase [cm/d]            */
            /*          (e.g. 1.0--3.0 cm/d depending on crop)              */
            /*    Root depth cannot exceed RDM (maximum rooting depth).     */
            /*    When FR = 0, root depth growth stops (assimilates are     */
            /*    fully allocated to above-ground organs).                  */
            /*-------------------------------------------------------------*/
            if (FR(i) == 0.0) {
                rs.RR(i) = 0.0;
            } else {
                real potential_growth = p.p.RRI;          // Daily max root depth increase [cm/d]
                real space_available = rs.RDM(i) - rs.RD(i);  // Remaining depth to fill [cm]
                if (space_available < 0.0) space_available = 0.0;

                rs.RR(i) = fmin(space_available, potential_growth);
            }
        });
    }

    /*=====================================================================*/
    /*  State Integration                                                  */
    /*=====================================================================*/

    /**
     * @brief Advance root state variables by one time step.
     *
     * Updates living root biomass (WRT), accumulates dead root biomass (DWRT),
     * recomputes total root biomass (TWRT), and advances the rooting depth (RD).
     *
     * @param[in,out] rs   Root dynamics state arrays
     * @param[in]     delt Time step size [d] (default 1.0)
     */
    void integrate(RootDynamicsState &rs, real delt = 1.0)
    {
        Kokkos::parallel_for("RootDynamics_Integrate", rs.nCells, KOKKOS_LAMBDA(const int i) {

            /*-------------------------------------------------------------*/
            /*  Update root biomass                                         */
            /*    WRT += GWRT * delt   (living root biomass [kg/ha])        */
            /*    DWRT += DRRT * delt  (dead root biomass [kg/ha])          */
            /*    TWRT = WRT + DWRT    (total root biomass [kg/ha])         */
            /*-------------------------------------------------------------*/
            rs.WRT(i) += rs.GWRT(i) * delt;
            rs.DWRT(i) += rs.DRRT(i) * delt;
            rs.TWRT(i) = rs.WRT(i) + rs.DWRT(i);

            /*-------------------------------------------------------------*/
            /*  Update rooting depth                                        */
            /*    RD += RR * delt                                           */
            /*    Root depth increases until it reaches RDM.                */
            /*    Units: [cm] += [cm/d] * [d]                               */
            /*-------------------------------------------------------------*/
            rs.RD(i) += rs.RR(i) * delt;
        });
    }
};

#endif
