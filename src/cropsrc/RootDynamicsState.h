/* -*- mode: c++; c-default-style: "linux" -*- */
/*******************************************************************************
 * @file    RootDynamicsState.h
 * @brief   WOFOST root dynamics state and rate variable arrays
 *
 * Declares all per-grid-cell state variables and rate variables used by the
 * RootDynamics compute kernel (RootDynamics.h).  Each member is a Kokkos
 * device-accessible 1-D array (realArr) sized to the number of simulation
 * cells (nCells).
 *
 * Variable naming follows WOFOST conventions:
 *   RDM   [cm]        Maximum attainable rooting depth
 *   RD    [cm]        Current rooting depth
 *   WRT   [kg/ha]     Living root dry biomass
 *   DWRT  [kg/ha]     Dead root dry biomass
 *   TWRT  [kg/ha]     Total root dry biomass (WRT + DWRT)
 *   RR    [cm/d]      Rate of root depth increase
 *   GRRT  [kg/(ha d)] Root biomass growth rate
 *   DRRT  [kg/(ha d)] Root biomass death rate
 *   GWRT  [kg/(ha d)] Net root biomass change rate
 *
 * Note on units: Rooting depth uses [cm] internally in WOFOST, following
 * the convention of the original FORTRAN implementation.  Some external
 * interfaces may convert to [m].
 ******************************************************************************/

#ifndef _ROOT_DYNAMICS_STATE_H_
#define _ROOT_DYNAMICS_STATE_H_

#include "../define.h"

/**
 * @class RootDynamicsState
 * @brief Holds all root-related state and rate arrays for the crop model.
 *
 * Arrays are allocated on the default Kokkos memory space and are accessed
 * inside parallel kernels.  The allocate() method must be called once before
 * any computation.
 */
class RootDynamicsState {
public:
    int nCells;  ///< Number of spatial grid cells

    /*=====================================================================*/
    /*  Rate Variables (per cell)                                          */
    /*=====================================================================*/

    realArr RR;    ///< Rate of root depth increase [cm/d]
                   ///<   RR = min(RDM - RD, RRI) when FR > 0, else 0.
                   ///<   RRI is the crop-specific daily maximum root depth
                   ///<   increase (typically 1.0--3.0 cm/d).
                   ///<   RR decreases as RD approaches RDM.

    realArr GRRT;  ///< Root biomass growth rate [kg/(ha d)]
                   ///<   GRRT = FR * DMI  (root partitioning fraction x total
                   ///<   dry matter increase).

    realArr DRRT;  ///< Root biomass death rate [kg/(ha d)]
                   ///<   DRRT = WRT * RDRRTB(DVS)  (living root biomass x
                   ///<   relative death rate from the RDRRTB table).

    realArr GWRT;  ///< Net change in living root biomass [kg/(ha d)]
                   ///<   GWRT = GRRT - DRRT.

    /*=====================================================================*/
    /*  State Variables (per cell)                                         */
    /*=====================================================================*/

    realArr RDM;   ///< Maximum attainable rooting depth [cm]
                   ///<   Determined at initialization as:
                   ///<     RDM = max(RDI, min(RDMCR, RDMSOL))
                   ///<   where RDI = initial depth, RDMCR = crop max depth,
                   ///<   RDMSOL = soil-limited max depth.
                   ///<   Typical range: 30--150 cm depending on crop and soil.

    realArr RD;    ///< Current rooting depth [cm]
                   ///<   Increases daily by RR until RDM is reached.
                   ///<   Starts at RDI (initial rooting depth at emergence).
                   ///<   Typical range: 10--150 cm.

    realArr WRT;   ///< Living root dry biomass [kg/ha]
                   ///<   WRT = TDWI * FR at initialization; updated by
                   ///<   WRT += GWRT * delt each time step.
                   ///<   Typical range: 0 -- 2000 kg/ha.

    realArr DWRT;  ///< Dead root dry biomass [kg/ha]
                   ///<   Accumulated dead root biomass: DWRT += DRRT * delt.
                   ///<   Monotonically increasing during the growing season.

    realArr TWRT;  ///< Total root dry biomass (living + dead) [kg/ha]
                   ///<   TWRT = WRT + DWRT.
                   ///<   Used for carbon balance verification.

    /*=====================================================================*/
    /*  Memory Allocation                                                  */
    /*=====================================================================*/

    /**
     * @brief Allocate all state and rate arrays for n grid cells.
     *
     * All arrays are 1-D and sized to n.  Zero-initialised by Kokkos.
     *
     * @param n  Number of spatial grid cells
     */
    void allocate(int n) {
        nCells = n;

        // Rate variables
        RR   = realArr("RR", n);
        GRRT = realArr("GRRT", n);
        DRRT = realArr("DRRT", n);
        GWRT = realArr("GWRT", n);

        // State variables
        RDM  = realArr("RDM", n);
        RD   = realArr("RD", n);
        WRT  = realArr("WRT", n);
        DWRT = realArr("DWRT", n);
        TWRT = realArr("TWRT", n);
    }
};

#endif
