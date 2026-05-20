/* -*- mode: c++; c-default-style: "linux" -*- */
/*******************************************************************************
 * @file    StemDynamicsState.h
 * @brief   WOFOST stem dynamics state and rate variable arrays
 *
 * Declares all per-grid-cell state variables and rate variables used by the
 * StemDynamics compute kernel (StemDynamics.h).  Each member is a Kokkos
 * device-accessible 1-D array (realArr) sized to the number of simulation
 * cells (nCells).
 *
 * Variable naming follows WOFOST conventions:
 *   WST  = living stem biomass         [kg/ha]
 *   DWST = dead stem biomass           [kg/ha]
 *   TWST = total stem biomass (W+DW)   [kg/ha]
 *   SAI  = stem area index             [-]  (dimensionless [ha/ha])
 *   GRST = stem growth rate            [kg/(ha d)]
 *   DRST = stem death rate             [kg/(ha d)]
 *   GWST = net stem biomass change     [kg/(ha d)]
 ******************************************************************************/

#ifndef _STEM_DYNAMICS_STATE_H_
#define _STEM_DYNAMICS_STATE_H_

#include "../define.h"

/**
 * @class StemDynamicsState
 * @brief Holds all stem-related state and rate arrays for the crop model.
 *
 * Arrays are allocated on the default Kokkos memory space and are accessed
 * inside parallel kernels.  The allocate() method must be called once before
 * any computation.
 */
class StemDynamicsState {
public:
    int nCells;  ///< Number of spatial grid cells

    /*=====================================================================*/
    /*  State Variables (per cell)                                         */
    /*=====================================================================*/

    realArr WST;   ///< Living stem dry biomass [kg/ha]
                   ///<   Typical range: 0 -- 8000 kg/ha depending on crop.

    realArr DWST;  ///< Dead stem dry biomass [kg/ha]
                   ///<   Cumulative biomass lost to senescence and death.
                   ///<   Monotonically increasing during the growing season.

    realArr TWST;  ///< Total stem dry biomass (living + dead) [kg/ha]
                   ///<   TWST = WST + DWST.
                   ///<   Used for carbon balance verification.

    realArr SAI;   ///< Stem area index [-]  (dimensionless [ha/ha])
                   ///<   Green area of stems per unit ground area.
                   ///<   Computed as WST * specific stem area (SSATB).
                   ///<   Typically much smaller than LAI (0 -- 1.5).

    /*=====================================================================*/
    /*  Rate Variables (per cell)                                          */
    /*=====================================================================*/

    realArr GRST;  ///< Stem biomass growth rate [kg/(ha d)]
                   ///<   GRST = ADMI * FS  (above-ground DM increase x stem
                   ///<   partitioning fraction).

    realArr DRST;  ///< Stem biomass death rate [kg/(ha d)]
                   ///<   DRST = WST * RDRSTB(DVS)  (living biomass x relative
                   ///<   death rate from the RDRSTB table).

    realArr GWST;  ///< Net change in living stem biomass [kg/(ha d)]
                   ///<   GWST = GRST - DRST - REALLOC_ST
                   ///<   Includes losses from reallocation to storage organs.

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

        // State variables
        WST  = realArr("WST", n);
        DWST = realArr("DWST", n);
        TWST = realArr("TWST", n);
        SAI  = realArr("SAI", n);

        // Rate variables
        GRST = realArr("GRST", n);
        DRST = realArr("DRST", n);
        GWST = realArr("GWST", n);
    }
};

#endif
