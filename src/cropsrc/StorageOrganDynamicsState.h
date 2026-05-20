/* -*- mode: c++; c-default-style: "linux" -*- */
/*******************************************************************************
 * @file    StorageOrganDynamicsState.h
 * @brief   WOFOST storage organ dynamics state and rate variable arrays
 *
 * Declares all per-grid-cell state variables and rate variables used by the
 * StorageOrganDynamics compute kernel (StorageOrganDynamics.h).  Each member
 * is a Kokkos device-accessible 1-D array (realArr) sized to the number of
 * simulation cells (nCells).
 *
 * Storage organs represent the economically harvested part of the crop:
 * grains in cereals, tubers in potatoes, pods in legumes, etc.
 *
 * Variable naming follows WOFOST conventions:
 *   WSO   [kg/ha]     Living storage organ dry biomass (grain/tuber yield)
 *   DWSO  [kg/ha]     Dead storage organ biomass (always 0 in WOFOST)
 *   TWSO  [kg/ha]     Total storage organ biomass (WSO + DWSO)
 *   PAI   [-]         Pod area index (pod area per unit ground area)
 *   GRSO  [kg/(ha d)] Storage organ growth rate
 *   DRSO  [kg/(ha d)] Storage organ death rate (always 0)
 *   GWSO  [kg/(ha d)] Net storage organ biomass change rate
 ******************************************************************************/

#ifndef _STORAGE_ORGAN_DYNAMICS_STATE_H_
#define _STORAGE_ORGAN_DYNAMICS_STATE_H_

#include "../define.h"

/**
 * @class StorageOrganDynamicsState
 * @brief Holds all storage-organ-related state and rate arrays for the crop model.
 *
 * Arrays are allocated on the default Kokkos memory space and are accessed
 * inside parallel kernels.  The allocate() method must be called once before
 * any computation.
 */
class StorageOrganDynamicsState {
public:
    int nCells;  ///< Number of spatial grid cells

    /*=====================================================================*/
    /*  State Variables (per cell)                                         */
    /*=====================================================================*/

    realArr WSO;   ///< Living storage organ dry biomass [kg/ha]
                   ///<   This is the primary crop yield variable.
                   ///<   WSO = TDWI * (1 - FR) * FO at initialization.
                   ///<   Updated by WSO += GWSO * delt each time step.
                   ///<   Typical range: 0 -- 12000 kg/ha at harvest
                   ///<   (depends on crop; wheat ~4000--8000, maize ~8000--12000).

    realArr DWSO;  ///< Dead storage organ dry biomass [kg/ha]
                   ///<   Always zero in WOFOST because storage organs do not
                   ///<   senesce during the simulation.  Included for
                   ///<   consistency with the organ state variable pattern.

    realArr TWSO;  ///< Total storage organ dry biomass (living + dead) [kg/ha]
                   ///<   TWSO = WSO + DWSO.  Effectively equal to WSO since
                   ///<   DWSO = 0.  Used for carbon balance verification.

    realArr PAI;   ///< Pod area index [-]  (dimensionless [ha/ha])
                   ///<   Green area of pods/storage organs per unit ground area.
                   ///<   Computed as WSO * SPA (specific pod area).
                   ///<   Typically small compared to LAI (0 -- 2.0).

    /*=====================================================================*/
    /*  Rate Variables (per cell)                                          */
    /*=====================================================================*/

    realArr GRSO;  ///< Storage organ biomass growth rate [kg/(ha d)]
                   ///<   GRSO = ADMI * FO  (above-ground DM increase x
                   ///<   storage organ partitioning fraction).  This is the
                   ///<   primary grain-filling rate.

    realArr DRSO;  ///< Storage organ biomass death rate [kg/(ha d)]
                   ///<   Always zero in WOFOST.  Storage organs do not die.

    realArr GWSO;  ///< Net change in storage organ biomass [kg/(ha d)]
                   ///<   GWSO = GRSO - DRSO + REALLOC_SO
                   ///<   Includes gains from reallocation of stem and leaf
                   ///<   reserves during grain filling.

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
        WSO  = realArr("WSO", n);
        DWSO = realArr("DWSO", n);
        TWSO = realArr("TWSO", n);
        PAI  = realArr("PAI", n);

        // Rate variables
        GRSO = realArr("GRSO", n);
        DRSO = realArr("DRSO", n);
        GWSO = realArr("GWSO", n);
    }
};

#endif
