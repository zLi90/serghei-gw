/* -*- mode: c++; c-default-style: "linux" -*- */
/*******************************************************************************
 * @file    PartitioningState.h
 * @brief   WOFOST biomass partitioning factor arrays
 *
 * Declares the per-grid-cell partitioning fraction arrays used by the
 * Partitioning compute kernel (Partitioning.h).  These four dimensionless
 * fractions control how daily assimilates are distributed among the four
 * organ pools (roots, leaves, stems, storage organs).
 *
 * Partitioning convention (WOFOST standard):
 *   - FR is the fraction of TOTAL dry matter allocated to roots.
 *   - FL, FS, FO are fractions of ABOVE-GROUND dry matter (1 - FR)
 *     allocated to leaves, stems, and storage organs respectively.
 *   - Mass balance: FR + (FL + FS + FO) * (1 - FR) = 1.0
 *
 * All fractions are functions of development stage (DVS) and are
 * interpolated from crop-specific AFGEN tables at each time step.
 ******************************************************************************/

#ifndef _PARTITIONING_STATE_H_
#define _PARTITIONING_STATE_H_

#include "../define.h"

/**
 * @class PartitioningState
 * @brief Holds the four biomass partitioning fraction arrays for the crop model.
 *
 * Arrays are allocated on the default Kokkos memory space and are accessed
 * inside parallel kernels.  The allocate() method must be called once before
 * any computation.
 */
class PartitioningState {
public:
    int nCells;  ///< Number of spatial grid cells

    /*=====================================================================*/
    /*  Partitioning Factors (per cell)                                    */
    /*=====================================================================*/

    realArr FR;   ///< Fraction of total dry matter allocated to roots [-]
                   ///<   Interpolated from the FRTB table at the current DVS.
                   ///<   Typical pattern: high at emergence (0.4--0.6),
                   ///<   decreasing through vegetative growth, near 0 after
                   ///<   flowering.  Range: [0, 1].

    realArr FL;   ///< Fraction of above-ground dry matter allocated to leaves [-]
                   ///<   Interpolated from the FLTB table at the current DVS.
                   ///<   Typical pattern: dominant early (0.6--0.8), drops
                   ///<   sharply around flowering (DVS ~ 1.0).  Range: [0, 1].

    realArr FS;   ///< Fraction of above-ground dry matter allocated to stems [-]
                   ///<   Interpolated from the FSTB table at the current DVS.
                   ///<   Typical pattern: moderate during vegetative growth
                   ///<   (0.2--0.5), drops to near 0 after flowering.  Range: [0, 1].

    realArr FO;   ///< Fraction of above-ground dry matter allocated to storage organs [-]
                   ///<   Interpolated from the FOTB table at the current DVS.
                   ///<   Typical pattern: zero before flowering, dominant after
                   ///<   flowering (0.8--1.0) during grain filling.  Range: [0, 1].

    /*=====================================================================*/
    /*  Memory Allocation                                                  */
    /*=====================================================================*/

    /**
     * @brief Allocate all partitioning factor arrays for n grid cells.
     *
     * All arrays are 1-D and sized to n.  Zero-initialised by Kokkos.
     *
     * @param n  Number of spatial grid cells
     */
    void allocate(int n) {
        nCells = n;
        FR = realArr("FR", n);
        FL = realArr("FL", n);
        FS = realArr("FS", n);
        FO = realArr("FO", n);
    }
};

#endif
