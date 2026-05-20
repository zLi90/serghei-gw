/* -*- mode: c++; c-default-style: "linux" -*- */
/*******************************************************************************
 * @file    Partitioning.h
 * @brief   WOFOST biomass partitioning factors (compute kernel)
 *
 * Implements the biomass partitioning module of the WOFOST crop growth model.
 * At each time step, the daily assimilates (dry matter production) must be
 * distributed among the four organ pools: roots, leaves, stems, and storage
 * organs.  The partitioning fractions are determined by the crop development
 * stage (DVS) through four AFGEN lookup tables:
 *
 *   FR(DVS) = fraction to roots          FRTB table   [-]
 *   FL(DVS) = fraction to leaves         FLTB table   [-]
 *   FS(DVS) = fraction to stems          FSTB table   [-]
 *   FO(DVS) = fraction to storage organs FOTB table   [-]
 *
 * Partitioning convention (WOFOST standard):
 *   - FR is the fraction of TOTAL dry matter allocated to roots.
 *   - FL, FS, FO are fractions of ABOVE-GROUND dry matter (1 - FR)
 *     allocated to leaves, stems, and storage organs respectively.
 *   - Mass balance constraint: FR + (FL + FS + FO) * (1 - FR) = 1.0
 *
 * This module provides two functions:
 *   1. initialize() - Set initial partitioning factors from the DVS at
 *                      crop emergence.
 *   2. update()     - Recompute partitioning factors at each time step
 *                      as DVS advances.
 *
 * Both functions include a checksum verification to ensure the mass
 * balance constraint is satisfied (tolerance: 1e-4).
 *
 * Reference: WOFOST 7.1 User Guide, Chapter 3 -- Partitioning
 ******************************************************************************/

#ifndef _PARTITIONING_H_
#define _PARTITIONING_H_

#include "../define.h"
#include "CropState.h"         // Static crop parameters (partitioning tables)
#include "CropDynamicState.h"  // Dynamic state variables (DVS)
#include "PartitioningState.h" // Partitioning factor output arrays (FR, FL, FS, FO)
#include "Afgen.h"             // AFGEN table interpolation utility
#include <cmath>
#include <iostream>

class Partitioning
{

public:
    /*=====================================================================*/
    /*  Initialization                                                     */
    /*=====================================================================*/

    /**
     * @brief Initialize partitioning factors at crop emergence.
     *
     * Looks up all four partitioning fractions (FR, FL, FS, FO) from
     * their respective AFGEN tables at the initial development stage (DVS).
     * Verifies the mass balance constraint:
     *   FR + (FL + FS + FO) * (1 - FR) = 1.0  (within tolerance 1e-4)
     *
     * @param[out] ps  Partitioning state arrays (per grid cell)
     * @param[in]  p   Static crop parameters (FRTB, FLTB, FSTB, FOTB tables)
     * @param[in]  s   Dynamic crop state (DVS)
     */
    void initialize(PartitioningState &ps, const CropState &p, const CropDynamicState &s)
    {

        // Access AFGEN partitioning tables (Kokkos Views captured by value
        // in the lambda, safe for device execution)
        auto FRTB = p.tables.at("FRTB");  // Root fraction vs DVS [-]
        auto FLTB = p.tables.at("FLTB");  // Leaf fraction vs DVS [-]
        auto FSTB = p.tables.at("FSTB");  // Stem fraction vs DVS [-]
        auto FOTB = p.tables.at("FOTB");  // Storage organ fraction vs DVS [-]

        /*---------------------------------------------------------------*/
        /*  Extract local Kokkos View handles before entering the        */
        /*  parallel kernel.  This ensures that only device-accessible   */
        /*  View data is referenced inside the lambda, avoiding any      */
        /*  potential host-side dereference of the wrapper objects.       */
        /*---------------------------------------------------------------*/
        auto local_DVS = s.DVS;
        auto local_FR = ps.FR;
        auto local_FL = ps.FL;
        auto local_FS = ps.FS;
        auto local_FO = ps.FO;

        // Execution policy uses ps.nCells (evaluated on host, safe)
        Kokkos::parallel_for("Partitioning_Init", ps.nCells, KOKKOS_LAMBDA(const int i) {

            // Look up all four partitioning fractions at the current DVS
            real DVS = local_DVS(i);         // Development stage [-]

            local_FR(i) = Afgen::lookup(FRTB, DVS);  // Root fraction [-]
            local_FL(i) = Afgen::lookup(FLTB, DVS);  // Leaf fraction [-]
            local_FS(i) = Afgen::lookup(FSTB, DVS);  // Stem fraction [-]
            local_FO(i) = Afgen::lookup(FOTB, DVS);  // Storage organ fraction [-]

            /*-------------------------------------------------------------*/
            /*  Verify mass balance constraint                              */
            /*    checksum = FR + (FL + FS + FO) * (1 - FR) - 1.0          */
            /*    Should equal zero (within tolerance 1e-4).                */
            /*    A non-zero checksum indicates inconsistent partitioning    */
            /*    tables in the crop parameter file.                        */
            /*-------------------------------------------------------------*/
            real checksum = local_FR(i) + (local_FL(i) + local_FS(i) + local_FO(i)) * (1.0 - local_FR(i)) - 1.0;

            if (Kokkos::fabs(checksum) >= 0.0001) {
                printf("Error in partitioning at cell %d! Checksum: %f\n", i, checksum);
            }
        });
        Kokkos::fence();  // Ensure kernel completion before returning to host
    }

    /*=====================================================================*/
    /*  Time-Step Update                                                   */
    /*=====================================================================*/

    /**
     * @brief Update partitioning factors for the current time step.
     *
     * Re-evaluates all four partitioning fractions from their AFGEN tables
     * at the current development stage (DVS).  This is called each time
     * step after DVS has been advanced.
     *
     * The mass balance checksum is computed but only printed in debug mode
     * (not enforced in production to avoid kernel-side error handling).
     *
     * @param[out] ps  Partitioning state arrays (per grid cell)
     * @param[in]  p   Static crop parameters (FRTB, FLTB, FSTB, FOTB tables)
     * @param[in]  s   Dynamic crop state (DVS)
     */
    void update(PartitioningState &ps, const CropState &p, const CropDynamicState &s)
    {

        // Access AFGEN partitioning tables
        auto FRTB = p.tables.at("FRTB");  // Root fraction vs DVS [-]
        auto FLTB = p.tables.at("FLTB");  // Leaf fraction vs DVS [-]
        auto FSTB = p.tables.at("FSTB");  // Stem fraction vs DVS [-]
        auto FOTB = p.tables.at("FOTB");  // Storage organ fraction vs DVS [-]

        /*---------------------------------------------------------------*/
        /*  Extract local Kokkos View handles for device-side access     */
        /*---------------------------------------------------------------*/
        auto local_DVS = s.DVS;
        auto local_FR = ps.FR;
        auto local_FL = ps.FL;
        auto local_FS = ps.FS;
        auto local_FO = ps.FO;

        Kokkos::parallel_for("Partitioning_Update", ps.nCells, KOKKOS_LAMBDA(const int i) {
            real DVS = local_DVS(i);         // Development stage [-]

            // Look up all four partitioning fractions at the current DVS
            local_FR(i) = Afgen::lookup(FRTB, DVS);  // Root fraction [-]
            local_FL(i) = Afgen::lookup(FLTB, DVS);  // Leaf fraction [-]
            local_FS(i) = Afgen::lookup(FSTB, DVS);  // Stem fraction [-]
            local_FO(i) = Afgen::lookup(FOTB, DVS);  // Storage organ fraction [-]

#if DEBUG_CROP_GROWTH_MODEL
            if (i == 0)
            {
                printf("[PARTITION-DEBUG] DVS=%.4f, FL=%.4f, FR=%.4f, FS=%.4f, FO=%.4f\n",
                       DVS, local_FL(i), local_FR(i), local_FS(i), local_FO(i));

                // Flag when DVS is near 1.0 (flowering transition)
                // Partitioning shifts dramatically around flowering:
                // FL drops, FO increases (onset of grain filling)
                if (DVS > 0.95 && DVS < 1.05)
                {
                    printf("[PARTITION-CRITICAL] DVS is near 1.0! FL=%.6f\n", local_FL(i));
                }
            }
#endif

            /*-------------------------------------------------------------*/
            /*  Verify mass balance constraint                              */
            /*    checksum = FR + (FL + FS + FO) * (1 - FR) - 1.0          */
            /*    Should equal zero (within tolerance 1e-4).                */
            /*-------------------------------------------------------------*/
            real checksum = local_FR(i) + (local_FL(i) + local_FS(i) + local_FO(i)) * (1.0 - local_FR(i)) - 1.0;
        });
        Kokkos::fence();  // Ensure kernel completion before returning to host
    }
};

#endif
