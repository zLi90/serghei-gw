/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _PARTITIONING_H_
#define _PARTITIONING_H_

#include "../define.h"
#include "CropState.h"         // Static parameters
#include "CropDynamicState.h"  // Dynamic variables (DVS)
#include "PartitioningState.h" // Output variables
#include "Afgen.h"             // Table lookup helper
#include <cmath>
#include <iostream>

class Partitioning
{

public:
    // Initialize partitioning factors based on initial DVS
    // 对应 Python: initialize()
    void initialize(PartitioningState &ps, const CropState &p, const CropDynamicState &s)
    {

        // Access Tables
        auto FRTB = p.tables.at("FRTB");
        auto FLTB = p.tables.at("FLTB");
        auto FSTB = p.tables.at("FSTB");
        auto FOTB = p.tables.at("FOTB");

        Kokkos::parallel_for("Partitioning_Init", ps.nCells, KOKKOS_LAMBDA(const int i) {
            real DVS = s.DVS(i);
            
            ps.FR(i) = Afgen::lookup(FRTB, DVS);
            ps.FL(i) = Afgen::lookup(FLTB, DVS);
            ps.FS(i) = Afgen::lookup(FSTB, DVS);
            ps.FO(i) = Afgen::lookup(FOTB, DVS);

            // DEBUG: Print partitioning factors
#if DEBUG_CROP_GROWTH_MODEL            
            printf("[PARTITIONING-INIT] cell=%d, DVS=%.2f, FR=%.4f, FL=%.4f, FS=%.4f, FO=%.4f\n", 
                   i, DVS, ps.FR(i), ps.FL(i), ps.FS(i), ps.FO(i));
#endif            



            // Check partitioning (inline for GPU compatibility)
            // checksum = FR+(FL+FS+FO)*(1.-FR) - 1.
            real checksum = ps.FR(i) + (ps.FL(i) + ps.FS(i) + ps.FO(i)) * (1.0 - ps.FR(i)) - 1.0;
            
            // Note: Printing from GPU/Device is possible but can be messy. 
            // In strict production code, we might set an error flag.
            // keeping simple print for now if error is large.
            if (std::abs(checksum) >= 0.0001) {
                // printf is supported in CUDA/Kokkos kernels usually
                printf("Error in partitioning at cell %d! Checksum: %f\n", i, checksum);
            } });
    }

    // Update partitioning factors (Integrate step in Python)
    // 对应 Python: integrate()
    // Note: Python code also has calc_rates, but it just returns state.
    // The actual update logic is in integrate.
    void update(PartitioningState &ps, const CropState &p, const CropDynamicState &s)
    {

        // Access Tables
        auto FRTB = p.tables.at("FRTB");
        auto FLTB = p.tables.at("FLTB");
        auto FSTB = p.tables.at("FSTB");
        auto FOTB = p.tables.at("FOTB");

        Kokkos::parallel_for("Partitioning_Update", ps.nCells, KOKKOS_LAMBDA(const int i) {
            real DVS = s.DVS(i);

            ps.FR(i) = Afgen::lookup(FRTB, DVS);
            ps.FL(i) = Afgen::lookup(FLTB, DVS);
            ps.FS(i) = Afgen::lookup(FSTB, DVS);
            ps.FO(i) = Afgen::lookup(FOTB, DVS);

#if DEBUG_CROP_GROWTH_MODEL

            if (i == 0)
            {

                printf("[PARTITION-DEBUG] DVS=%.4f, FL=%.4f, FR=%.4f, FS=%.4f, FO=%.4f\n",
                       DVS, ps.FL(i), ps.FR(i), ps.FS(i), ps.FO(i));

                if (DVS > 0.95 && DVS < 1.05)
                {
                    printf("[PARTITION-CRITICAL] DVS is near 1.0! FL=%.6f\n", ps.FL(i));
                }
            }
#endif

            // Check partitioning
            real checksum = ps.FR(i) + (ps.FL(i) + ps.FS(i) + ps.FO(i)) * (1.0 - ps.FR(i)) - 1.0;

            // if (std::abs(checksum) >= 0.0001) { ... }
        });
    }
};

#endif
