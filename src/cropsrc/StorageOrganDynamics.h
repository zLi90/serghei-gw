/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _STORAGE_ORGAN_DYNAMICS_H_
#define _STORAGE_ORGAN_DYNAMICS_H_

#include "../define.h"
#include "CropState.h"        // Static parameters
#include "CropDynamicState.h" // Dynamic variables (ADMI, FO, FR)
#include "StorageOrganDynamicsState.h" // Output variables
#include <cmath>

class StorageOrganDynamics {

public:
    // Initialize storage organ states
    // 对应 Python: initialize()
    void initialize(StorageOrganDynamicsState &sods, const CropState &p, 
                    const realArr &FR, const realArr &FO) {
        
        Kokkos::parallel_for("StorageOrganDynamics_Init", sods.nCells, KOKKOS_LAMBDA(const int i) {
            
            real TDWI = p.p.TDWI;
            real SPA = p.p.SPA;
            real fr = FR(i);
            real fo = FO(i);
            
            // Initial storage organ biomass
            // WSO  = (params.TDWI * (1-FR)) * FO
            sods.WSO(i) = (TDWI * (1.0 - fr)) * fo;
            sods.DWSO(i) = 0.0;
            sods.TWSO(i) = sods.WSO(i) + sods.DWSO(i);
            
            // Initial Pod Area Index
            // PAI = WSO * params.SPA
            sods.PAI(i) = sods.WSO(i) * SPA;
        });
    }

    // Calculate rates
    // 对应 Python: calc_rates()
    // Need ADMI (Above-ground dry matter increase), FO (Fraction to organs), REALLOC_SO (Reallocation)
    void calc_rates(StorageOrganDynamicsState &sods, const CropState &p,
                    const realArr &ADMI, const realArr &FO, const realArr &REALLOC_SO) {
        
        Kokkos::parallel_for("StorageOrganDynamics_Rates", sods.nCells, KOKKOS_LAMBDA(const int i) {
            
            real admi = ADMI(i);
            real fo = FO(i);
            real realloc_so = REALLOC_SO(i);

            // Growth rate organs
            // rates.GRSO = ADMI * FO
            sods.GRSO(i) = admi * fo;
            
            // Death rate organs (Always 0 in WOFOST)
            // rates.DRSO = 0.0
            sods.DRSO(i) = 0.0;
            
            // Net change
            // rates.GWSO = rates.GRSO - rates.DRSO + k.REALLOC_SO
            // Note: REALLOC_SO is added here (gain from stem/leaf translocation)
            sods.GWSO(i) = sods.GRSO(i) - sods.DRSO(i) + realloc_so;
        });
    }

    // Integrate states
    // 对应 Python: integrate()
    void integrate(StorageOrganDynamicsState &sods, const CropState &p, real delt = 1.0) {
        
        Kokkos::parallel_for("StorageOrganDynamics_Integrate", sods.nCells, KOKKOS_LAMBDA(const int i) {
            
            real SPA = p.p.SPA;

            // Storage organ biomass (living, dead, total)
            sods.WSO(i) += sods.GWSO(i) * delt;
            sods.DWSO(i) += sods.DRSO(i) * delt;
            sods.TWSO(i) = sods.WSO(i) + sods.DWSO(i);

            // Calculate Pod Area Index (PAI)
            // states.PAI = states.WSO * params.SPA
            sods.PAI(i) = sods.WSO(i) * SPA;
        });
    }
};

#endif