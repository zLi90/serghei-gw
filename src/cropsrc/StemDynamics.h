/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _STEM_DYNAMICS_H_
#define _STEM_DYNAMICS_H_

#include "../define.h"
#include "CropState.h"        // Static parameters
#include "CropDynamicState.h" // Dynamic variables (DVS, ADMI, FS, FR)
#include "StemDynamicsState.h" // Output variables
#include "Afgen.h"            // Table lookup helper
#include <cmath>

class StemDynamics {

public:
    // Initialize stem states
    // 对应 Python: initialize()
    void initialize(StemDynamicsState &sds, const CropState &p, const CropDynamicState &s, 
                    const realArr &FR, const realArr &FS) {
        
        auto SSATB = p.tables.at("SSATB");

        Kokkos::parallel_for("StemDynamics_Init", sds.nCells, KOKKOS_LAMBDA(const int i) {
            
            real TDWI = p.p.TDWI;
            real fr = FR(i);
            real fs = FS(i);
            
            // Set initial stem biomass
            // WST  = (params.TDWI * (1-FR)) * FS
            sds.WST(i) = (TDWI * (1.0 - fr)) * fs;
            sds.DWST(i) = 0.0;
            sds.TWST(i) = sds.WST(i) + sds.DWST(i);
            
            // Initial Stem Area Index
            // SAI = WST * params.SSATB(DVS)
            real DVS = s.DVS(i);
            real ssa = Afgen::lookup(SSATB, DVS);
            sds.SAI(i) = sds.WST(i) * ssa;
        });
    }

    // Calculate rates
    // 对应 Python: calc_rates()
    // Need ADMI (Above-ground dry matter increase), FS (Fraction to stems), REALLOC_ST (Reallocation)
    void calc_rates(StemDynamicsState &sds, const CropState &p, const CropDynamicState &s,
                    const realArr &ADMI, const realArr &FS, const realArr &REALLOC_ST) {
        
        auto RDRSTB = p.tables.at("RDRSTB");

        Kokkos::parallel_for("StemDynamics_Rates", sds.nCells, KOKKOS_LAMBDA(const int i) {
            
            real DVS = s.DVS(i);
            real admi = ADMI(i);
            real fs = FS(i);
            real realloc_st = REALLOC_ST(i);

            // Growth rate stems
            // rates.GRST = ADMI * FS
            sds.GRST(i) = admi * fs;
            
            // Death rate stems
            // rates.DRST = params.RDRSTB(DVS) * states.WST
            real rdr = Afgen::lookup(RDRSTB, DVS);
            sds.DRST(i) = rdr * sds.WST(i);
            
            // Net change
            // rates.GWST = rates.GRST - rates.DRST - k.REALLOC_ST
            sds.GWST(i) = sds.GRST(i) - sds.DRST(i) - realloc_st;
        });
    }

    // Integrate states
    // 对应 Python: integrate()
    void integrate(StemDynamicsState &sds, const CropState &p, const CropDynamicState &s, real delt = 1.0) {
        
        auto SSATB = p.tables.at("SSATB");

        Kokkos::parallel_for("StemDynamics_Integrate", sds.nCells, KOKKOS_LAMBDA(const int i) {
            
            // Stem biomass (living, dead, total)
            sds.WST(i) += sds.GWST(i) * delt;
            sds.DWST(i) += sds.DRST(i) * delt;
            sds.TWST(i) = sds.WST(i) + sds.DWST(i);

            // Calculate Stem Area Index (SAI)
            // states.SAI = states.WST * params.SSATB(DVS)
            real DVS = s.DVS(i);
            real ssa = Afgen::lookup(SSATB, DVS);
            sds.SAI(i) = sds.WST(i) * ssa;
        });
    }
};

#endif