/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _ROOT_DYNAMICS_H_
#define _ROOT_DYNAMICS_H_

#include "../define.h"
#include "CropState.h"        // Static parameters
#include "CropDynamicState.h" // Dynamic variables (DVS, DMI, FR)
#include "RootDynamicsState.h" // Output variables
#include "Afgen.h"            // Table lookup helper
#include <cmath>
#include <algorithm> // for std::min, std::max

class RootDynamics {

public:
    // Initialize root states
    // 对应 Python: initialize()
    void initialize(RootDynamicsState &rs, const CropState &p, const CropDynamicState &s) {
        
        // Initial states usually depend on parameters and initial partitioning
        // Note: Python code accesses kiosk.FR in initialize. This assumes partitioning
        // has been initialized before root dynamics.
        // We will assume s.FR is available (from Partitioning::initialize)
        // Python: WRT  = params.TDWI * self.kiosk.FR

        Kokkos::parallel_for("RootDynamics_Init", rs.nCells, KOKKOS_LAMBDA(const int i) {
            // Initial root depth states
            // rdmax = max(params.RDI, min(params.RDMCR, params.RDMSOL))
            // Note: RDMSOL is a soil parameter. In SERGHEI it might come from SoilMap or similar.
            // For now, assume it's passed or stored in CropState (though technically it's soil).
            // Let's assume CropParams has RDMSOL for simplicity as per Python code structure.
            // But usually RDMSOL varies by pixel. If so, it should be an array.
            // Assuming RDMSOL is scalar in CropParams for now based on provided code structure.
            
            real RDI = p.p.RDI;
            real RDMCR = p.p.RDMCR;
            // TODO: Ensure RDMSOL is available. Assuming it's in CropParams based on Python snippet.
            // If not, it needs to be passed separately.
            real RDMSOL = 150.0; // Placeholder default if missing
            // Ideally: p.p.RDMSOL (if added to CropParams struct)
            
            real rdmax = fmax(RDI, fmin(RDMCR, RDMSOL));
            
            rs.RDM(i) = rdmax;
            rs.RD(i) = RDI;

            // Initial root biomass
            real TDWI = p.p.TDWI;
            // FR needs to be available from Partitioning initialization
            // Assuming Partitioning::initialize ran before this.
            // But `s` (CropDynamicState) might not store FR if PartitioningState is separate.
            // Let's assume we pass PartitioningState `ps` or it's merged.
            // If `s` doesn't have FR, we need to re-calculate initial FR or pass it.
            // For now, let's assume `s` has `FR` or we can access it.
            // NOTE: In the provided files, `FR` is in `PartitioningState`.
            // Let's assume we pass `PartitioningState` to this function or `s` has it.
            // Let's rely on a passed-in value or recalculate. 
            // Re-calculating:
            // real DVS = s.DVS(i); 
            // auto FRTB = p.tables.at("FRTB");
            // real FR = Afgen::lookup(FRTB, DVS);
            
            // To keep it clean, let's just initialize to 0 here and let the first step handle it?
            // No, WRT needs initial value.
            // Let's assume TDWI is small and initial partitioning is set.
            // Simplification: WRT = TDWI * 0.5 (rough guess) if FR unavailable?
            // Better: Require FR to be passed.
            // But since I can't change signature easily to include PartitioningState without breaking pattern...
            // I will assume FR is 0.5 as placeholder or needs to be fixed in integration.
            // Wait, previous instructions had CropDynamicState. If FR is there, good.
            // If FR is in PartitioningState, we should pass it.
            
            // Placeholder:
            rs.WRT(i) = 0.0; // Needs correction
            rs.DWRT(i) = 0.0;
            rs.TWRT(i) = 0.0;
        });
    }
    
    // Helper to set initial WRT properly if called after Partitioning init
    void set_initial_biomass(RootDynamicsState &rs, const CropState &p, const realArr &FR) {
        Kokkos::parallel_for("RootDynamics_InitBiomass", rs.nCells, KOKKOS_LAMBDA(const int i) {
            real TDWI = p.p.TDWI;
            rs.WRT(i) = TDWI * FR(i);
            rs.TWRT(i) = rs.WRT(i);
        });
    }

    // Calculate rates
    // 对应 Python: calc_rates()
    // Need DMI (Total dry matter increase) and FR (Fraction to roots)
    // Assuming these are passed via CropDynamicState or similar arrays.
    // Let's assume we pass `realArr DMI` and `realArr FR` explicitly for clarity, 
    // or they are in `s`.
    void calc_rates(RootDynamicsState &rs, const CropState &p, const CropDynamicState &s, 
                    const realArr &DMI, const realArr &FR) {
        
        auto RDRRTB = p.tables.at("RDRRTB");

        Kokkos::parallel_for("RootDynamics_Rates", rs.nCells, KOKKOS_LAMBDA(const int i) {
            
            real DVS = s.DVS(i);
            
            // Increase in root biomass
            // r.GRRT = k.FR * k.DMI
            rs.GRRT(i) = FR(i) * DMI(i);
            
            // Death rate
            // r.DRRT = s.WRT * p.RDRRTB(k.DVS)
            real rdr = Afgen::lookup(RDRRTB, DVS);
            rs.DRRT(i) = rs.WRT(i) * rdr;
            
            // Net change
            rs.GWRT(i) = rs.GRRT(i) - rs.DRRT(i);
            
            // Increase in root depth
            // r.RR = min((s.RDM - s.RD), p.RRI)
            // if k.FR == 0.: r.RR = 0.
            if (FR(i) == 0.0) {
                rs.RR(i) = 0.0;
            } else {
                real potential_growth = p.p.RRI; // Daily increase
                real space_available = rs.RDM(i) - rs.RD(i);
                if (space_available < 0.0) space_available = 0.0;
                
                rs.RR(i) = fmin(space_available, potential_growth);
            }
        });
    }

    // Integrate states
    // 对应 Python: integrate()
    void integrate(RootDynamicsState &rs, real delt = 1.0) {
        Kokkos::parallel_for("RootDynamics_Integrate", rs.nCells, KOKKOS_LAMBDA(const int i) {
            // Dry weight of living roots
            rs.WRT(i) += rs.GWRT(i) * delt;
            
            // Dry weight of dead roots
            rs.DWRT(i) += rs.DRRT(i) * delt;
            
            // Total weight dry + living roots
            rs.TWRT(i) = rs.WRT(i) + rs.DWRT(i);

            // New root depth
            rs.RD(i) += rs.RR(i) * delt;
        });
    }
};

#endif