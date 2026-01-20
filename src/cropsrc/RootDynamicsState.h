/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _ROOT_DYNAMICS_STATE_H_
#define _ROOT_DYNAMICS_STATE_H_

#include "../define.h"

class RootDynamicsState {
public:
    int nCells;

    // Rate Variables specific to this module (states are in CropDynamicState usually, but kept here for modularity logic)
    // Note: In Wofost72.h, we link these back to CropDynamicState or use them as intermediates.
    // If these duplicate CropDynamicState, they act as the "Local State" buffer.
    
    realArr RR;   // Growth rate root depth [cm d-1]
    realArr GRRT; // Growth rate root biomass [kg ha-1 d-1]
    realArr DRRT; // Death rate root biomass [kg ha-1 d-1]
    realArr GWRT; // Net change in root biomass [kg ha-1 d-1]

    // Additional state variables needed by the code
    realArr RDM;  // Maximum rooting depth [cm]
    realArr RD;   // Current rooting depth [cm]
    realArr WRT;  // Living root biomass [kg ha-1]
    realArr DWRT; // Dead root biomass [kg ha-1]
    realArr TWRT; // Total root biomass [kg ha-1]

    // Allocation
    void allocate(int n) {
        nCells = n;
        RR = realArr("RR", n);
        GRRT = realArr("GRRT", n);
        DRRT = realArr("DRRT", n);
        GWRT = realArr("GWRT", n);
        
        // Additional variables
        RDM = realArr("RDM", n);
        RD = realArr("RD", n);
        WRT = realArr("WRT", n);
        DWRT = realArr("DWRT", n);
        TWRT = realArr("TWRT", n);
    }
};

#endif
