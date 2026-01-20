/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _EVAPOTRANSPIRATION_STATE_H_
#define _EVAPOTRANSPIRATION_STATE_H_

#include "../define.h"

class EvapotranspirationState {
public:
    int nCells;

    // Rate Variables (per day)
    realArr TRA;   // Actual transpiration rate [cm/day]
    realArr TRAMX; // Maximum transpiration rate [cm/day]
    realArr EVS;   // Actual evaporation rate from soil surface [cm/day]
    realArr E0;    // Potential evaporation from open water [cm/day]
    realArr ES0;   // Potential evaporation from bare soil [cm/day]
    realArr ET0;   // Potential evapotranspiration from reference crop [cm/day]
    realArr RFTRA; // Reduction factor for transpiration (water stress) [-]

    // Allocation
    void allocate(int n) {
        nCells = n;
        TRA   = realArr("TRA", n);
        TRAMX = realArr("TRAMX", n);
        EVS   = realArr("EVS", n);
        E0    = realArr("E0", n); 
        ES0   = realArr("ES0", n); 
        ET0   = realArr("ET0", n); 
        RFTRA = realArr("RFTRA", n);
    }
};

#endif