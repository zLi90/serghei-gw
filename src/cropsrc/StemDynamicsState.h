/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _STEM_DYNAMICS_STATE_H_
#define _STEM_DYNAMICS_STATE_H_

#include "../define.h"

// 存储茎动态相关的状态变量
class StemDynamicsState {
public:
    int nCells;

    // State Variables
    realArr WST;  // Weight of living stems [kg ha-1]
    realArr DWST; // Weight of dead stems [kg ha-1]
    realArr TWST; // Total weight of stems [kg ha-1]
    realArr SAI;  // Stem Area Index [-]

    // Rate Variables (Stored here for output or intermediate usage)
    realArr GRST; // Growth rate stem biomass [kg ha-1 d-1]
    realArr DRST; // Death rate stem biomass [kg ha-1 d-1]
    realArr GWST; // Net change in stem biomass [kg ha-1 d-1]

    // Allocation
    void allocate(int n) {
        nCells = n;
        WST  = realArr("WST", n);
        DWST = realArr("DWST", n);
        TWST = realArr("TWST", n);
        SAI  = realArr("SAI", n);

        GRST = realArr("GRST", n);
        DRST = realArr("DRST", n);
        GWST = realArr("GWST", n);
    }
};

#endif