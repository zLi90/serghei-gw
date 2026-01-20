/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _RESPIRATION_STATE_H_
#define _RESPIRATION_STATE_H_

#include "../define.h"

class RespirationState {
public:
    int nCells;

    // Rate Variables
    realArr PMRES; // Potential maintenance respiration rate

    // Allocation
    void allocate(int n) {
        nCells = n;
        PMRES = realArr("PMRES", n);
    }
};

#endif