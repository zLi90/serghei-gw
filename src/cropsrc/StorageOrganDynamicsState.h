/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _STORAGE_ORGAN_DYNAMICS_STATE_H_
#define _STORAGE_ORGAN_DYNAMICS_STATE_H_

#include "../define.h"

// 存储储藏器官动态相关的状态变量
class StorageOrganDynamicsState {
public:
    int nCells;

    // State Variables
    realArr WSO;  // Weight of living storage organs [kg ha-1]
    realArr DWSO; // Weight of dead storage organs [kg ha-1]
    realArr TWSO; // Total weight of storage organs [kg ha-1]
    realArr PAI;  // Pod Area Index [-]

    // Rate Variables (Stored here for output or intermediate usage)
    realArr GRSO; // Growth rate storage organs [kg ha-1 d-1]
    realArr DRSO; // Death rate storage organs [kg ha-1 d-1]
    realArr GWSO; // Net change in storage organ biomass [kg ha-1 d-1]

    // Allocation
    void allocate(int n) {
        nCells = n;
        WSO  = realArr("WSO", n);
        DWSO = realArr("DWSO", n);
        TWSO = realArr("TWSO", n);
        PAI  = realArr("PAI", n);

        GRSO = realArr("GRSO", n);
        DRSO = realArr("DRSO", n);
        GWSO = realArr("GWSO", n);
    }
};

#endif