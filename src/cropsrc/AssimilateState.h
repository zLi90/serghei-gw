/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _ASSIMILATE_STATE_H_
#define _ASSIMILATE_STATE_H_

#include "../define.h"

// 存储同化过程相关的速率变量
class AssimilateState {
public:
    int nCells;

    // Output Rate Variables
    realArr PGASS; // Potential assimilation rate [kg CH2O ha-1 d-1]
    
    // Allocation
    void allocate(int n) {
        nCells = n;
        PGASS = realArr("PGASS", n);
    }
};

// 辅助结构：存储 Astro 函数的返回值
// SArray 用于这种小型、固定大小的结构非常合适，如果我们需要在 View 中存储它。
// 但这里它只是函数返回值，使用 struct 即可。
struct AstroVars {
    real DAYL;   
    real DAYLP;  
    real SINLD;  
    real COSLD;  
    real DIFPP;  
    real ATMTR;  
    real DSINBE; 
    real ANGOT;  
};

#endif