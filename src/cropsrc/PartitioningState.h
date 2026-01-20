/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _PARTITIONING_STATE_H_
#define _PARTITIONING_STATE_H_

#include "../define.h"

class PartitioningState {
public:
    int nCells;

    // Partitioning Factors
    realArr FR; 
    realArr FL; 
    realArr FS; 
    realArr FO; 

    // Allocation
    void allocate(int n) {
        nCells = n;
        FR = realArr("FR", n);
        FL = realArr("FL", n);
        FS = realArr("FS", n);
        FO = realArr("FO", n);
    }
};

#endif