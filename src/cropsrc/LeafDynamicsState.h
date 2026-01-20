/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _LEAF_DYNAMICS_STATE_H_
#define _LEAF_DYNAMICS_STATE_H_

#include "../define.h"

// 存储叶片动态相关的状态变量
class LeafDynamicsState {
public:
    int nCells;
    static const int MAX_LEAF_CLASSES = 365; 

    // Scalar State Variables (per cell) - pointers to global CropDynamicState or independent?
    // Independent copies for module-specific internal state
    realArr LAIEM;  // LAI at emergence [-]
    realArr LASUM;  // Total leaf area (sum of LV*SLA) [-]
    realArr LAIEXP; // LAI under exponential growth [-]
    realArr LAIMAX; // Maximum LAI reached [-]
    
    // Additional state variables needed by the code
    realArr WLV;    // Living leaf biomass [kg ha-1]
    realArr DWLV;   // Dead leaf biomass [kg ha-1]
    realArr TWLV;   // Total leaf biomass [kg ha-1]
    realArr LAI;    // Leaf Area Index [-]

    // Arrays for leaf classes
    intArr num_leaf_classes; 

    // 使用 realArr2 (defined in define.h as Kokkos::View<real**, ...>)
    realArr2 LV;    // Leaf biomass per class [kg ha-1]
    realArr2 SLA;   // Specific leaf area per class [ha kg-1]
    realArr2 LVAGE; // Leaf age per class [d]

    // Local Rate Variables
    realArr GRLV;   
    realArr DSLV1;  
    realArr DSLV2;  
    realArr DSLV3;  
    realArr DSLV;   
    realArr DALV;   
    realArr DRLV;   
    realArr SLAT;   
    realArr FYSAGE; 
    realArr GLAIEX; 
    realArr GLASOL; 

    // Allocation
    void allocate(int n) {
        nCells = n;
        LAIEM = realArr("LAIEM", n);
        LASUM = realArr("LASUM", n);
        LAIEXP = realArr("LAIEXP", n);
        LAIMAX = realArr("LAIMAX", n);
        
        // Additional variables
        WLV = realArr("WLV", n);
        DWLV = realArr("DWLV", n);
        TWLV = realArr("TWLV", n);
        LAI = realArr("LAI", n);

        num_leaf_classes = intArr("num_leaf_classes", n);
        
        // Allocate 2D arrays: [nCells, MAX_LEAF_CLASSES]
        // Note: LayoutRight is default in define.h, so (n, m) is row-major. 
        // This is good for accessing leaf classes of a single cell contiguously? 
        // No, LayoutRight means stride is on last dim. (i, j) -> i*m + j. 
        // So leaf classes for cell i are contiguous. Good for CPU cache, maybe coalescing issues on GPU if threads access (i, c).
        // SERGHEI define.h uses LayoutRight.
        LV = realArr2("LV_Classes", n, MAX_LEAF_CLASSES);
        SLA = realArr2("SLA_Classes", n, MAX_LEAF_CLASSES);
        LVAGE = realArr2("LVAGE_Classes", n, MAX_LEAF_CLASSES);

        GRLV = realArr("GRLV", n);
        DSLV1 = realArr("DSLV1", n);
        DSLV2 = realArr("DSLV2", n);
        DSLV3 = realArr("DSLV3", n);
        DSLV = realArr("DSLV", n);
        DALV = realArr("DALV", n);
        DRLV = realArr("DRLV", n);
        SLAT = realArr("SLAT", n);
        FYSAGE = realArr("FYSAGE", n);
        GLAIEX = realArr("GLAIEX", n);
        GLASOL = realArr("GLASOL", n);
    }
};

#endif
