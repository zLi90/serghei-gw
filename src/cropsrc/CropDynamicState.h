/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _CROP_DYNAMIC_STATE_H_
#define _CROP_DYNAMIC_STATE_H_

#include "../define.h"

namespace CropStage {
    const int EMERGING = 0;
    const int VEGETATIVE = 1;
    const int REPRODUCTIVE = 2;
    const int MATURE = 3;
}

class CropDynamicState {
public:
    int nCells;

    // --- 核心状态变量 (积分变量) ---
    // Phenology
    realArr DVS;   // Development stage [-]
    realArr TSUM;  // Temperature sum [C d]
    realArr TSUME; // Temperature sum for emergence [C d]
    intArr STAGE;  // Current stage ID
    realArr VERN;  // Vernalisation days
    intArr ISVERNALISED; // Flag
    realArr VERNR; // Vernalisation Rate (Intermediate)
    realArr DTSUM; // TSUM increase rate (Intermediate)
    realArr DVR;   // Development rate (Intermediate)

    // Biomass States (from Leaf, Stem, Root, Storage Dynamics)
    realArr WLV; // Living leaf biomass (from LeafDynamicsState)
    realArr DWLV; // Dead leaf biomass (from LeafDynamicsState)
    realArr TWLV; // Total leaf biomass (from LeafDynamicsState)
    realArr LAI; // Leaf Area Index (from LeafDynamicsState)
    
    realArr WST; // Living stem biomass (from StemDynamicsState)
    realArr DWST; // Dead stem biomass (from StemDynamicsState)
    realArr TWST; // Total stem biomass (from StemDynamicsState)
    realArr SAI; // Stem Area Index (from StemDynamicsState)

    realArr WRT; // Living root biomass (from RootDynamicsState)
    realArr DWRT; // Dead root biomass (from RootDynamicsState)
    realArr TWRT; // Total root biomass (from RootDynamicsState)
    realArr RD;    // Rooting Depth [cm]

    realArr WSO; // Living storage organ biomass (from StorageOrganDynamicsState)
    realArr DWSO; // Dead storage organ biomass (from StorageOrganDynamicsState)
    realArr TWSO; // Total storage organ biomass (from StorageOrganDynamicsState)
    realArr PAI; // Pod Area Index (from StorageOrganDynamicsState)

    // Main Rates & Reductions
    realArr GASS;      // Actual assimilation rate corrected for water stress
    realArr MRES;      // Actual maintenance respiration rate
    realArr ASRC;      // Net available assimilates
    realArr DMI;       // Total dry matter increase
    realArr ADMI;      // Aboveground dry matter increase
    realArr RFTRA;     // Water stress reduction factor for transpiration
    realArr REALLOC_LV; // Reallocation from leaves
    realArr REALLOC_ST; // Reallocation from stems
    realArr REALLOC_SO; // Reallocation to storage organs
    
    // Partitioning factors
    realArr FL;        // Fraction of biomass to leaves
    realArr FR;        // Fraction of biomass to roots
    
    // Internal control flag
    intArr skip_other_modules; 

    // Allocation using Kokkos Views (realArr defined in define.h)
    void allocate(int n) {
        nCells = n;
        // Phenology
        DVS   = realArr("DVS", n); TSUM  = realArr("TSUM", n); TSUME = realArr("TSUME", n);
        STAGE = intArr("STAGE", n); VERN = realArr("VERN", n); ISVERNALISED = intArr("ISVERNALISED", n);
        VERNR = realArr("VERNR", n); DTSUM = realArr("DTSUM", n); DVR = realArr("DVR", n);
        
        // Biomass
        WLV = realArr("WLV", n); DWLV = realArr("DWLV", n); TWLV = realArr("TWLV", n); LAI = realArr("LAI", n);
        WST = realArr("WST", n); DWST = realArr("DWST", n); TWST = realArr("TWST", n); SAI = realArr("SAI", n);
        WRT = realArr("WRT", n); DWRT = realArr("DWRT", n); TWRT = realArr("TWRT", n); RD = realArr("RD", n);
        WSO = realArr("WSO", n); DWSO = realArr("DWSO", n); TWSO = realArr("TWSO", n); PAI = realArr("PAI", n);

        // Rates and reductions
        GASS = realArr("GASS", n); MRES = realArr("MRES", n); ASRC = realArr("ASRC", n);
        DMI = realArr("DMI", n); ADMI = realArr("ADMI", n); RFTRA = realArr("RFTRA", n);
        REALLOC_LV = realArr("REALLOC_LV", n); REALLOC_ST = realArr("REALLOC_ST", n); REALLOC_SO = realArr("REALLOC_SO", n);
        
        // Partitioning factors
        FL = realArr("FL", n); FR = realArr("FR", n);
        
        skip_other_modules = intArr("skip_other_modules", n);
    }
};
#endif
