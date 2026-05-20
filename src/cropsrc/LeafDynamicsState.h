/* -*- mode: c++; c-default-style: "linux" -*- */
/*******************************************************************************
 * @file    LeafDynamicsState.h
 * @brief   WOFOST leaf dynamics state and rate variable arrays
 *
 * Declares all per-grid-cell state variables and rate variables used by the
 * LeafDynamics compute kernel (LeafDynamics.h).  Each member is a Kokkos
 * device-accessible 1-D or 2-D array (realArr / realArr2) sized to the
 * number of simulation cells (nCells).
 *
 * Variable naming follows WOFOST conventions:
 *   W  = living (wet) biomass          [kg/ha]
 *   DW = dead biomass                  [kg/ha]
 *   TW = total biomass (W + DW)        [kg/ha]
 *   LAI = leaf area index              [ha/ha]  (dimensionless, expressed as [-])
 *   GR = growth rate                   [kg/(ha d)]
 *   DR = death rate                    [kg/(ha d)]
 *   GW = net biomass change (GR - DR)  [kg/(ha d)]
 *
 * The leaf age-class system tracks cohorts of leaves created each day.
 * Each class stores its own biomass, specific leaf area, and physiological
 * age.  This allows the model to senesce the oldest leaves first and to
 * apply ageing-based death (leaves older than SPAN days die).
 ******************************************************************************/

#ifndef _LEAF_DYNAMICS_STATE_H_
#define _LEAF_DYNAMICS_STATE_H_

#include "../define.h"

/**
 * @class LeafDynamicsState
 * @brief Holds all leaf-related state and rate arrays for the crop model.
 *
 * Arrays are allocated on the default Kokkos memory space and are accessed
 * inside parallel kernels.  The allocate() method must be called once before
 * any computation.
 */
class LeafDynamicsState {
public:
    int nCells;                                  ///< Number of spatial grid cells

    static const int MAX_LEAF_CLASSES = 365;     ///< Maximum number of leaf age
                                                 ///< classes (one per day; 365 d
                                                 ///< covers a full growing season)

    /*=====================================================================*/
    /*  Scalar LAI State Variables (per cell)                              */
    /*=====================================================================*/

    realArr LAIEM;   ///< LAI at crop emergence [-]
                     ///<   Prescribed minimum LAI at initialisation (e.g. 0.07).
                     ///<   Typical range: 0.01 -- 0.15.

    realArr LASUM;   ///< Total leaf area from all age classes [-]
                     ///<   Computed as sum( LV[c] * SLA[c] ).
                     ///<   This is the leaf-only contribution to the full LAI.

    realArr LAIEXP;  ///< LAI under exponential growth regime [-]
                     ///<   Tracks LAI while the crop is in the exponential
                     ///<   expansion phase (LAIEXP < 6.0).  Once LAIEXP >= 6.0,
                     ///<   the source-limited growth pathway takes over.

    realArr LAIMAX;  ///< Maximum LAI reached during the growing season [-]
                     ///<   Monotonically increasing; never decreases.
                     ///<   Typical range: 1.0 -- 8.0 depending on crop.

    /*=====================================================================*/
    /*  Biomass State Variables (per cell)                                 */
    /*=====================================================================*/

    realArr WLV;     ///< Living leaf dry biomass [kg/ha]
                     ///<   Sum of biomass across all living leaf age classes.
                     ///<   Typical range: 0 -- 5000 kg/ha.

    realArr DWLV;    ///< Dead leaf dry biomass [kg/ha]
                     ///<   Cumulative biomass removed by senescence.
                     ///<   Monotonically increasing (dead leaves are not recycled).

    realArr TWLV;    ///< Total leaf dry biomass (living + dead) [kg/ha]
                     ///<   TWLV = WLV + DWLV.
                     ///<   Useful for carbon budget verification.

    realArr LAI;     ///< Leaf area index [-]  (green area index from leaves)
                     ///<   LAI = LASUM + SAI + PAI (leaf + stem + pod areas).
                     ///<   Dimensionless [ha/ha].  Typical range: 0 -- 8.

    /*=====================================================================*/
    /*  Leaf Age-Class Arrays (2-D: nCells x MAX_LEAF_CLASSES)            */
    /*=====================================================================*/

    intArr num_leaf_classes;  ///< Number of active leaf age classes per cell [-]
                              ///<   Ranges from 1 at emergence up to
                              ///<   MAX_LEAF_CLASSES.  Classes are indexed 0
                              ///<   (newest) to n-1 (oldest).

    realArr2 LV;     ///< Leaf biomass per age class [kg/ha]
                     ///<   LV(i, c) = biomass of leaf class c in cell i.
                     ///<   Dimension 2 layout: [nCells, MAX_LEAF_CLASSES].

    realArr2 SLA;    ///< Specific leaf area per age class [ha/kg]
                     ///<   SLA(i, c) converts biomass to area for class c.
                     ///<   Typical range: 0.0015 -- 0.0030 ha/kg
                     ///<   (equivalent to 15 -- 30 m^2/kg).

    realArr2 LVAGE;  ///< Physiological age per age class [-]
                     ///<   Accumulated physiological ageing (in effective days).
                     ///<   Leaves with LVAGE > SPAN are considered senescent.

    /*=====================================================================*/
    /*  Rate Variables (per cell)                                          */
    /*=====================================================================*/

    realArr GRLV;    ///< Leaf biomass growth rate [kg/(ha d)]
                     ///<   GRLV = ADMI * FL  (above-ground DM increase x leaf
                     ///<   partitioning fraction).

    realArr DSLV1;   ///< Leaf death rate due to water stress [kg/(ha d)]
                     ///<   DSLV1 = WLV * (1 - RFTRA) * PERDL.
                     ///<   RFTRA: transpiration reduction factor (0 = full
                     ///<   stress, 1 = no stress).  PERDL: max relative death
                     ///<   rate per day (typically 0.03).

    realArr DSLV2;   ///< Leaf death rate due to self-shading [kg/(ha d)]
                     ///<   Activated when LAI exceeds a critical threshold
                     ///<   (LAICR = 3.2 / KDIF).  Maximum rate 3% of WLV per day.

    realArr DSLV3;   ///< Leaf death rate due to frost [kg/(ha d)]
                     ///<   Placeholder; currently always 0 in this implementation.

    realArr DSLV;    ///< Maximum stress-related leaf death rate [kg/(ha d)]
                     ///<   DSLV = max(DSLV1, DSLV2, DSLV3).

    realArr DALV;    ///< Leaf death rate due to ageing [kg/(ha d)]
                     ///<   Sum of biomass in all classes with LVAGE > SPAN.

    realArr DRLV;    ///< Total leaf death rate [kg/(ha d)]
                     ///<   DRLV = max(DSLV, DALV).  The larger of stress or
                     ///<   ageing death is applied to avoid double-counting.

    realArr SLAT;    ///< Effective specific leaf area for new growth [ha/kg]
                     ///<   Interpolated from the SLATB table at the current DVS,
                     ///<   then adjusted if exponential growth limits expansion.

    realArr FYSAGE;  ///< Daily physiological ageing increment [-]
                     ///<   FYSAGE = max(0, (TEMP - TBASE) / (35 - TBASE)).
                     ///<   Accumulated in LVAGE during integrate().

    realArr GLAIEX;  ///< Potential LAI increase from exponential growth [-/d]
                     ///<   GLAIEX = LAIEXP * RGRLAI * dteff (when LAIEXP < 6).

    realArr GLASOL;  ///< Potential LAI increase from source-limited growth [-/d]
                     ///<   GLASOL = GRLV * SLAT.

    /*=====================================================================*/
    /*  Memory Allocation                                                  */
    /*=====================================================================*/

    /**
     * @brief Allocate all state and rate arrays for n grid cells.
     *
     * 1-D arrays are sized to n.  2-D leaf class arrays are sized to
     * n x MAX_LEAF_CLASSES.  All arrays are zero-initialised by Kokkos.
     *
     * @param n  Number of spatial grid cells
     */
    void allocate(int n) {
        nCells = n;

        // Scalar LAI state
        LAIEM  = realArr("LAIEM", n);
        LASUM  = realArr("LASUM", n);
        LAIEXP = realArr("LAIEXP", n);
        LAIMAX = realArr("LAIMAX", n);

        // Biomass state
        WLV  = realArr("WLV", n);
        DWLV = realArr("DWLV", n);
        TWLV = realArr("TWLV", n);
        LAI  = realArr("LAI", n);

        // Number of active leaf classes per cell
        num_leaf_classes = intArr("num_leaf_classes", n);

        // 2-D leaf class arrays: [nCells, MAX_LEAF_CLASSES]
        // LayoutRight (row-major): cell i's classes are stored contiguously,
        // which is cache-friendly for the per-cell class iteration patterns
        // used in calc_rates() and integrate().
        LV    = realArr2("LV_Classes", n, MAX_LEAF_CLASSES);
        SLA   = realArr2("SLA_Classes", n, MAX_LEAF_CLASSES);
        LVAGE = realArr2("LVAGE_Classes", n, MAX_LEAF_CLASSES);

        // Rate variables
        GRLV   = realArr("GRLV", n);
        DSLV1  = realArr("DSLV1", n);
        DSLV2  = realArr("DSLV2", n);
        DSLV3  = realArr("DSLV3", n);
        DSLV   = realArr("DSLV", n);
        DALV   = realArr("DALV", n);
        DRLV   = realArr("DRLV", n);
        SLAT   = realArr("SLAT", n);
        FYSAGE = realArr("FYSAGE", n);
        GLAIEX = realArr("GLAIEX", n);
        GLASOL = realArr("GLASOL", n);
    }
};

#endif
