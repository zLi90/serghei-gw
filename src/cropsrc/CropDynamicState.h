/**
 * @file CropDynamicState.h
 * @brief Dynamic (time-varying) crop state variables for the WOFOST model.
 *
 * This file defines:
 *   - CropStage namespace: enumeration of crop growth stages (emerging,
 *     vegetative, reproductive, mature).
 *   - CropDynamicState class: per-cell dynamic state arrays (Kokkos Views)
 *     that evolve during the simulation. These include phenological state,
 *     area indices, carbon fluxes, and water-stress reduction factors.
 *
 * All state arrays are allocated with extent nCells (one value per surface
 * grid cell) and are accessed in parallel kernels.
 *
 * Design note: biomass organ states (WLV, WST, WRT, WSO, etc.) are stored
 * in their respective organ dynamics modules (LeafDynamicsState,
 * StemDynamicsState, RootDynamicsState, StorageOrganDynamicsState) to
 * avoid duplication. This class retains only the shared cross-module
 * states (LAI, SAI, PAI, RD) and the main carbon balance fluxes.
 *
 * @see Wofost72.h for the orchestrator that drives these states.
 */

/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _CROP_DYNAMIC_STATE_H_
#define _CROP_DYNAMIC_STATE_H_

#include "../define.h"

/**
 * @namespace CropStage
 * @brief Integer constants representing crop phenological growth stages.
 *
 * These match the WOFOST convention:
 *   - EMERGING     (0): seed has been sown but has not yet emerged.
 *   - VEGETATIVE   (1): crop is in the vegetative phase (DVS 0--1).
 *   - REPRODUCTIVE (2): crop is in the reproductive phase (DVS 1--2).
 *   - MATURE       (3): crop has reached physiological maturity.
 */
namespace CropStage {
    const int EMERGING     = 0;  /**< Pre-emergence: no assimilation, growth, or ET. */
    const int VEGETATIVE   = 1;  /**< Vegetative growth: leaf expansion, tillering. */
    const int REPRODUCTIVE = 2;  /**< Generative growth: flowering, grain fill. */
    const int MATURE       = 3;  /**< Physiological maturity: growth stops. */
}

/**
 * @class CropDynamicState
 * @brief Per-cell arrays of dynamic crop state variables.
 *
 * Each Kokkos View has extent nCells (one element per surface grid cell).
 * The allocate() method creates all views at once.
 *
 * Variable naming follows WOFOST conventions (e.g., DVS = Development Stage).
 *
 * Units legend:
 *   [-]       = dimensionless
 *   [deg.C.d] = degree-days (temperature sum)
 *   [d]       = days
 *   [ha/ha]   = leaf area index
 *   [cm]      = rooting depth
 *   [kg/ha/d] = dry matter flux per day
 */
class CropDynamicState {
public:
    int nCells;  /**< Number of surface grid cells with active crop simulation. */

    /* ================================================================== */
    /* Phenological State Variables                                       */
    /* ================================================================== */

    realArr DVS;           /**< Development stage [-].
                                Range: 0.0 (emergence) to ~2.0 (maturity).
                                Vegetative phase: 0 < DVS <= 1.
                                Reproductive phase: 1 < DVS <= 2. */
    realArr TSUM;          /**< Accumulated effective temperature sum [deg.C.d].
                                Drives development from emergence through
                                anthesis (TSUM1) and maturity (TSUM1+TSUM2). */
    realArr TSUME;         /**< Accumulated temperature sum for emergence
                                [deg.C.d]. Crop emerges when TSUME >= TSUMEM. */
    intArr  STAGE;         /**< Current phenological stage ID [-].
                                See CropStage namespace for values. */
    realArr VERN;          /**< Accumulated vernalisation days [d].
                                Relevant only for winter crops (e.g., winter wheat). */
    intArr  ISVERNALISED;  /**< Vernalisation completion flag [-].
                                0 = not yet vernalised, 1 = fully vernalised. */

    /* ================================================================== */
    /* Phenological Rate Variables (intermediate, per-timestep)           */
    /* ================================================================== */

    realArr VERNR;         /**< Vernalisation rate [d/d]. Increment to VERN. */
    realArr DTSUM;         /**< Increment to temperature sum [deg.C.d/d]. */
    realArr DVR;           /**< Development rate [d-1]. Increment to DVS. */

    /*
     * Note: Biomass organ states (WLV, WST, WRT, WSO, DWLV, DWST, DWRT,
     * DWSO, TWLV, TWST, TWRT, TWSO) are stored exclusively in their
     * respective organ dynamics modules (LeafDynamicsState, StemDynamicsState,
     * etc.) to avoid data duplication. Only the shared area-index and
     * rooting-depth variables remain here.
     */

    /* ================================================================== */
    /* Area Index and Rooting Depth (shared across modules)               */
    /* ================================================================== */

    realArr LAI;  /**< Leaf Area Index [ha/ha].
                        Range: 0.0 -- ~10.0. Used by Assimilation,
                        Evapotranspiration, and LeafDynamics modules. */
    realArr SAI;  /**< Stem Area Index [ha/ha].
                        Range: 0.0 -- ~2.0. Used by LeafDynamics for
                        green-area computation. */
    realArr PAI;  /**< Pod Area Index [ha/ha].
                        Range: 0.0 -- ~1.0. Used by LeafDynamics for
                        green-area computation. */
    realArr RD;   /**< Rooting depth [cm].
                        Range: RDI -- RDMCR. Used by RootDynamics and
                        Evapotranspiration to define the root zone. */

    /* ================================================================== */
    /* Main Carbon Fluxes and Reduction Factors                           */
    /* ================================================================== */

    realArr GASS;       /**< Actual gross assimilation rate (water-stress
                             corrected) [kg CH2O/ha/d]. */
    realArr MRES;       /**< Actual maintenance respiration rate
                             [kg CH2O/ha/d]. Limited to GASS if GASS < PMRES. */
    realArr ASRC;       /**< Net available assimilates [kg CH2O/ha/d].
                             ASRC = GASS - MRES. */
    realArr DMI;        /**< Total dry matter increase [kg DM/ha/d].
                             DMI = CVF * ASRC, where CVF is the weighted
                             conversion efficiency. */
    realArr ADMI;       /**< Above-ground dry matter increase [kg DM/ha/d].
                             ADMI = (1 - FR) * DMI. */
    realArr RFTRA;      /**< Water stress reduction factor for transpiration
                             [-]. Range: 0.0 -- 1.0. 1.0 = no stress. */

    /* ================================================================== */
    /* Reallocation Fluxes (organ-to-organ biomass transfer)              */
    /* ================================================================== */

    realArr REALLOC_LV;  /**< Reallocation flux from leaves to storage organs
                              [kg DM/ha/d]. Currently initialised to 0. */
    realArr REALLOC_ST;  /**< Reallocation flux from stems to storage organs
                              [kg DM/ha/d]. Currently initialised to 0. */
    realArr REALLOC_SO;  /**< Reallocation flux to storage organs
                              [kg DM/ha/d]. Currently initialised to 0. */

    /* ================================================================== */
    /* Internal Control                                                   */
    /* ================================================================== */

    intArr skip_other_modules;  /**< Per-cell flag: 1 = skip assimilation,
                                     partitioning, and organ dynamics
                                     (used during pre-emergence). 0 = normal
                                     operation. */

    /* ================================================================== */
    /* Memory Allocation                                                  */
    /* ================================================================== */

    /**
     * @brief Allocate all dynamic state arrays for the given number of cells.
     *
     * Each Kokkos View is created with extent n. Must be called once before
     * the simulation begins.
     *
     * @param n  Number of surface grid cells (nCells).
     */
    void allocate(int n) {
        nCells = n;

        /* Phenological states */
        DVS   = realArr("DVS", n);
        TSUM  = realArr("TSUM", n);
        TSUME = realArr("TSUME", n);
        STAGE = intArr("STAGE", n);
        VERN  = realArr("VERN", n);
        ISVERNALISED = intArr("ISVERNALISED", n);

        /* Phenological rates (intermediate) */
        VERNR = realArr("VERNR", n);
        DTSUM = realArr("DTSUM", n);
        DVR   = realArr("DVR", n);

        /*
         * Biomass organ states are allocated by their respective modules.
         * Only the shared cross-module states are allocated here.
         */
        LAI = realArr("LAI", n);
        SAI = realArr("SAI", n);
        PAI = realArr("PAI", n);
        RD  = realArr("RD", n);

        /* Carbon fluxes and reduction factors */
        GASS  = realArr("GASS", n);
        MRES  = realArr("MRES", n);
        ASRC  = realArr("ASRC", n);
        DMI   = realArr("DMI", n);
        ADMI  = realArr("ADMI", n);
        RFTRA = realArr("RFTRA", n);

        /* Reallocation fluxes */
        REALLOC_LV = realArr("REALLOC_LV", n);
        REALLOC_ST = realArr("REALLOC_ST", n);
        REALLOC_SO = realArr("REALLOC_SO", n);

        /* Control flag */
        skip_other_modules = intArr("skip_other_modules", n);
    }
};

#endif /* _CROP_DYNAMIC_STATE_H_ */
