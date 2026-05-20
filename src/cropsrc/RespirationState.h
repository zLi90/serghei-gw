/* -*- mode: c++; c-default-style: "linux" -*- */

/*******************************************************************************
 * @file    RespirationState.h
 * @brief   State container for crop maintenance respiration variables.
 *
 * @details This header defines the RespirationState class, which holds the
 *          output arrays produced by the Respiration module each simulation
 *          day. It is designed for use with Kokkos device memory, and provides
 *          an allocate() method to size the arrays for the computational grid.
 *
 *          The primary output variable is PMRES (Potential Maintenance
 *          Respiration rate), which represents the total carbohydrate cost
 *          of maintaining all living crop tissue, after accounting for
 *          senescence reduction and temperature effects.
 *
 *          In the WOFOST assimilation balance, PMRES is subtracted from
 *          gross assimilation (PGASS) to determine the net available
 *          carbohydrates for growth:
 *
 *            ASRC = PGASS - PMRES    (assimilation surplus, may be negative)
 *
 *          If ASRC < 0, the crop must mobilize reserve carbohydrates
 *          (e.g., from stems or leaves) to meet its maintenance demand.
 *
 * @note    Units for PMRES: [kg CH2O / ha / d]
 *          (kilograms of carbohydrate equivalents per hectare per day)
 *
 * @see     Respiration.h  for the calculation routine that fills PMRES
 ******************************************************************************/

#ifndef _RESPIRATION_STATE_H_
#define _RESPIRATION_STATE_H_

#include "../define.h"

/**
 * @class RespirationState
 * @brief Kokkos-based state container for respiration rate variables.
 *
 * This class manages device-side memory for the potential maintenance
 * respiration rate array. It is allocated once per simulation and updated
 * each day by the Respiration::calc_rates() kernel.
 */
class RespirationState {
public:

    /* ------------------------------------------------------------------- *
     * Grid dimensions
     * ------------------------------------------------------------------- */

    /** @brief Number of computational cells in the grid. */
    int nCells;

    /* ------------------------------------------------------------------- *
     * Rate Variables (output)
     * ------------------------------------------------------------------- */

    /**
     * @brief Potential maintenance respiration rate [kg CH2O / ha / d].
     *
     * This is the total carbohydrate expenditure required to maintain all
     * living crop tissue (roots, leaves, stems, storage organs) at the
     * current temperature, after applying the senescence reduction factor.
     *
     * Physical meaning:
     *   - Represents the biochemical cost of protein turnover, ion
     *     gradient maintenance, and cellular repair processes
     *   - Scales with living biomass and temperature
     *   - Must be satisfied before any assimilates are allocated to growth
     *
     * Typical ranges (wheat, temperate climate):
     *   Early growth (DVS < 0.5):  5 - 50 kg CH2O/ha/d
     *   Peak growth  (DVS ~ 1.0):  50 - 200 kg CH2O/ha/d
     *   Late growth  (DVS > 1.5):  10 - 80 kg CH2O/ha/d (reduced by senescence)
     */
    realArr PMRES;

    /* ------------------------------------------------------------------- *
     * Memory Allocation
     * ------------------------------------------------------------------- */

    /**
     * @brief Allocate device memory for all respiration state arrays.
     *
     * @param n  Number of computational cells to allocate for.
     *
     * This method must be called once before the first call to
     * Respiration::calc_rates(). It allocates a Kokkos View of size n
     * on the default device memory space.
     */
    void allocate(int n) {
        nCells = n;
        PMRES = realArr("PMRES", n);
    }
};

#endif
