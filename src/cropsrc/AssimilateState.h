/* -*- mode: c++; c-default-style: "linux" -*- */
/**
 * @file AssimilateState.h
 * @brief Assimilation output state variables and astronomical data structure
 *        for the WOFOST crop growth model.
 *
 * This file defines two key data structures used by the assimilation module:
 *
 *   1. AssimilateState - Stores per-cell output rate variables from the
 *      photosynthesis/assimilation calculations. Contains the potential gross
 *      assimilation rate (PGASS) array, which holds the daily potential
 *      carbohydrate production for each computational cell.
 *
 *   2. AstroVars - A plain struct (not Kokkos View-based) that holds the
 *      astronomical variables computed by the ASTRO routine for a single
 *      location/day. These are scalar return values from the astro() function
 *      in Assimilation.h, used to pass solar geometry and radiation
 *      partitioning data between the astronomical and assimilation routines.
 *
 * Memory management:
 *   - AssimilateState uses Kokkos Views (realArr) for GPU-compatible
 *     parallel computation across all grid cells.
 *   - AstroVars is a lightweight scalar struct, created on the stack within
 *     Kokkos kernels for each cell's astronomical calculations.
 *
 * @see Assimilation.h for the routines that populate these structures
 * @see ../define.h    for the realArr typedef (Kokkos View<real*>)
 */

#ifndef _ASSIMILATE_STATE_H_
#define _ASSIMILATE_STATE_H_

#include "../define.h"

/*======================================================================*
 *  ASSIMILATION OUTPUT STATE                                            *
 *======================================================================*/

/**
 * @class AssimilateState
 * @brief Stores output rate variables from the crop assimilation (photosynthesis) module.
 *
 * Each member array is indexed by cell number (0 to nCells-1) and stores
 * the computed rate values for one simulation time step (typically one day).
 *
 * Array members:
 *   - PGASS: Potential gross canopy assimilation rate [kg CH2O ha^-1 d^-1].
 *     This is the total daily carbohydrate production assuming no water stress,
 *     nutrient limitations, or pest damage.
 *     Typical range: 0 (no canopy/night) to ~400 (full canopy, clear sky, warm).
 *     This value is later reduced by water stress, respiration, and partitioning
 *     factors to obtain actual growth rates.
 */
class AssimilateState {
public:
    int nCells;  ///< Number of computational cells in the simulation domain

    /* ---- Output Rate Variables ---- */

    /**
     * @brief Potential gross canopy assimilation rate [kg CH2O ha^-1 d^-1].
     *
     * Computed by the Assimilation::calc_rates() function. Represents the
     * maximum daily carbohydrate production by the crop canopy, integrating
     * the leaf-level photosynthesis model over canopy depth and daylight hours.
     * Units are in carbohydrate (CH2O) equivalents, converted from CO2 using
     * the molecular weight ratio 30/44.
     */
    realArr PGASS;

    /* ---- Memory Allocation ---- */

    /**
     * @brief Allocate memory for all state arrays.
     *
     * @param n  Number of computational cells
     *
     * Allocates Kokkos Views on the default memory space (device or host,
     * depending on the Kokkos configuration). All arrays are initialized
     * to zero by Kokkos.
     */
    void allocate(int n) {
        nCells = n;
        PGASS = realArr("PGASS", n);
    }
};

/*======================================================================*
 *  ASTRONOMICAL VARIABLES STRUCTURE                                     *
 *======================================================================*/

/**
 * @struct AstroVars
 * @brief Astronomical variables returned by the ASTRO routine for a single
 *        day and location.
 *
 * This struct is used as a lightweight return value from the
 * Assimilation::astro() function. It is created on the stack within
 * Kokkos parallel kernels, avoiding dynamic memory allocation on the GPU.
 *
 * All members are scalar values (not arrays) for a single cell/day.
 *
 * Members:
 *   - DAYL   : Astronomical daylength, sun above horizon [h]. Range: 0-24.
 *              Used as the integration period for daily assimilation.
 *
 *   - DAYLP  : Photoperiodic daylength, sun above -4 degrees [h]. Range: 0-24.
 *              The -4 degree angle corresponds to civil twilight, used in
 *              photoperiod-sensitive phenology calculations (e.g., flowering).
 *
 *   - SINLD  : sin(latitude) * sin(declination) [-].
 *              Component of the solar elevation equation due to latitude and
 *              season. Range: -0.4 to +0.4 approximately.
 *
 *   - COSLD  : cos(latitude) * cos(declination) [-].
 *              Component of the solar elevation equation due to the diurnal
 *              cycle. Range: 0.0 to 1.0. Always non-negative.
 *
 *   - DIFPP  : Diffuse photosynthetically active radiation (PAR) at canopy top
 *              [J m^-2 s^-1]. This is the maximum diffuse radiation available
 *              for canopy photosynthesis. Typical range: 0-150.
 *
 *   - ATMTR  : Atmospheric transmission coefficient [-].
 *              Ratio of measured global radiation to extraterrestrial radiation.
 *              Range: 0 (heavy overcast) to ~0.75 (clear sky).
 *
 *   - DSINBE : Daily integral of sin(beta) * (1 + 0.4*sin(beta)) over daylight
 *              hours [s]. Used to distribute daily radiation into instantaneous
 *              values for the diurnal Gaussian integration in TOTASS7.
 *              Typical range: 0 to ~120000 s.
 *
 *   - ANGOT  : Angot extraterrestrial radiation [J m^-2 d^-1].
 *              The theoretical maximum radiation at the top of the atmosphere
 *              for the given day and latitude. Typical range: 0 to ~40e6.
 *              Used as the denominator for computing atmospheric transmission.
 */
struct AstroVars {
    real DAYL;    ///< Astronomical daylength (sun above horizon) [h], range: 0-24
    real DAYLP;   ///< Photoperiodic daylength (sun above -4 deg) [h], range: 0-24
    real SINLD;   ///< sin(lat) * sin(DEC) [-], range: ~-0.4 to +0.4
    real COSLD;   ///< cos(lat) * cos(DEC) [-], range: 0.0 to 1.0
    real DIFPP;   ///< Diffuse PAR at canopy top [J m^-2 s^-1], range: 0-150
    real ATMTR;   ///< Atmospheric transmission coefficient [-], range: 0-0.75
    real DSINBE;  ///< Daily integral of sin(beta)*(1+0.4*sin(beta)) [s]
    real ANGOT;   ///< Angot extraterrestrial radiation [J m^-2 d^-1]
};

#endif
