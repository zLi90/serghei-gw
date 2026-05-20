/* -*- mode: c++; c-default-style: "linux" -*- */
/**
 * @file EvapotranspirationState.h
 * @brief Evapotranspiration state arrays for the WOFOST crop growth model.
 *
 * This file defines the EvapotranspirationState class, which stores per-cell
 * evapotranspiration rate variables computed by the Evapotranspiration module.
 * Each member is a Kokkos 1D View (realArr) indexed by cell number (0 to nCells-1).
 *
 * The state variables represent the daily water balance components for the
 * soil-crop-atmosphere system:
 *
 *   - Reference ET values (E0, ES0, ET0) copied from MeteoState for convenience
 *   - Maximum (unstressed) rates accounting for canopy effects (EVS, TRAMX)
 *   - Actual rates reduced by water/oxygen stress (TRA)
 *   - Stress reduction factors used by other modules (RFTRA)
 *   - Root-zone soil moisture snapshot used in the stress calculations (RZSM)
 *
 * All rate variables have units of [cm d^-1] following WOFOST convention.
 * Conversion to [mm d^-1] requires multiplication by 10.
 *
 * Data flow:
 *   MeteoState (E0, ES0, ET0) -> Evapotranspiration::calc_rates()
 *     -> EvapotranspirationState (all fields below)
 *     -> used by Respiration, Growth, and WaterBalance modules
 *
 * @see Evapotranspiration.h  for the calculations that populate these arrays
 * @see ../define.h           for the realArr typedef (Kokkos View<real*>)
 */

#ifndef _EVAPOTRANSPIRATION_STATE_H_
#define _EVAPOTRANSPIRATION_STATE_H_

#include "../define.h"

/**
 * @class EvapotranspirationState
 * @brief Per-cell evapotranspiration rate variables for one simulation time step.
 *
 * All arrays are allocated via the allocate() method and populated by
 * Evapotranspiration::calc_rates() on each simulation day.
 */
class EvapotranspirationState {
public:
    int nCells;  ///< Number of computational cells in the simulation domain

    /*==================================================================*
     *  RATE VARIABLES (per day)                                         *
     *==================================================================*/

    /**
     * @name Actual and Potential Evapotranspiration Rates
     * All rates are in [cm d^-1] (WOFOST convention).
     * To convert to [mm d^-1], multiply by 10.
     * @{ */

    /**
     * @brief Actual crop transpiration rate [cm d^-1].
     *
     * The water actually transpired by the crop, after accounting for water
     * stress and oxygen stress. Computed as TRAMX * RFTRA.
     * Typical range: 0 (stressed/dormant) to ~0.8 (well-watered, full canopy).
     * This value is used by the water balance and growth respiration modules.
     */
    realArr TRA;

    /**
     * @brief Maximum (potential) crop transpiration rate [cm d^-1].
     *
     * The transpiration rate assuming no water or oxygen stress, accounting
     * only for canopy cover (LAI) effects on energy partitioning.
     * Computed as CFET * ET0 * (1 - exp(-KGLOB * LAI)).
     * Typical range: 0 (bare soil) to ~0.8 (full canopy, high ET0).
     */
    realArr TRAMX;

    /**
     * @brief Actual evaporation rate from the soil surface [cm d^-1].
     *
     * Maximum soil evaporation below the canopy, reduced by canopy shading
     * via the exponential attenuation of radiation through the canopy.
     * Computed as ES0 * exp(-KGLOB * LAI).
     * Typical range: 0 (full canopy) to ~0.5 (bare soil, high ET0).
     */
    realArr EVS;

    /**
     * @brief Potential evaporation from open water surface [cm d^-1].
     *
     * Copied from MeteoState for convenient per-cell access. Computed using
     * the Penman equation for a free water surface.
     * Typical range: 0.1-1.5 depending on climate.
     */
    realArr E0;

    /**
     * @brief Potential evaporation from bare soil surface [cm d^-1].
     *
     * Copied from MeteoState. Computed using the Penman equation with soil
     * surface resistance adjustments.
     * Typical range: 0.05-1.2 depending on climate.
     */
    realArr ES0;

    /**
     * @brief Potential evapotranspiration from reference crop [cm d^-1].
     *
     * Copied from MeteoState. Computed using the Penman-Monteith equation
     * for the FAO reference crop (well-watered grass, 12 cm height, albedo = 0.23).
     * Typical range: 0.05-1.0 depending on climate.
     */
    realArr ET0;

    /** @} */

    /*==================================================================*
     *  STRESS REDUCTION FACTORS                                         *
     *==================================================================*/

    /**
     * @brief Combined reduction factor for transpiration due to water/oxygen stress [-].
     *
     * Dimensionless factor in the range [0, 1]:
     *   - 1.0: no stress, crop transpires at the maximum rate (TRAMX)
     *   - 0.0: complete stress, transpiration ceases
     *
     * Computed as: RFTRA = RFOS * RFWS, where:
     *   RFOS = oxygen stress factor (waterlogging)
     *   RFWS = water stress factor (drought)
     *
     * This factor is used by the growth module to reduce assimilation rates
     * under water-limited conditions.
     */
    realArr RFTRA;

    /*==================================================================*
     *  SOIL MOISTURE STATE                                              *
     *==================================================================*/

    /**
     * @brief Root-zone soil moisture content [cm^3 cm^-3].
     *
     * Volumetric water content averaged over the root zone depth, provided
     * by the SERGHEI hydrological model. Used in water stress calculations
     * to determine RFWS and RFOS.
     * Typical range: SMW (wilting point, ~0.05-0.15) to SM0 (saturation, ~0.35-0.55).
     */
    realArr RZSM;

    /*==================================================================*
     *  MEMORY ALLOCATION                                                *
     *==================================================================*/

    /**
     * @brief Allocate memory for all evapotranspiration state arrays.
     *
     * Creates Kokkos Views on the default memory space (device or host,
     * depending on the Kokkos configuration). All arrays have the same
     * length n (number of computational cells).
     *
     * @param[in] n  Number of computational cells
     *
     * Arrays are labeled with descriptive names for debugging and profiling
     * with Kokkos tools.
     */
    void allocate(int n) {
        nCells = n;
        TRA   = realArr("TRA", n);    // Actual transpiration [cm d^-1]
        TRAMX = realArr("TRAMX", n);  // Maximum transpiration [cm d^-1]
        EVS   = realArr("EVS", n);    // Maximum soil evaporation [cm d^-1]
        E0    = realArr("E0", n);     // Potential open water evaporation [cm d^-1]
        ES0   = realArr("ES0", n);    // Potential bare soil evaporation [cm d^-1]
        ET0   = realArr("ET0", n);    // Potential reference crop ET [cm d^-1]
        RFTRA = realArr("RFTRA", n);  // Transpiration reduction factor [-]
        RZSM  = realArr("RZSM", n);   // Root-zone soil moisture [cm^3 cm^-3]
    }
};

#endif
