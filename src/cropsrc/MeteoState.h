/* -*- mode: c++; c-default-style: "linux" -*- */
/**
 * @file MeteoState.h
 * @brief Meteorological forcing data storage for the WOFOST crop growth model.
 *
 * This file defines the MeteoState class, which stores all meteorological
 * time-series data required to drive the WOFOST crop simulation. The data
 * includes site metadata (location, elevation), simulation timing information,
 * and daily weather variables stored as Kokkos Views for GPU-compatible
 * parallel access.
 *
 * Data is organized as 1D arrays indexed by time step (day number), where
 * each element corresponds to one day in the simulation period. The class
 * provides a utility method (get_index_at_time) to map simulation time in
 * seconds to the corresponding array index.
 *
 * Required input variables (read from external weather files):
 *   - time  : Time offset from simulation start [s], step = 86400 s (1 day)
 *   - irrad : Daily total global solar radiation [J m^-2 d^-1]
 *   - tmin  : Daily minimum air temperature [C]
 *   - tmax  : Daily maximum air temperature [C]
 *   - vap   : Daily mean vapor pressure [hPa]
 *   - wind  : Daily mean wind speed [m s^-1]
 *   - rain  : Daily total precipitation [mm d^-1]
 *
 * Derived variables (computed during initialization from the above inputs):
 *   - e0    : Potential evaporation from open water surface [cm d^-1]
 *   - es0   : Potential evaporation from bare soil surface [cm d^-1]
 *   - et0   : Potential evapotranspiration from reference crop [cm d^-1]
 *
 * Site metadata:
 *   - lat, lon : Geographic coordinates [degrees]
 *   - elev     : Elevation above sea level [m]
 *   - angstA, angstB : Angstrom formula coefficients for radiation estimation
 *
 * @note Units for E0, ES0, ET0 are in [cm d^-1] following WOFOST convention,
 *       not [mm d^-1]. The conversion is handled during MeteoInit.
 *
 * @see ../define.h for the realArr typedef (Kokkos View<real*>)
 */

#ifndef _METEO_STATE_H_
#define _METEO_STATE_H_

#include "../define.h"

class MeteoState {

public:

    /*==================================================================*
     *  SITE METADATA                                                    *
     *==================================================================*/

    real lat;    ///< Geographic latitude [degrees], positive = North. Range: -90 to +90.
    real lon;    ///< Geographic longitude [degrees]. Range: -180 to +180.
    real elev;   ///< Elevation above sea level [m]. Affects air pressure and ET calculations.

    /**
     * @brief Angstrom formula coefficients for estimating solar radiation
     *        from sunshine duration.
     *
     * The Angstrom-Prescott formula: AVRAD = (angstA + angstB * n/N) * ANGOT
     * where n = actual sunshine hours, N = maximum possible sunshine hours,
     * and ANGOT = extraterrestrial radiation.
     *
     * Typical ranges:
     *   - angstA: 0.15-0.30 (intercept, overcast-day transmission)
     *   - angstB: 0.40-0.60 (slope, clear-sky increment)
     */
    real angstA;
    real angstB;

    /*==================================================================*
     *  SIMULATION TIMING                                                *
     *==================================================================*/

    int start_year;   ///< Year of simulation start (e.g., 2020)
    int start_month;  ///< Month of simulation start (1-12)
    int start_day;    ///< Day of month of simulation start (1-31)
    int start_doy;    ///< Day of year of simulation start (1-365/366)

    /*==================================================================*
     *  DATA DIMENSIONS                                                  *
     *==================================================================*/

    size_t num_records; ///< Total number of daily records in the time series [-]

    /*==================================================================*
     *  DAILY WEATHER TIME SERIES (Kokkos Views)                         *
     *==================================================================*/

    /**
     * @name Primary Meteorological Variables
     * Daily time-series data, one value per day, stored as Kokkos 1D Views.
     * Array index corresponds to day number (0-based).
     * @{ */

    realArr time;   ///< Time offset from simulation start [s]. Step = 86400 s (1 day).
                    ///< Example: day 0 = 0.0, day 1 = 86400.0, day 2 = 172800.0, ...

    realArr irrad;  ///< Daily total global solar radiation [J m^-2 d^-1].
                    ///< Typical range: 0 (polar night) to ~40e6 (tropical clear sky).

    realArr tmin;   ///< Daily minimum air temperature [C].
                    ///< Measured at screen height (~2 m). Typical range: -40 to +30.

    realArr tmax;   ///< Daily maximum air temperature [C].
                    ///< Measured at screen height (~2 m). Typical range: -20 to +55.

    realArr vap;    ///< Daily mean vapor pressure [hPa].
                    ///< Indicates atmospheric humidity. Typical range: 1 (cold/arid) to 40 (tropical/humid).

    realArr wind;   ///< Daily mean wind speed at 10 m height [m s^-1].
                    ///< Typical range: 0.5 (sheltered) to 15 (exposed/coastal).

    realArr rain;   ///< Daily total precipitation [mm d^-1].
                    ///< Water equivalent. Typical range: 0 (dry) to 100+ (extreme events).

    /** @} */

    /*==================================================================*
     *  DERIVED EVAPOTRANSPIRATION VARIABLES                             *
     *==================================================================*/

    /**
     * @name Computed Evapotranspiration Variables
     * These are derived from the primary meteorological variables using
     * Penman and Penman-Monteith equations during MeteoInit.
     * Units follow WOFOST convention: [cm d^-1], not [mm d^-1].
     * @{ */

    realArr e0;     ///< Potential evaporation from open water surface [cm d^-1].
                    ///< Computed using the Penman equation. Represents water loss
                    ///< from a free water surface (lake, pan). Typical range: 0.1-1.5.

    realArr es0;    ///< Potential evaporation from bare (unvegetated) soil [cm d^-1].
                    ///< Lower than E0 due to soil surface resistance.
                    ///< Typical range: 0.05-1.2.

    realArr et0;    ///< Potential evapotranspiration from a reference crop [cm d^-1].
                    ///< Computed using the Penman-Monteith equation for a hypothetical
                    ///< reference crop (well-watered grass, 12 cm height, albedo = 0.23).
                    ///< Typical range: 0.05-1.0.

    /** @} */

    /*==================================================================*
     *  UTILITY METHODS                                                  *
     *==================================================================*/

    /**
     * @brief Get the array index corresponding to a given simulation time.
     *
     * Converts simulation time (in seconds since start) to the corresponding
     * daily record index by integer division by 86400 (seconds per day).
     * Implements a step-function (piecewise-constant) interpolation: each
     * day's weather values are used for the full 24-hour period.
     *
     * @param[in] current_time  Simulation time since start [s].
     *                          Must be >= 0.
     * @return    Array index into the daily time-series [0, num_records-1].
     *
     * @note Assumes uniform daily time steps (86400 s). Does not handle
     *       sub-daily or irregular time steps.
     * @note Boundary handling: clamps the index to [0, num_records-1].
     */
    KOKKOS_INLINE_FUNCTION
    int get_index_at_time(real current_time) const {
        // Compute day index by integer division (floor)
        // Each day = 86400 seconds
        int idx = static_cast<int>(current_time / 86400.0);

        // Clamp to valid range
        if (idx < 0) idx = 0;
        // Note: num_records is size_t; cast to int for safe comparison
        if (idx >= static_cast<int>(num_records)) idx = static_cast<int>(num_records) - 1;

        return idx;
    }

    /*==================================================================*
     *  MEMORY ALLOCATION                                                *
     *==================================================================*/

    /**
     * @brief Allocate memory for all meteorological time-series arrays.
     *
     * Creates Kokkos Views on the default memory space (device or host,
     * depending on the Kokkos configuration). All arrays have the same
     * length n (number of daily records).
     *
     * @param[in] n  Number of daily records to allocate [size_t].
     *               Must match the length of the input weather data.
     *
     * Arrays are labeled with descriptive names for debugging and profiling.
     */
    void allocate(size_t n) {
        num_records = n;

        // Primary meteorological variables
        time  = realArr("meteo_time", n);    // [s]
        irrad = realArr("meteo_irrad", n);   // [J m^-2 d^-1]
        tmin  = realArr("meteo_tmin", n);    // [C]
        tmax  = realArr("meteo_tmax", n);    // [C]
        vap   = realArr("meteo_vap", n);     // [hPa]
        wind  = realArr("meteo_wind", n);    // [m s^-1]
        rain  = realArr("meteo_rain", n);    // [mm d^-1]

        // Derived evapotranspiration variables
        e0    = realArr("meteo_e0", n);      // [cm d^-1]
        es0   = realArr("meteo_es0", n);     // [cm d^-1]
        et0   = realArr("meteo_et0", n);     // [cm d^-1]
    }
};

#endif
