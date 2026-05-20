/* -*- mode: c++; c-default-style: "linux" -*- */
/**
 * @file Phenology.h
 * @brief Crop phenology development module for the WOFOST crop growth model.
 *
 * This file implements the phenological development of a crop through its life
 * cycle, calculating the daily development rate (DVR) based on temperature,
 * vernalisation (for winter crops), and photoperiod sensitivity. The crop
 * progresses through four growth stages:
 *
 *   - EMERGING   : Sowing to emergence (DVS: -0.1 -> 0.0)
 *   - VEGETATIVE : Emergence to anthesis (DVS: 0.0 -> 1.0)
 *   - REPRODUCTIVE : Anthesis to maturity (DVS: 1.0 -> DVSEND, typically 2.0)
 *   - MATURE     : Post-maturity; development stops
 *
 * Key concepts:
 *   - DVS (Development Stage): A dimensionless scalar tracking crop maturity.
 *     Typical range: -0.1 (sowing) to 2.0 (maturity).
 *   - TSUM (Temperature Sum): Accumulated effective temperature [C d] driving
 *     development from one stage to the next.
 *   - Vernalisation: A cold-temperature requirement that must be satisfied
 *     before rapid vegetative development can proceed (e.g., winter wheat).
 *     Controlled by parameters VERNBASE, VERNSAT, and the VERNRTB table.
 *   - Photoperiod sensitivity: Daylength-dependent development rate reduction,
 *     applicable when IDSL >= 1. DLO (optimum daylength) and DLC (critical
 *     daylength) define the response curve.
 *
 * References:
 *   - WOFOST 7.1 User Guide, Alterra, Wageningen UR
 *   - Boogaard et al. (2014), WOFOST Control Centre 2.1
 *
 * @see CropState.h        for static crop parameters (temperature thresholds, TSUM values)
 * @see CropDynamicState.h for dynamic state arrays (DVS, STAGE, TSUM, VERN)
 * @see MeteoState.h       for meteorological forcing data (TMIN, TMAX)
 * @see Afgen.h            for table-lookup interpolation utilities
 */

#ifndef _PHENOLOGY_H_
#define _PHENOLOGY_H_

#include "../define.h"
#include "CropState.h"        // Static crop parameters (temperature thresholds, TSUM targets)
#include "CropDynamicState.h" // Dynamic state arrays (current development stage, DVS, accumulated TSUM)
#include "MeteoState.h"       // Meteorological forcing (daily min/max temperature [C])
#include "Afgen.h"            // Piecewise linear interpolation for temperature-response tables
#include <cmath>

class Phenology {
public:

    /*====================================================================*
     *  INITIALIZATION                                                    *
     *====================================================================*/

    /**
     * @brief Initialize crop phenological state at the start of simulation.
     *
     * Sets the initial development stage (DVS) and growth stage (STAGE) for
     * each computational cell, based on the crop parameter DVSI (initial DVS
     * value). Also initializes accumulated temperature sums and vernalisation
     * state to zero.
     *
     * @param[out] s  Dynamic crop state to initialize
     * @param[in]  p  Static crop parameters (uses p.DVSI for initial DVS)
     *
     * Initialization logic:
     *   - If DVSI < 0: crop starts from sowing (STAGE = EMERGING, DVS = DVSI,
     *     typically -0.1). The crop must emerge before vegetative growth begins.
     *   - If DVSI >= 0: crop is assumed already emerged (STAGE = VEGETATIVE,
     *     DVS = DVSI, typically 0.0). Used for re-starting simulations or when
     *     transplanting established seedlings.
     */
    void initialize(CropDynamicState &s, const CropState &p) {
        Kokkos::parallel_for("Pheno_Init", s.nCells, KOKKOS_LAMBDA(const int i) {
            if (p.p.DVSI < 0.0) {
                // Crop starts from sowing: must pass through emergence phase
                s.STAGE(i) = CropStage::EMERGING;  // Emergence phase (DVS: -0.1 -> 0.0)
                s.DVS(i) = p.p.DVSI;               // Initial DVS, e.g. -0.1 at sowing
            } else {
                // Crop is already emerged: begin vegetative growth directly
                s.STAGE(i) = CropStage::VEGETATIVE; // Vegetative phase (DVS: 0.0 -> 1.0)
                s.DVS(i) = p.p.DVSI;               // Initial DVS, e.g. 0.0 at emergence
            }

            // Reset accumulated state variables
            s.TSUM(i) = 0.0;    // Accumulated effective temperature for vegetative/reproductive phase [C d]
            s.TSUME(i) = 0.0;   // Accumulated effective temperature for emergence phase [C d] (reserved for future use)
            s.VERN(i) = 0.0;    // Accumulated vernalisation units [dimensionless, typically degree-days]
            s.ISVERNALISED(i) = 0; // Vernalisation completion flag: 0 = not vernalised, 1 = vernalised
        });
    }

    /*====================================================================*
     *  RATE CALCULATION                                                  *
     *====================================================================*/

    /**
     * @brief Calculate daily phenological development rates.
     *
     * Computes the daily development rate (DVR) [d^-1] for each cell, which
     * governs how fast the crop advances through its life cycle. The rate is
     * determined by:
     *   1. Mean daily air temperature (via the DTSMTB lookup table)
     *   2. Vernalisation status (winter crops only, IDSL >= 2)
     *   3. Photoperiod sensitivity (daylength crops, IDSL >= 1)
     *
     * @param[in,out] s  Dynamic crop state (reads DVS, STAGE; writes DVR, DTSUM, VERNR)
     * @param[in]     p  Static crop parameters (IDSL, TBASEM, TSUMEM, TSUM1, TSUM2, etc.)
     * @param[in]     m  Meteorological forcing (TMIN, TMAX)
     * @param[in]     drv_idx          Index into meteorological time-series arrays
     * @param[in]     current_day_of_year  Day of year (1-365), used for daylength calculation
     *
     * Temperature processing:
     *   - TEMP = (TMIN + TMAX) / 2.0 : daily mean air temperature [C]
     *     Typical range: -10 to +40 C depending on climate and season.
     *
     * Vernalisation (IDSL >= 2):
     *   - Applicable only during the VEGETATIVE stage before vernalisation
     *     is complete.
     *   - The vernalisation rate VERNR is looked up from VERNRTB(TEMP).
     *   - VERNFAC scales the development rate from 0 (no vernalisation yet)
     *     to 1 (fully vernalised), based on accumulated vernalisation (VERN)
     *     relative to VERNBASE and VERNSAT thresholds.
     *   - If DVS >= VERNDVS (0.3), vernalisation is forced complete.
     *
     * Photoperiod (IDSL >= 1):
     *   - DVRED = (DAYLP - DLC) / (DLO - DLC), clamped to [0, 1].
     *   - DLO: optimum daylength for development [h], typically 12-18 h.
     *   - DLC: critical daylength below which development is suppressed [h].
     *   - DAYLP: photoperiodic daylength [h], calculated with sun angle = -4 degrees.
     */
    void calc_rates(CropDynamicState &s, const CropState &p, const MeteoState &m, int drv_idx, int current_day_of_year) {
        // Retrieve phenological lookup tables
        auto DTSMTB = p.tables.at("DTSMTB");    // Temperature -> effective temperature sum rate [C d d^-1]
        auto VERNRTB = p.hasTable("VERNRTB") ? p.tables.at("VERNRTB") : realArr(); // Temperature -> vernalisation rate table (optional)

        Kokkos::parallel_for("Pheno_Rates", s.nCells, KOKKOS_LAMBDA(const int i) {

            /* ---- Step 1: Compute daily mean air temperature ---- */
            real TMIN = m.tmin(drv_idx);           // Daily minimum temperature [C]
            real TMAX = m.tmax(drv_idx);           // Daily maximum temperature [C]
            real TEMP = (TMIN + TMAX) / 2.0;      // Daily mean temperature [C]

            /* ---- Step 2: Vernalisation factor (VERNFAC) ---- */
            // Vernalisation is required for winter crops (IDSL >= 2) and applies
            // only during the vegetative stage before the vernalisation requirement
            // has been satisfied.
            real VERNFAC = 1.0; // Vernalisation factor: 0 = no development, 1 = full rate [-]
            if (p.p.IDSL >= 2) { // IDSL >= 2 indicates the crop has a vernalisation requirement
                if (s.STAGE(i) == CropStage::VEGETATIVE && s.ISVERNALISED(i) == 0) {
                    // Critical DVS threshold beyond which vernalisation is forced complete
                    real VERNDVS = 0.3; // Typical value: 0.3 (early vegetative stage)

                    if (s.DVS(i) < VERNDVS) {
                        // Look up daily vernalisation rate from temperature response table
                        real rate = Afgen::lookup(VERNRTB, TEMP);
                        s.VERNR(i) = rate; // Store daily vernalisation rate [vern units d^-1]

                        // Compute vernalisation factor from accumulated vernalisation.
                        // VERNFAC increases linearly from 0 at VERNBASE to 1 at VERNSAT.
                        real req_base = p.p.VERNBASE; // Minimum vernalisation for any development [vern units]
                        real req_sat = p.p.VERNSAT;   // Vernalisation needed for full development [vern units]
                        real r = (s.VERN(i) - req_base) / (req_sat - req_base);
                        r = std::max(0.0, std::min(1.0, r)); // Clamp to [0, 1]
                        VERNFAC = r;
                    } else {
                        // DVS has exceeded the critical threshold: force vernalisation complete
                        // to prevent the crop from being permanently stuck at low VERNFAC
                        s.VERNR(i) = 0.0;
                        VERNFAC = 1.0;
                    }
                } else {
                    // Not in the correct stage, or already vernalised
                    s.VERNR(i) = 0.0;
                    VERNFAC = 1.0;
                }
            }

            /* ---- Step 3: Photoperiod reduction factor (DVRED) ---- */
            // Photoperiod sensitivity is applicable when IDSL >= 1.
            // Short-day crops (e.g., soybean) and long-day crops (e.g., wheat)
            // differ in their DLO/DLC parameter values.
            real DVRED = 1.0; // Photoperiod reduction factor: 0 = no development, 1 = optimum [-]
            if (p.p.IDSL >= 1) {
                // Calculate photoperiodic daylength (sun angle = -4 degrees for civil twilight)
                real DAYLP = this->calculate_daylength(current_day_of_year, m.lat);

                // Photoperiod reduction: linear interpolation between critical and optimum daylength
                // DLO: optimum daylength for maximum development rate [h]
                // DLC: critical daylength below which development is suppressed [h]
                real DLO = p.p.DLO;
                real DLC = p.p.DLC;

                real ratio = (DAYLP - DLC) / (DLO - DLC);

                // Clamp ratio to [0, 1]
                if (ratio < 0.0) ratio = 0.0;
                if (ratio > 1.0) ratio = 1.0;

                DVRED = ratio;
            }

            /* ---- Step 4: Compute development rate (DVR) by growth stage ---- */
            // DVR [d^-1] determines how many DVS units are accumulated per day.
            // The calculation differs by stage:
            if (s.STAGE(i) == CropStage::EMERGING) {
                // Emergence phase (DVS: -0.1 -> 0.0)
                // Effective temperature above the base temperature for emergence
                real dtsume = TEMP - p.p.TBASEM; // TBASEM: base temperature for emergence [C], typically 0-10 C
                if (dtsume < 0.0) dtsume = 0.0;  // No emergence progress below base temperature
                // DVR = 0.1 * effective_temperature / TSUMEM
                // The 0.1 factor arises because DVS must increase by exactly 0.1 (from -0.1 to 0.0)
                // over the total temperature sum TSUMEM [C d].
                s.DVR(i) = 0.1 * dtsume / p.p.TSUMEM; // DVR [d^-1]; TSUMEM: temperature sum for emergence [C d]
            }
            else if (s.STAGE(i) == CropStage::VEGETATIVE) {
                // Vegetative phase (DVS: 0.0 -> 1.0)
                // Effective daily temperature sum from lookup table, modified by vernalisation and photoperiod
                s.DTSUM(i) = Afgen::lookup(DTSMTB, TEMP) * VERNFAC * DVRED; // Effective temperature sum rate [C d d^-1]
                s.DVR(i) = s.DTSUM(i) / p.p.TSUM1; // DVR [d^-1]; TSUM1: temperature sum from emergence to anthesis [C d]
            }
            else if (s.STAGE(i) == CropStage::REPRODUCTIVE) {
                // Reproductive phase (DVS: 1.0 -> DVSEND, typically 2.0)
                // No vernalisation or photoperiod correction in the reproductive phase
                s.DTSUM(i) = Afgen::lookup(DTSMTB, TEMP); // Effective temperature sum rate [C d d^-1]
                s.DVR(i) = s.DTSUM(i) / p.p.TSUM2; // DVR [d^-1]; TSUM2: temperature sum from anthesis to maturity [C d]
            }
            else if (s.STAGE(i) == CropStage::MATURE) {
                // Maturity: no further phenological development
                s.DTSUM(i) = 0.0;
                s.DVR(i) = 0.0;
            }
            // Error handling: undefined stage
            else {
                printf("Error: no stage defined in phenology calculation.\n");
            }

        });
    }

    /*====================================================================*
     *  STATE INTEGRATION                                                 *
     *====================================================================*/

    /**
     * @brief Integrate (advance) phenological state by one time step.
     *
     * Updates the accumulated vernalisation units, development stage (DVS),
     * accumulated temperature sums (TSUM), and triggers stage transitions
     * when threshold DVS values are reached.
     *
     * @param[in,out] s  Dynamic crop state (updates DVS, TSUM, VERN, STAGE, ISVERNALISED)
     * @param[in]     p  Static crop parameters (uses VERNSAT, DVSEND)
     * @param[in]     delt  Time step duration [d], default = 1.0 (daily)
     *
     * Stage transitions:
     *   - EMERGING -> VEGETATIVE : when DVS reaches 0.0
     *   - VEGETATIVE -> REPRODUCTIVE : when DVS reaches 1.0 (anthesis)
     *   - REPRODUCTIVE -> MATURE : when DVS reaches DVSEND (typically 2.0)
     */
    void integrate(CropDynamicState &s, const CropState &p, real delt = 1.0) {
        Kokkos::parallel_for("Pheno_Integrate", s.nCells, KOKKOS_LAMBDA(const int i) {

            /* ---- Step 1: Integrate vernalisation (IDSL >= 2 only) ---- */
            if (p.p.IDSL >= 2) {
                if (s.STAGE(i) == CropStage::VEGETATIVE) {
                    s.VERN(i) += s.VERNR(i) * delt; // Accumulate vernalisation units [vern units]

                    // Check if vernalisation requirement has been satisfied
                    if (s.VERN(i) >= p.p.VERNSAT) {
                        // Accumulated vernalisation exceeds saturation threshold
                        s.ISVERNALISED(i) = 1; // Mark crop as fully vernalised
                    } else if (s.DVS(i) >= 0.3 && s.ISVERNALISED(i) == 0) {
                        // Failsafe: if DVS has reached the critical threshold (0.3)
                        // but vernalisation is still incomplete, force it complete
                        // to prevent the crop from being permanently development-suppressed
                        s.ISVERNALISED(i) = 1;
                    }
                }
            }

            /* ---- Step 2: Integrate development stage and trigger transitions ---- */
            if (s.STAGE(i) == CropStage::EMERGING) {
                // Advance DVS during emergence phase (DVS: -0.1 -> 0.0)
                s.DVS(i) += s.DVR(i) * delt;
                if (s.DVS(i) >= 0.0) {
                    // Emergence complete: transition to vegetative stage
                    s.STAGE(i) = CropStage::VEGETATIVE;
                    s.DVS(i) = 0.0;   // Reset DVS to start of vegetative phase
                    s.TSUM(i) = 0.0;  // Reset accumulated temperature sum
                }
            }
            else if (s.STAGE(i) == CropStage::VEGETATIVE) {
                // Advance DVS during vegetative phase (DVS: 0.0 -> 1.0)
                s.TSUM(i) += s.DTSUM(i) * delt; // Accumulate effective temperature [C d]
                s.DVS(i) += s.DVR(i) * delt;    // Advance development stage
                if (s.DVS(i) >= 1.0) {
                    // Anthesis reached: transition to reproductive stage
                    s.STAGE(i) = CropStage::REPRODUCTIVE;
                    s.DVS(i) = 1.0; // Clamp DVS to exactly 1.0 at anthesis
                }
            }
            else if (s.STAGE(i) == CropStage::REPRODUCTIVE) {
                // Advance DVS during reproductive phase (DVS: 1.0 -> DVSEND)
                s.TSUM(i) += s.DTSUM(i) * delt; // Accumulate effective temperature [C d]
                s.DVS(i) += s.DVR(i) * delt;    // Advance development stage
                if (s.DVS(i) >= p.p.DVSEND) {
                    // Maturity reached: transition to mature stage
                    s.STAGE(i) = CropStage::MATURE;
                    s.DVS(i) = p.p.DVSEND; // Clamp DVS to end value (typically 2.0)
                }
            }
            else if (s.STAGE(i) == CropStage::MATURE) {
                // No further phenological changes after maturity
            }
            // Error handling: undefined stage
            else {
                printf("Error: no stage defined in phenology integration.\n");
            }
        });
    }

    /*====================================================================*
     *  ASTRONOMICAL UTILITIES                                           *
     *====================================================================*/

    /**
     * @brief Calculate photoperiodic daylength for a given day and latitude.
     *
     * Computes the length of the day (in hours) during which the sun is above
     * a specified elevation angle. By default, the angle is set to -4 degrees,
     * which corresponds to the threshold for civil twilight (used in WOFOST for
     * photoperiod-sensitive development calculations).
     *
     * The algorithm is derived from the WOFOST ASTRO.FOR routine and uses
     * standard astronomical formulae for solar declination and the hour angle
     * at sunrise/sunset.
     *
     * @param[in] day_of_year  Day of year (1 = Jan 1, 365 = Dec 31)
     * @param[in] latitude     Geographic latitude [degrees], positive = North
     * @param[in] angle        Sun elevation angle threshold [degrees],
     *                         default = -4.0 (civil twilight boundary)
     * @return    DAYLP        Photoperiodic daylength [h], range [0, 24]
     *
     * Physical basis:
     *   - Solar declination (DEC) is approximated assuming perihelion at
     *     day 10 (January 10).
     *   - The hour angle at sunrise/sunset is determined by the condition
     *     sin(sun_elevation) = sin(angle).
     *   - Daylength = (2/15) * arccos(AOB) * (180/pi), expressed as 12*(1 + 2*asin(AOB)/PI).
     *   - At high latitudes during summer, the sun may never set (DAYLP = 24);
     *     during winter, it may never rise (DAYLP = 0).
     */
    KOKKOS_INLINE_FUNCTION
    real calculate_daylength(int day_of_year, real latitude, real angle = -4.0) const {
        const real RAD = 0.017453292519943295; // Degrees to radians conversion factor (PI / 180.0)

        // Solar declination (DEC) [radians]
        // Approximation: DEC oscillates between +/- 23.45 degrees over the year
        // Offset of +10 days accounts for the fact that perihelion (closest
        // approach to the sun) occurs around January 3-5
        real sin_obliquity = sin(23.45 * RAD);              // sin(obliquity of ecliptic)
        real cos_pos = cos(2.0 * PI * (real(day_of_year) + 10.0) / 365.0);
        real DEC = -asin(sin_obliquity * cos_pos);          // Declination [radians], range [-23.45, +23.45] deg

        // Sine and cosine of the product of latitude and declination
        real SINLD = sin(RAD * latitude) * sin(DEC);  // sin(lat)*sin(DEC) [-]
        real COSLD = cos(RAD * latitude) * cos(DEC);  // cos(lat)*cos(DEC) [-]

        // Hour angle at sunrise/sunset for the given elevation angle
        // AOB = (-sin(angle) + sin(lat)*sin(DEC)) / (cos(lat)*cos(DEC))
        real AOB = (-sin(angle * RAD) + SINLD) / COSLD;

        // Compute daylength from the hour angle
        real DAYLP = 0.0;
        if (std::abs(AOB) <= 1.0) {
            // Normal case: sun rises and sets
            DAYLP = 12.0 * (1.0 + 2.0 * asin(AOB) / PI); // Daylength [h]
        } else if (AOB > 1.0) {
            // Polar day: sun never sets (midnight sun)
            DAYLP = 24.0;
        } else {
            // Polar night: sun never rises
            DAYLP = 0.0;
        }

        return DAYLP;
    }

};

#endif
