/* -*- mode: c++; c-default-style: "linux" -*- */
/**
 * @file Assimilation.h
 * @brief Photosynthesis and CO2 assimilation calculations for the WOFOST crop growth model.
 *
 * This file implements the three core routines of the WOFOST assimilation module:
 *
 *   1. ASTRO   - Astronomical calculations (daylength, solar declination,
 *                atmospheric transmission, diffuse radiation fraction).
 *   2. ASSIM7  - Instantaneous canopy CO2 assimilation using a three-point
 *                Gaussian integration over canopy depth (Spitters, 1986).
 *   3. TOTASS7 - Daily total gross canopy CO2 assimilation using a three-point
 *                Gaussian integration over the daylength period.
 *
 * The main entry point is calc_rates(), which orchestrates the full assimilation
 * pipeline: astronomical variables -> instantaneous assimilation -> daily total ->
 * conversion to carbohydrate equivalents (PGASS).
 *
 * Key physical concepts:
 *   - Gross CO2 assimilation (DTGA) is computed in [kg CO2 ha^-1 d^-1] and then
 *     converted to carbohydrate equivalents [kg CH2O ha^-1 d^-1] using the
 *     molecular weight ratio 30/44 (CH2O/CO2).
 *   - The canopy is treated as a multilayer system with separate sunlit and
 *     shaded leaf fractions for realistic light-response modeling.
 *   - A Gaussian quadrature scheme (3 points) is used for numerical integration
 *     over both canopy depth and diurnal time, providing a good balance between
 *     accuracy and computational cost.
 *
 * References:
 *   - Spitters, C.J.T. (1986). Separating the diffuse and direct component of
 *     global radiation and its implications for model canopy photosynthesis.
 *     Agricultural and Forest Meteorology, 38:217-229.
 *   - van Laar, H.H., et al. (1997). Simulation of Crop Growth: WOFOST 7.1.
 *     Ab-Tera, Wageningen.
 *
 * @see AssimilateState.h  for output state (PGASS) and AstroVars structure
 * @see CropState.h        for static crop parameters and AFGEN lookup tables
 * @see CropDynamicState.h for dynamic crop state (DVS, LAI)
 * @see MeteoState.h       for meteorological forcing (TMIN, TMAX, IRRAD)
 */

#ifndef _ASSIMILATION_H_
#define _ASSIMILATION_H_

#include "../define.h"
#include "CropState.h"
#include "CropDynamicState.h"
#include "AssimilateState.h"
#include "MeteoState.h"
#include "Afgen.h"
#include <cmath>
#include <cstdio>

class Assimilation
{

private:
    /** @brief Degrees-to-radians conversion constant. */
    static constexpr real RAD = 0.017453292519943295; // PI / 180.0

public:

    /*====================================================================*
     *  ASTRONOMICAL CALCULATIONS (ASTRO)                                 *
     *====================================================================*/

    /**
     * @brief Calculate astronomical daylength and radiation characteristics.
     *
     * Implements the WOFOST ASTRO routine, which computes solar geometry
     * variables needed for the assimilation calculations. These include
     * the actual and photoperiodic daylength, sine/cosine of latitude *
     * declination products, atmospheric transmission, and the diffuse
     * radiation component.
     *
     * @param[in] day_of_year  Day of year (1-365)
     * @param[in] latitude     Geographic latitude [degrees], positive = North.
     *                         Typical range: -60 to +60 for major agricultural regions.
     * @param[in] radiation    Daily total global solar radiation [J m^-2 d^-1].
     *                         Typical range: 0 (polar night) to 40e6 (tropical desert).
     * @param[in] debug        Enable debug output for the first cell
     * @return    AstroVars    Structure containing all computed astronomical variables:
     *   - DAYL   : Astronomical daylength (sun above horizon) [h]. Range: 0-24.
     *   - DAYLP  : Photoperiodic daylength (sun above -4 deg) [h]. Range: 0-24.
     *   - SINLD  : sin(latitude) * sin(declination) [-]. Range: -0.4 to +0.4.
     *   - COSLD  : cos(latitude) * cos(declination) [-]. Range: 0.0 to 1.0.
     *   - DIFPP  : Diffuse photosynthetically active radiation (PAR) [J m^-2 s^-1].
     *   - ATMTR  : Atmospheric transmission coefficient [-]. Range: 0-1.
     *   - DSINBE : Integral of sin(solar_elevation) * (1 + 0.4*sin(elevation))
     *              over the day [s]. Used for daily PAR distribution.
     *   - ANGOT  : Angot extraterrestrial radiation [J m^-2 d^-1].
     *              Typical range: 0 to ~40e6.
     *
     * Physical basis:
     *   - Solar declination varies annually between +/- 23.45 degrees.
     *   - Solar constant is ~1370 W m^-2, with a small annual oscillation
     *     due to Earth's elliptical orbit (eccentricity correction).
     *   - The diffuse fraction (FRDIF) of radiation increases under cloudy
     *     conditions (low atmospheric transmission) and is nearly zero
     *     under clear skies (high transmission).
     */
    KOKKOS_INLINE_FUNCTION
    AstroVars astro(int day_of_year, real latitude, real radiation, bool debug) const
    {
        AstroVars av;

        real LAT = latitude;            // Latitude [degrees]
        real IDAY = (real)day_of_year;  // Day of year [-]
        real AVRAD = radiation;         // Global solar radiation [J m^-2 d^-1]
        real ANGLE = -4.0;              // Sun elevation angle for photoperiodic daylength [degrees]

        /* ---- Solar declination and solar constant ---- */
        // Declination (DEC): angle between the sun's rays and the equatorial plane [radians]
        // Varies from -23.45 deg (southern solstice) to +23.45 deg (northern solstice)
        real sin_2345 = sin(23.45 * RAD);
        real cos_time = cos(2.0 * PI * (IDAY + 10.0) / 365.0);
        real DEC = -asin(sin_2345 * cos_time);

        // Solar constant with eccentricity correction [J m^-2 s^-1]
        // Base value: ~1370 W m^-2; annual variation amplitude: ~3.3%
        real SC = 1370.0 * (1.0 + 0.033 * cos(2.0 * PI * IDAY / 365.0));

        /* ---- Latitude-declination geometry ---- */
        real sin_lat = sin(RAD * LAT);
        real cos_lat = cos(RAD * LAT);
        real sin_dec = sin(DEC);
        real cos_dec = cos(DEC);

        // sin(lat)*sin(DEC): component of solar elevation due to latitude and season
        av.SINLD = sin_lat * sin_dec;
        // cos(lat)*cos(DEC): component of solar elevation due to diurnal cycle
        av.COSLD = cos_lat * cos_dec;

        /* ---- Daylength and integrated solar elevation ---- */
        // AOB: tangent of the hour angle at sunrise/sunset (sun above horizon, angle=0)
        real AOB = av.SINLD / av.COSLD;
        real DSINB = 0.0;   // Integrated sin(solar_elevation) over the day [s]

        if (std::abs(AOB) <= 1.0)
        {
            // Normal case: sun rises and sets
            av.DAYL = 12.0 * (1.0 + 2.0 * asin(AOB) / PI); // Astronomical daylength [h]
            real sqrt_1_AOB2 = sqrt(1.0 - AOB * AOB);
            // DSINB: integral of sin(beta) over daylight hours, where beta is solar elevation
            // Used to compute total daily extraterrestrial radiation
            DSINB = 3600.0 * (av.DAYL * av.SINLD + 24.0 * av.COSLD * sqrt_1_AOB2 / PI);
            // DSINBE: like DSINB but weighted by (1 + 0.4*sin(beta)) for diffuse sky radiation
            // This weighting accounts for the higher diffuse radiation near the horizon
            av.DSINBE = 3600.0 * (av.DAYL * (av.SINLD + 0.4 * (av.SINLD * av.SINLD + av.COSLD * av.COSLD * 0.5)) +
                                  12.0 * av.COSLD * (2.0 + 3.0 * 0.4 * av.SINLD) * sqrt_1_AOB2 / PI);
        }
        else
        {
            // Polar day or polar night
            if (AOB > 1.0)
                av.DAYL = 24.0;  // Midnight sun: sun never sets
            if (AOB < -1.0)
                av.DAYL = 0.0;   // Polar night: sun never rises
            DSINB = 3600.0 * (av.DAYL * av.SINLD);
            av.DSINBE = 3600.0 * (av.DAYL * (av.SINLD + 0.4 * (av.SINLD * av.SINLD + av.COSLD * av.COSLD * 0.5)));
        }

        /* ---- Photoperiodic daylength (with civil twilight, angle = -4 deg) ---- */
        real sin_angle = sin(ANGLE * RAD);
        real AOB_CORR = (-sin_angle + av.SINLD) / av.COSLD;

        if (std::abs(AOB_CORR) <= 1.0)
        {
            av.DAYLP = 12.0 * (1.0 + 2.0 * asin(AOB_CORR) / PI); // Photoperiodic daylength [h]
        }
        else if (AOB_CORR > 1.0)
        {
            av.DAYLP = 24.0; // Sun always above -4 deg
        }
        else
        {
            av.DAYLP = 0.0;  // Sun never reaches -4 deg
        }

        /* ---- Atmospheric transmission and diffuse radiation ---- */
        // Angot extraterrestrial radiation: theoretical maximum radiation at the
        // top of the atmosphere for the given day and latitude [J m^-2 d^-1]
        av.ANGOT = SC * DSINB;

        // Atmospheric transmission: ratio of measured to extraterrestrial radiation [-]
        // Range: 0 (heavy overcast) to ~0.75 (clear sky)
        if (av.DAYL > 0.0 && av.ANGOT > 0.0)
        {
            av.ATMTR = AVRAD / av.ANGOT;
        }
        else
        {
            av.ATMTR = 0.0;
        }

        // Diffuse radiation fraction (FRDIF) as a function of atmospheric transmission.
        // Empirical relationships from Spitters et al. (1986):
        //   - Clear sky (ATMTR > 0.75): ~23% diffuse
        //   - Partly cloudy (0.35 < ATMTR < 0.75): increasing diffuse fraction
        //   - Overcast (ATMTR < 0.07): nearly 100% diffuse
        real FRDIF;
        if (av.ATMTR > 0.75)
        {
            FRDIF = 0.23;
        }
        else if (av.ATMTR > 0.35)
        {
            FRDIF = 1.33 - 1.46 * av.ATMTR;
        }
        else if (av.ATMTR > 0.07)
        {
            FRDIF = 1.0 - 2.3 * (av.ATMTR - 0.07) * (av.ATMTR - 0.07);
        }
        else
        {
            FRDIF = 1.0;
        }

        // Diffuse PAR at the top of the canopy [J m^-2 s^-1]
        // FRDIF * ATMTR * 0.5 * SC: diffuse fraction of the transmitted extraterrestrial
        // radiation, with 0.5 factor for the PAR portion of the solar spectrum (400-700 nm)
        av.DIFPP = FRDIF * av.ATMTR * 0.5 * SC;

#if DEBUG_CROP_GROWTH_MODEL
        printf("[ASSIM-DEBUG-ASTRO] IDAY=%.1f, SC=%.2f, DEC=%.4f, SINLD=%.4f, COSLD=%.4f, DAYL=%.4f, DSINBE=%.4f, ANGOT=%.4f, ATMTR=%.4f, FRDIF=%.4f, DIFPP=%.4f\n",
               IDAY, SC, DEC, av.SINLD, av.COSLD, av.DAYL, av.DSINBE, av.ANGOT, av.ATMTR, FRDIF, av.DIFPP);
#endif

        return av;
    }

    /*====================================================================*
     *  INSTANTANEOUS CANOPY ASSIMILATION (ASSIM7)                        *
     *====================================================================*/

    /**
     * @brief Calculate instantaneous gross canopy CO2 assimilation rate.
     *
     * Implements the WOFOST ASSIM7 routine using a three-point Gaussian
     * integration over canopy depth. The canopy is divided into sunlit and
     * shaded leaf fractions, each with distinct light-response curves.
     * This is a key routine that couples the radiation regime to the
     * leaf-level photosynthesis model.
     *
     * @param[in] AMAX   Maximum leaf CO2 assimilation rate [kg CO2 ha^-1 h^-1].
     *                    Typical range: 40-100 for C3 crops, 70-100 for C4 crops.
     * @param[in] EFF    Initial light use efficiency of a single leaf
     *                    [kg CO2 ha^-1 h^-1 / (J m^-2 s^-1)].
     *                    Typical range: 0.4-0.5 for C3, 0.45-0.55 for C4.
     * @param[in] LAI    Leaf Area Index [m^2 leaf m^-2 ground]. Range: 0-12.
     * @param[in] KDIF   Extinction coefficient for diffuse visible radiation [-].
     *                    Typical range: 0.5-0.8.
     * @param[in] SINB   Sine of solar elevation [-]. Range: 0-1.
     * @param[in] PARDIR Direct photosynthetically active radiation (PAR) [J m^-2 s^-1].
     * @param[in] PARDIF Diffuse PAR [J m^-2 s^-1].
     * @param[in] debug  Enable debug output
     * @return    FGROS  Instantaneous gross canopy CO2 assimilation [kg CO2 ha^-1 h^-1].
     *
     * Algorithm:
     *   1. Compute canopy reflection coefficients for horizontal (REFH) and
     *      spherical (REFS) leaf angle distributions.
     *   2. Compute extinction coefficients for direct beam on horizontal (KDIRBL)
     *      and spherical (KDIRT) leaf surfaces.
     *   3. At each of 3 Gaussian integration depths through the canopy:
     *      a. Compute absorbed diffuse, direct total, and direct beam radiation.
     *      b. Compute shaded leaf assimilation (FGRSH) using a negative exponential
     *         light-response curve.
     *      c. Compute sunlit leaf assimilation (FGRSUN) including the direct beam
     *         contribution, using the full non-linear response.
     *      d. Weight sunlit (FSLLA) and shaded (1-FSLLA) contributions.
     *   4. Integrate across depths using Gaussian weights.
     *
     * Gaussian quadrature points and weights (3-point):
     *   XGAUSS = {0.1127, 0.5000, 0.8873}  (canonical positions on [0,1])
     *   WGAUSS = {0.2778, 0.4444, 0.2778}  (corresponding weights)
     */
    KOKKOS_INLINE_FUNCTION
    real assim7(real AMAX, real EFF, real LAI, real KDIF, real SINB, real PARDIR, real PARDIF, bool debug) const
    {
        // Three-point Gaussian quadrature nodes and weights for integration over canopy depth
        const real XGAUSS[3] = {0.1127017, 0.5000000, 0.8872983}; // Canonical positions on [0, 1]
        const real WGAUSS[3] = {0.2777778, 0.4444444, 0.2777778}; // Corresponding integration weights

        // SCV: scattering coefficient of leaves for visible radiation [-]
        // Typical value: 0.2 (20% of incident PAR is scattered by reflection/transmission)
        const real SCV = 0.2;

        /* ---- Canopy reflection and extinction coefficients ---- */
        // Reflection coefficient for a horizontal leaf angle distribution
        real sqrt_1_SCV = sqrt(1.0 - SCV);
        real REFH = (1.0 - sqrt_1_SCV) / (1.0 + sqrt_1_SCV);

        // Reflection coefficient for a spherical leaf angle distribution
        // (lower than REFH because more leaves face away from direct beam)
        real REFS = REFH * 2.0 / (1.0 + 1.6 * SINB);

        // Extinction coefficient for direct beam radiation on a horizontal surface
        // KDIRBL increases as solar elevation decreases (lower sun = more extinction per unit LAI)
        real KDIRBL = (0.5 / max(0.001, SINB)) * KDIF / (0.8 * sqrt_1_SCV);

        // Extinction coefficient for total direct radiation (including scattered component)
        real KDIRT = KDIRBL * sqrt_1_SCV;

        /* ---- Gaussian integration over canopy depth ---- */
        real FGROS = 0.0; // Accumulated instantaneous gross assimilation [kg CO2 ha^-1 h^-1]

        for (int i = 0; i < 3; ++i)
        {
            // Cumulative LAI at this Gaussian integration point
            real LAIC = LAI * XGAUSS[i]; // LAI from canopy top to this depth [m^2 leaf m^-2 ground]

            // Exponential attenuation factors for each radiation component
            real exp_KDIF = exp(-KDIF * LAIC);    // Fraction of diffuse PAR remaining above this depth
            real exp_KDIRT = exp(-KDIRT * LAIC);  // Fraction of total direct PAR remaining
            real exp_KDIRBL = exp(-KDIRBL * LAIC);// Fraction of direct beam (no scatter) remaining

            // Absorbed radiation fluxes at this canopy depth [J m^-2 s^-1]
            // VISDF: absorbed diffuse radiation per unit leaf area at this depth
            real VISDF = (1.0 - REFS) * PARDIF * KDIF * exp_KDIF;
            // VIST: absorbed total direct radiation (beam + scattered) per unit leaf area
            real VIST = (1.0 - REFS) * PARDIR * KDIRT * exp_KDIRT;
            // VISD: absorbed direct beam radiation only (no scatter) per unit leaf area
            real VISD = (1.0 - SCV) * PARDIR * KDIRBL * exp_KDIRBL;

            // Shaded leaf assimilation: uses only indirect radiation (diffuse + scattered direct)
            real VISSHD = VISDF + VIST - VISD; // Absorbed radiation by shaded leaves [J m^-2 s^-1]

            // Negative exponential light-response curve for shaded leaves:
            // FGRSH = AMAX * (1 - exp(-VISSHD * EFF / AMAX))
            real AMAX_eff = max(2.0, AMAX); // Ensure AMAX_eff >= 2 to avoid division issues
            real FGRSH = AMAX * (1.0 - exp(-VISSHD * EFF / AMAX_eff)); // Shaded leaf rate [kg CO2 ha^-1 h^-1]

            // Direct radiation flux perpendicular to the sun direction (for sunlit leaves)
            real VISPP = (1.0 - SCV) * PARDIR / max(0.001, SINB); // [J m^-2 s^-1]

            // Sunlit leaf assimilation using the full analytical solution of the
            // negative exponential light-response curve (Spitters, 1986)
            real FGRSUN;
            if (VISPP <= 0.0)
            {
                // No direct beam: sunlit leaves behave same as shaded
                FGRSUN = FGRSH;
            }
            else
            {
                // Full sunlit leaf rate: integrates the light response over the
                // leaf angle distribution within the sunlit fraction
                FGRSUN = AMAX * (1.0 - (AMAX - FGRSH) * (1.0 - exp(-VISPP * EFF / AMAX_eff)) / (EFF * VISPP));
            }

            // Fraction of sunlit leaves at this canopy depth
            real FSLLA = exp_KDIRBL; // Decreases exponentially with canopy depth

            // Combined assimilation at this depth: weighted sum of sunlit and shaded
            real FGL = FSLLA * FGRSUN + (1.0 - FSLLA) * FGRSH; // [kg CO2 ha^-1 h^-1]

            // Accumulate weighted contribution across Gaussian integration points
            FGROS += FGL * WGAUSS[i];

#if DEBUG_CROP_GROWTH_MODEL
            printf("[ASSIM-DEBUG-ASSIM7] i=%d, LAIC=%.6f, VISDF=%.6f, VIST=%.6f, VISD=%.6f, VISSHD=%.6f, FGRSH=%.6f, FGRSUN=%.6f, FSLLA=%.6f, FGL=%.6f\n",
                   i, LAIC, VISDF, VIST, VISD, VISSHD, FGRSH, FGRSUN, FSLLA, FGL);
#endif
        }

        // Scale by total LAI to convert from per-unit-LAI to per-hectare rate
        FGROS = FGROS * LAI; // [kg CO2 ha^-1 h^-1]
        return FGROS;
    }

    /*====================================================================*
     *  DAILY TOTAL ASSIMILATION (TOTASS7)                                *
     *====================================================================*/

    /**
     * @brief Calculate daily total gross canopy CO2 assimilation.
     *
     * Integrates the instantaneous assimilation rate (assim7) over the daylight
     * period using a three-point Gaussian quadrature scheme. At each time point,
     * the solar elevation is computed to determine the instantaneous direct and
     * diffuse PAR, which are then passed to assim7.
     *
     * @param[in] DAYL   Astronomical daylength [h]. Range: 0-24.
     * @param[in] AMAX   Maximum leaf CO2 assimilation rate [kg CO2 ha^-1 h^-1].
     * @param[in] EFF    Light use efficiency [kg CO2 ha^-1 h^-1 / (J m^-2 s^-1)].
     * @param[in] LAI    Leaf Area Index [m^2 m^-2].
     * @param[in] KDIF   Diffuse light extinction coefficient [-].
     * @param[in] AVRAD  Daily total global solar radiation [J m^-2 d^-1].
     * @param[in] DIFPP  Diffuse PAR at canopy top [J m^-2 s^-1].
     * @param[in] DSINBE Daily integral of sin(beta)*(1+0.4*sin(beta)) [s].
     * @param[in] SINLD  sin(latitude)*sin(declination) [-].
     * @param[in] COSLD  cos(latitude)*cos(declination) [-].
     * @param[in] debug  Enable debug output.
     * @return    DTGA   Daily total gross canopy CO2 assimilation [kg CO2 ha^-1 d^-1].
     *
     * Algorithm:
     *   1. For each of 3 Gaussian integration points across the daylength:
     *      a. Compute the hour of day and solar elevation (SINB).
     *      b. Partition PAR into diffuse (PARDIF) and direct (PARDIR) components.
     *      c. Call assim7 to compute instantaneous assimilation at this time.
     *   2. Integrate using Gaussian weights, then multiply by daylength.
     */
    KOKKOS_INLINE_FUNCTION
    real totass7(real DAYL, real AMAX, real EFF, real LAI, real KDIF, real AVRAD,
                 real DIFPP, real DSINBE, real SINLD, real COSLD, bool debug) const
    {
        // Three-point Gaussian quadrature nodes and weights for integration over the day
        const real XGAUSS[3] = {0.1127017, 0.5000000, 0.8872983};
        const real WGAUSS[3] = {0.2777778, 0.4444444, 0.2777778};

        real DTGA = 0.0; // Accumulated daily total gross assimilation [kg CO2 ha^-1 d^-1]

        // Only compute if there is daylight, leaf area, and photosynthetic capacity
        if (AMAX > 0.0 && LAI > 0.0 && DAYL > 0.0)
        {
            for (int i = 0; i < 3; ++i)
            {
                // Compute the hour of day for this Gaussian point
                // HOUR ranges from 12 - DAYL/2 (sunrise) to 12 + DAYL/2 (sunset)
                real HOUR = 12.0 + 0.5 * DAYL * XGAUSS[i]; // Hour of day [h]

                // Compute sin(solar elevation) at this hour
                // The cosine term represents the diurnal variation of solar elevation
                real cos_term = cos(2.0 * PI * (HOUR + 12.0) / 24.0);
                real SINB = max(0.0, SINLD + COSLD * cos_term); // sin(solar elevation) [-], clamped >= 0

                // Compute total photosynthetically active radiation (PAR) at this instant
                // PAR = 0.5 * AVRAD * SINB * (1 + 0.4*SINB) / DSINBE
                // The factor 0.5 accounts for the PAR fraction of total solar radiation.
                // The weighting (1 + 0.4*SINB) approximates the non-linear relationship
                // between solar elevation and atmospheric transmission.
                real safe_DSINBE = max(1.0, DSINBE); // Prevent division by zero
                real PAR = 0.5 * AVRAD * SINB * (1.0 + 0.4 * SINB) / safe_DSINBE; // Total PAR [J m^-2 s^-1]

                // Partition PAR into diffuse and direct components
                // Diffuse PAR cannot exceed the product of sin(elevation) and DIFPP
                real PARDIF = min(PAR, SINB * DIFPP); // Diffuse PAR [J m^-2 s^-1]
                real PARDIR = PAR - PARDIF;            // Direct PAR [J m^-2 s^-1]

#if DEBUG_CROP_GROWTH_MODEL
                printf("[ASSIM-DEBUG-TOTASS7] i=%d, HOUR=%.4f, SINB=%.6f, PAR=%.6f, PARDIF=%.6f, PARDIR=%.6f\n",
                       i, HOUR, SINB, PAR, PARDIF, PARDIR);
#endif

                // Compute instantaneous gross assimilation at this hour
                real FGROS = assim7(AMAX, EFF, LAI, KDIF, SINB, PARDIR, PARDIF, debug);

                // Accumulate weighted contribution
                DTGA += FGROS * WGAUSS[i];

#if DEBUG_CROP_GROWTH_MODEL
                printf("[ASSIM-DEBUG-TOTASS7] i=%d, FGROS=%.6f, Accumulated DTGA=%.6f\n", i, FGROS, DTGA);
#endif
            }
        }

        // Multiply by daylength to convert from hourly rate to daily total
        // FGROS is in [kg CO2 ha^-1 h^-1], so DTGA = sum(FGROS * WGAUSS) * DAYL [kg CO2 ha^-1 d^-1]
        DTGA *= DAYL;
        return DTGA;
    }

    /*====================================================================*
     *  MAIN RATE CALCULATION ENTRY POINT                                 *
     *====================================================================*/

    /**
     * @brief Calculate daily potential crop assimilation for all cells.
     *
     * This is the main entry point for the assimilation module. It orchestrates
     * the full computation pipeline for each computational cell:
     *
     *   1. Retrieve meteorological variables and crop state (DVS, LAI).
     *   2. Compute astronomical variables (daylength, radiation partitioning).
     *   3. Look up photosynthetic parameters from crop parameter tables:
     *      - AMAXTB: max assimilation rate vs. DVS [kg CO2 ha^-1 h^-1]
     *      - EFFTB:  light use efficiency vs. temperature [kg CO2 ha^-1 h^-1 / (J m^-2 s^-1)]
     *      - KDIFTB: diffuse extinction coefficient vs. DVS [-]
     *      - TMPFTB: temperature correction factor for AMAX vs. temperature [-]
     *      - TMNFTB: low-temperature reduction factor vs. Tmin [-]
     *   4. Compute daily total gross CO2 assimilation (DTGA).
     *   5. Apply low-temperature stress correction (TMNF).
     *   6. Convert CO2 to CH2O (carbohydrate) equivalents (30/44 ratio).
     *
     * @param[out] as                   Assimilation output state (PGASS written)
     * @param[in]  p                    Static crop parameters and AFGEN tables
     * @param[in]  s                    Dynamic crop state (DVS, LAI)
     * @param[in]  m                    Meteorological forcing data
     * @param[in]  drv_idx              Index into meteorological time-series arrays
     * @param[in]  current_day_of_year  Day of year (1-365) for astronomical calculations
     *
     * Output:
     *   - PGASS: Potential gross assimilation rate [kg CH2O ha^-1 d^-1].
     *     Typical range: 0 (no canopy) to ~400 (full canopy, clear sky, warm).
     *
     * Temperature conventions:
     *   - TEMP  = (TMIN + TMAX) / 2 : 24-hour mean temperature [C]
     *   - DTEMP = TMIN + 0.75*(TMAX-TMIN) : daytime mean temperature [C],
     *     weighted toward TMAX because photosynthesis occurs during daylight
     *     when temperatures are closer to the daily maximum.
     */
    void calc_rates(AssimilateState &as, const CropState &p, const CropDynamicState &s,
                    const MeteoState &m, int drv_idx, int current_day_of_year)
    {
        // Retrieve crop parameter lookup tables
        auto AMAXTB = p.tables.at("AMAXTB"); // DVS -> max CO2 assimilation rate [kg CO2 ha^-1 h^-1]
        auto EFFTB = p.tables.at("EFFTB");   // Temperature -> light use efficiency [kg CO2 ha^-1 h^-1 / (J m^-2 s^-1)]
        auto KDIFTB = p.tables.at("KDIFTB"); // DVS -> diffuse light extinction coefficient [-]
        auto TMPFTB = p.tables.at("TMPFTB"); // Temperature -> AMAX correction factor [-], range [0, 1]
        auto TMNFTB = p.tables.at("TMNFTB"); // Tmin -> low-temperature reduction factor [-], range [0, 1]

        Kokkos::parallel_for("Assim_Rates", as.nCells, KOKKOS_LAMBDA(const int i) {
            bool debug = (i == 0); // Debug output for the first cell only

            // Current crop state
            real DVS = s.DVS(i);  // Development Stage [-], range: -0.1 to ~2.0
            real LAI = s.LAI(i);  // Leaf Area Index [m^2 m^-2], range: 0 to ~12

            // Meteorological forcing for this day
            real TMIN = m.tmin(drv_idx); // Daily minimum temperature [C]
            real TMAX = m.tmax(drv_idx); // Daily maximum temperature [C]

            // Compute two temperature metrics:
            // TEMP: 24-hour average temperature, used for general reference
            real TEMP = (TMIN + TMAX) / 2.0;          // 24-hour mean [C]
            // DTEMP: daytime-weighted average temperature, used for photosynthesis
            // calculations because photosynthesis occurs during daylight when
            // temperatures are higher (TMAX-weighted by factor 0.75)
            real DTEMP = TMIN + (TMAX - TMIN) * 0.75; // Daytime mean [C]

            real IRRAD = m.irrad(drv_idx); // Daily global solar radiation [J m^-2 d^-1]
            real LAT = m.lat;              // Latitude [degrees]
#if DEBUG_CROP_GROWTH_MODEL
            printf("[ASSIM-DEBUG] cell=%d, DVS=%.4f, LAI=%.6f, TEMP(24h)=%.2f, DTEMP(Day)=%.2f, TMIN=%.2f, IRRAD=%.1f, LAT=%.2f, lon=%.2f\n",
                   i, DVS, LAI, TEMP, DTEMP, TMIN, IRRAD, LAT, m.lon);
#endif
            // Note: TMINRA should ideally be a 7-day rolling average of daily minimum
            // temperature, but is here simplified to the current day's TMIN.
            real TMINRA = TMIN; // 7-day running average of Tmin [C] (simplified)

            /* ---- Astronomical calculations ---- */
            AstroVars av = astro(current_day_of_year, LAT, IRRAD, debug);

            /* ---- Photosynthetic parameter lookup ---- */
            // AMAX_BASE: base maximum assimilation rate, dependent on crop development stage
            real AMAX_BASE = Afgen::lookup(AMAXTB, DVS); // [kg CO2 ha^-1 h^-1]

            // TMPF: temperature correction factor for AMAX.
            // Uses daytime temperature (DTEMP) because photosynthesis is a daytime process.
            // Range: 0 (extreme temperatures) to 1 (optimal temperature)
            real TMPF = Afgen::lookup(TMPFTB, DTEMP); // [-]
            real AMAX = AMAX_BASE * TMPF; // Corrected maximum assimilation rate [kg CO2 ha^-1 h^-1]

            // KDIF: diffuse light extinction coefficient, depends on canopy structure
            // which changes with development stage
            real KDIF = Afgen::lookup(KDIFTB, DVS); // [-], range: 0.4-1.0

            // EFF: initial light use efficiency of leaves at low light intensity.
            // Temperature-dependent; uses daytime temperature (DTEMP).
            // Typical range: 0.4-0.55 [kg CO2 ha^-1 h^-1 / (J m^-2 s^-1)]
            real EFF = Afgen::lookup(EFFTB, DTEMP);

#if DEBUG_CROP_GROWTH_MODEL
            printf("[ASSIM-DEBUG] cell=%d, DVS=%.4f, LAI=%.6f, TEMP(24h)=%.2f, DTEMP(Day)=%.2f, TMIN=%.2f, IRRAD=%.1f, LAT=%.2f\n",
                   i, DVS, LAI, TEMP, DTEMP, TMIN, IRRAD, LAT);
            printf("[ASSIM-DEBUG] AMAX_BASE=%.4f, TMPF=%.4f, AMAX=%.4f, KDIF=%.4f, EFF=%.4f\n", AMAX_BASE, TMPF, AMAX, KDIF, EFF);
            printf("[ASSIM-DEBUG] DAYL=%.4f, SINLD=%.4f, COSLD=%.4f, DIFPP=%.4f, DSINBE=%.4f\n",
                   av.DAYL, av.SINLD, av.COSLD, av.DIFPP, av.DSINBE);
#endif

            /* ---- Daily total gross CO2 assimilation ---- */
            real DTGA = totass7(av.DAYL, AMAX, EFF, LAI, KDIF, IRRAD, av.DIFPP, av.DSINBE, av.SINLD, av.COSLD, debug);
            // DTGA: daily total gross canopy CO2 assimilation [kg CO2 ha^-1 d^-1]

            /* ---- Low-temperature stress correction ---- */
            // TMNF: reduction factor accounting for sub-optimal nighttime temperatures.
            // Frost or very low temperatures can reduce photosynthesis even during the day.
            // Range: 0 (complete inhibition) to 1 (no stress)
            real TMNF = Afgen::lookup(TMNFTB, TMINRA);
            DTGA *= TMNF; // Apply low-temperature correction

            /* ---- Convert CO2 to carbohydrate (CH2O) equivalents ---- */
            // Molecular weight conversion: CO2 (44 g/mol) -> CH2O (30 g/mol)
            // This accounts for the fact that the carbohydrate produced by
            // photosynthesis (CH2O) is lighter than the CO2 consumed.
            as.PGASS(i) = DTGA * 30.0 / 44.0; // Potential gross assimilation [kg CH2O ha^-1 d^-1]

#if DEBUG_CROP_GROWTH_MODEL
            printf("[ASSIM-DEBUG] DTGA_FINAL=%.6f, TMNF=%.4f, PGASS=%.6f\n", DTGA, TMNF, as.PGASS(i));
#endif
        });
    }
};

#endif
