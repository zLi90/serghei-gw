/* -*- mode: c++; c-default-style: "linux" -*- */
/**
 * @file Evapotranspiration.h
 * @brief Evapotranspiration calculations for the WOFOST crop growth model.
 *
 * This file implements the computation of potential and actual evaporation
 * and transpiration rates, including water-stress and oxygen-stress reduction
 * factors that regulate crop water use.
 *
 * The module computes the following quantities for each computational cell:
 *
 *   - EVS   : Maximum evaporation rate from a wet soil surface below the
 *             crop canopy, reduced by canopy shading [cm d^-1].
 *   - TRAMX : Maximum (potential) transpiration rate from the crop canopy [cm d^-1].
 *   - RFTRA : Combined reduction factor for transpiration due to water and/or
 *             oxygen stress [-]. Range: 0 (complete stress) to 1 (no stress).
 *   - TRA   : Actual transpiration rate = TRAMX * RFTRA [cm d^-1].
 *
 * The E0, ES0, and ET0 reference values are read from MeteoState (pre-computed
 * during MeteoInit using Penman and Penman-Monteith equations).
 *
 * Water stress (RFWS):
 *   - Occurs when root-zone soil moisture (SM) drops below a critical threshold (SMCR).
 *   - SMCR is derived from soil hydraulic properties (field capacity SMFCF, wilting
 *     point SMW) and the crop's rooting depth factor (DEPNR) using the SWEAF function.
 *   - RFWS decreases linearly from 1.0 at SMCR to 0.0 at SMW (wilting point).
 *
 * Oxygen stress (RFOS):
 *   - Occurs when root-zone soil moisture approaches saturation, limiting
 *     oxygen availability for root respiration.
 *   - Controlled by parameters SM0 (saturated moisture content) and CRAIRC
 *     (critical air content for root aeration).
 *   - Only applied when the oxygen stress flag is enabled (IOX = 1).
 *
 * References:
 *   - WOFOST 7.1 User Guide, Alterra, Wageningen UR
 *   - Supit, I., et al. (1994). System Description of the WOFOST 6.0 Crop
 *     Growth Simulation Model. EUR 14656 EN.
 *
 * @see EvapotranspirationState.h for output state arrays
 * @see CropState.h              for static crop parameters (CFET, SMFCF, SMW, DEPNR, IOX)
 * @see CropDynamicState.h       for dynamic crop state (DVS, LAI)
 * @see MeteoState.h             for meteorological forcing (E0, ES0, ET0)
 */

#ifndef _EVAPOTRANSPIRATION_H_
#define _EVAPOTRANSPIRATION_H_

#include "../define.h"
#include "CropState.h"               // Static crop parameters
#include "CropDynamicState.h"        // DVS, LAI, SM (from SERGHEI hydrological model)
#include "MeteoState.h"              // TMIN, TMAX, IRRAD, VAP, WIND (for ET0 calculation)
#include "EvapotranspirationState.h" // Output TRA, TRAMX, EVS, RFTRA
#include "Afgen.h"                   // For KDIFTB lookup
#include <cmath>
#include <algorithm> // For fmin, fmax

/*======================================================================*
 *  SOIL WATER EVAPORATION AVAILABILITY FACTOR (SWEAF)                   *
 *======================================================================*/

/**
 * @brief Compute the Soil Water Evaporative Availability Factor (SWEAF).
 *
 * This function calculates the fraction of available soil water that can be
 * freely evaporated or transpired before water stress occurs. The factor
 * depends on the potential evapotranspiration rate (ET0) and the crop's
 * rooting depth indicator (DEPNR). Deeper-rooted crops can access more soil
 * water before experiencing stress.
 *
 * @param[in] ET0_val  Potential evapotranspiration from reference crop [cm d^-1].
 *                     Typical range: 0.05-1.0.
 * @param[in] DEPNR    Crop rooting depth factor [-].
 *                     Range: 1 (shallow-rooted) to 5 (deep-rooted).
 *                     Typical values: 1.5 (lettuce), 3.0 (wheat), 4.5 (sugar beet).
 * @return    sweaf    Soil water availability factor [-]. Clamped to [0.10, 0.95].
 *
 * Physical basis:
 *   - Empirical coefficients A = 0.76 and B = 1.5 are derived from
 *     WOFOST calibration against field measurements.
 *   - For shallow-rooted crops (DEPNR < 3), an additional correction term
 *     accounts for the reduced water buffer in the shallow root zone.
 *   - The result is bounded to prevent physically unrealistic values:
 *     minimum 0.10 (always some water available) and maximum 0.95
 *     (never 100% available, some water is always tightly bound to soil).
 */
KOKKOS_INLINE_FUNCTION
real SWEAF(real ET0_val, real DEPNR)
{
    const real A = 0.76;   // Empirical intercept coefficient [-]
    const real B = 1.5;    // Empirical slope coefficient [-]

    // Base SWEAF calculation: inverse function of ET0, reduced by deeper rooting
    real sweaf = 1.0 / (A + B * ET0_val) - (5.0 - DEPNR) * 0.10;

    // Additional correction for shallow-rooted crops (DEPNR < 3)
    // These crops deplete their smaller root zone faster at high ET0
    if (DEPNR < 3.0)
    {
        sweaf += (ET0_val - 0.6) / (DEPNR * (DEPNR + 3.0));
    }

    // Clamp to physically reasonable range
    if (sweaf < 0.10) sweaf = 0.10;  // Minimum: at least 10% available
    if (sweaf > 0.95) sweaf = 0.95;  // Maximum: no more than 95% available

    return sweaf;
}

/*======================================================================*
 *  EVAPOTRANSPIRATION CLASS                                             *
 *======================================================================*/

class Evapotranspiration
{

public:

    /*==================================================================*
     *  RATE CALCULATION                                                 *
     *==================================================================*/

    /**
     * @brief Calculate potential and actual evaporation/transpiration rates.
     *
     * Computes the daily evapotranspiration components for each computational
     * cell, including maximum soil evaporation, maximum canopy transpiration,
     * water-stress and oxygen-stress reduction factors, and actual transpiration.
     *
     * @param[out] ets            Evapotranspiration state (TRA, TRAMX, EVS, RFTRA, E0, ES0, ET0)
     * @param[in]  p              Static crop parameters (CFET, SMFCF, SMW, DEPNR, IOX, SM0, CRAIRC)
     * @param[in]  s              Dynamic crop state (DVS for table lookups, LAI for shading)
     * @param[in]  m              Meteorological forcing (pre-computed E0, ES0, ET0)
     * @param[in]  drv_idx        Index into meteorological time-series arrays
     * @param[in]  root_zone_SM   Average volumetric soil moisture content in the root zone
     *                            [cm^3 cm^-3], provided by the SERGHEI hydrological model.
     *                            Range: SMW (wilting point, ~0.05-0.15) to SM0 (saturation, ~0.35-0.55).
     *
     * Calculation sequence:
     *   1. Read E0, ES0, ET0 from MeteoState (pre-computed by MeteoInit)
     *   2. Apply crop-specific correction: ET0_CROP = CFET * ET0
     *   3. Compute canopy light extinction: KGLOB = 0.75 * KDIF(DVS)
     *   4. Compute maximum soil evaporation: EVS = ES0 * exp(-KGLOB * LAI)
     *   5. Compute maximum transpiration: TRAMX = ET0_CROP * (1 - exp(-KGLOB * LAI))
     *   6. Compute critical soil moisture: SMCR using SWEAF function
     *   7. Compute water stress factor: RFWS = (SM - SMW) / (SMCR - SMW), clamped [0, 1]
     *   8. Compute oxygen stress factor: RFOS (if IOX = 1)
     *   9. Compute total reduction: RFTRA = RFOS * RFWS
     *  10. Compute actual transpiration: TRA = TRAMX * RFTRA
     */
    void calc_rates(EvapotranspirationState &ets, const CropState &p, const CropDynamicState &s,
                    const MeteoState &m, int drv_idx, const realArr &root_zone_SM)
    {
        // Retrieve the diffuse light extinction coefficient table
        auto KDIFTB = p.tables.at("KDIFTB");

        Kokkos::parallel_for("EvapoTrans_Rates", ets.nCells, KOKKOS_LAMBDA(const int i) {

            /* ---- Current crop and soil state ---- */
            real DVS = s.DVS(i);     // Development Stage [-], range: -0.1 to ~2.0
            real LAI = s.LAI(i);     // Leaf Area Index [m^2 m^-2], range: 0-12
            real SM = root_zone_SM(i); // Volumetric soil moisture in root zone [cm^3 cm^-3]

            /* ---- Step 1: Read pre-computed reference ET values ---- */
            // These were calculated in MeteoInit using Penman/Penman-Monteith equations.
            // Units: cm/day (WOFOST convention, not mm/day).
            ets.E0(i) = m.e0(drv_idx);    // Potential evaporation from open water [cm d^-1]
            ets.ES0(i) = m.es0(drv_idx);  // Potential evaporation from bare soil [cm d^-1]
            ets.ET0(i) = m.et0(drv_idx);  // Potential evapotranspiration from reference crop [cm d^-1]

            /* ---- Step 2: Crop-specific ET correction ---- */
            // CFET is a crop-specific correction factor applied to the reference ET0.
            // CFET = 1.0 for the standard reference crop; may differ for other crops.
            // Range: CFET > 0. Typical values: 0.8-1.3.
            real ET0_CROP = fmax(0.0, p.p.CFET * ets.ET0(i)); // Crop-corrected potential ET [cm d^-1]

            /* ---- Step 3: Canopy radiation extinction ---- */
            // KGLOB: global extinction coefficient for total (direct + diffuse) radiation [-].
            // Derived from the diffuse extinction coefficient KDIF, scaled by 0.75.
            // KGLOB determines how effectively the canopy intercepts radiation and,
            // consequently, how much radiation reaches the soil surface.
            real KGLOB = 0.75 * Afgen::lookup(KDIFTB, DVS);
            real EKL = exp(-KGLOB * LAI); // Fraction of radiation reaching the soil surface [-]

            /* ---- Step 4: Maximum soil evaporation (EVS) ---- */
            // Soil evaporation decreases exponentially with LAI because the canopy
            // shades the soil surface. With full canopy closure (high LAI), very
            // little radiation reaches the soil, and evaporation is minimal.
            ets.EVS(i) = ets.ES0(i) * EKL; // Maximum soil evaporation [cm d^-1]

            /* ---- Step 5: Maximum transpiration (TRAMX) ---- */
            // The canopy transpires the fraction of ET0 not allocated to soil evaporation.
            // (1 - EKL) represents the fraction of energy intercepted by the canopy.
            ets.TRAMX(i) = ET0_CROP * (1.0 - EKL); // Maximum transpiration [cm d^-1]

#if DEBUG_CROP_GROWTH_MODEL
            printf("Evapotranspiration: i=%d, DVS=%.2f, LAI=%.2f, SM=%.2f, E0=%.6f, ES0=%.6f, ET0=%.6f, EVS=%.6f, TRAMX=%.6f\n",
                   i, DVS, LAI, SM, ets.E0(i), ets.ES0(i), ets.ET0(i), ets.EVS(i), ets.TRAMX(i));
#endif

            /* ---- Step 6: Critical soil moisture for water stress ---- */
            // SMFCF: soil moisture at field capacity [cm^3 cm^-3].
            //   The water content after drainage has slowed to negligible rates.
            //   Typical range: 0.20-0.45 depending on soil texture.
            // SMW: soil moisture at wilting point [cm^3 cm^-3].
            //   The water content below which plant roots cannot extract water.
            //   Typical range: 0.03-0.20 depending on soil texture.
            // DEPNR: crop rooting depth factor [-].
            //   Indicates the effective rooting depth for water extraction.
            //   Range: 1 (shallow) to 5 (deep). Typical: 3.0 for cereals.
            real SMFCF = p.p.SMFCF;  // Field capacity [cm^3 cm^-3]
            real SMW = p.p.SMW;      // Wilting point [cm^3 cm^-3]
            real DEPNR = p.p.DEPNR;  // Rooting depth factor [-]

#if DEBUG_CROP_GROWTH_MODEL
            printf("Evapotranspiration Debug: i=%d, SMFCF=%.4f, SMW=%.4f, DEPNR=%.4f\n", i, SMFCF, SMW, DEPNR);
#endif

            // SWDEP: soil water depletion factor from the SWEAF function [-].
            // Determines what fraction of the available water (between SMW and SMFCF)
            // can be extracted before water stress begins.
            real SWDEP = SWEAF(ET0_CROP, DEPNR);

            // SMCR: critical soil moisture content [cm^3 cm^-3].
            // Below SMCR, the crop experiences water stress.
            // SMCR = (1 - SWDEP) * (SMFCF - SMW) + SMW
            real SMCR = (1.0 - SWDEP) * (SMFCF - SMW) + SMW;

            /* ---- Step 7: Water stress reduction factor (RFWS) ---- */
            // Linear reduction from 1.0 (no stress) at SM >= SMCR to 0.0 (full stress)
            // at SM <= SMW (wilting point).
            //   SM >= SMCR : RFWS = 1.0 (no water stress, sufficient soil moisture)
            //   SM <= SMW  : RFWS = 0.0 (severe water stress, transpiration ceases)
            //   SMW < SM < SMCR : RFWS = (SM - SMW) / (SMCR - SMW) (linear interpolation)
            real RFWS = (SM - SMW) / (SMCR - SMW);
            if (RFWS < 0.0) RFWS = 0.0;  // Cannot be negative
            if (RFWS > 1.0) RFWS = 1.0;  // Cannot exceed 1.0

            /* ---- Step 8: Oxygen stress reduction factor (RFOS) ---- */
            // Oxygen stress occurs in waterlogged soils when the air-filled pore
            // space drops below a critical threshold (CRAIRC), depriving roots of
            // oxygen for respiration. Only computed when the oxygen stress flag is
            // enabled (IOX = 1).
            real RFOS = 1.0; // Default: no oxygen stress
            if (p.p.IOX == 1) {
                // SM0: saturated soil moisture content [cm^3 cm^-3].
                //   The maximum water-holding capacity of the soil. All pores are filled.
                //   Typical range: 0.35-0.55 depending on soil texture.
                // CRAIRC: critical air content for root aeration [cm^3 cm^-3].
                //   The minimum air-filled pore space needed for root respiration.
                //   Typical range: 0.03-0.10.
                real SM0 = p.p.SM0;
                real CRAIRC = p.p.CRAIRC;

                // RFOSMX: ratio of actual air content to critical air content [-].
                //   Air content = SM0 - SM (porosity minus actual water content).
                //   When SM approaches SM0 (saturation), air content approaches 0
                //   and RFOS approaches 0 (complete oxygen stress, transpiration stops).
                //   When SM <= SM0 - CRAIRC, air content exceeds the critical threshold
                //   and RFOS = 1.0 (no oxygen stress).
                real RFOSMX = (SM0 - SM) / CRAIRC;
                RFOS = std::max(0.0, std::min(1.0, RFOSMX)); // Clamp to [0, 1]
            }

            /* ---- Step 9: Combined transpiration reduction factor ---- */
            // The total reduction factor is the product of water stress and oxygen stress.
            // This means either type of stress alone can reduce transpiration.
            ets.RFTRA(i) = RFOS * RFWS; // Combined reduction factor [-], range [0, 1]

#if DEBUG_CROP_GROWTH_MODEL
            printf("Evapotranspiration Debug: i=%d, SM=%.4f, SMCR=%.4f, RFWS=%.4f, RFOS=%.4f, RFTRA=%.4f\n", i, SM, SMCR, RFWS, RFOS, ets.RFTRA(i));
#endif

            /* ---- Step 10: Actual transpiration rate ---- */
            // Actual transpiration is the maximum rate reduced by the combined
            // stress factor. This is the water actually used by the crop.
            ets.TRA(i) = ets.TRAMX(i) * ets.RFTRA(i); // Actual transpiration [cm d^-1]
        });
    }
};

#endif
