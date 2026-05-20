/* -*- mode: c++; c-default-style: "linux" -*- */

/*******************************************************************************
 * @file    Respiration.h
 * @brief   Maintenance respiration calculations for the WOFOST crop growth model.
 *
 * @details This module computes the potential maintenance respiration rate (PMRES)
 *          for each computational cell, accounting for:
 *            1. Organ-specific maintenance respiration coefficients
 *            2. Senescence reduction via a development-stage-dependent lookup table
 *            3. Temperature correction using the Q10 exponential model
 *
 *          The calculation follows the WOFOST (World Food Studies) approach:
 *
 *            RMRES = (RMR*WRT + RML*WLV + RMS*WST + RMO*WSO) * RFSETB(DVS)
 *            TEFF  = Q10 ^ ((TEMP - 25.0) / 10.0)
 *            PMRES = RMRES * TEFF
 *
 *          where:
 *            - RMRES is the temperature-uncorrected maintenance respiration [kg CH2O/ha/d]
 *            - TEFF  is the dimensionless temperature correction factor
 *            - PMRES is the final potential maintenance respiration rate [kg CH2O/ha/d]
 *            - Reference temperature for Q10 is 25 degrees Celsius
 *
 * @note    All biomass quantities use units of [kg dry matter / ha].
 *          Growth respiration (conversion of assimilates to structural biomass)
 *          is handled separately and is NOT part of this module.
 *
 * @see     RespirationState.h  for the output state container (PMRES array)
 * @see     CropState.h         for static crop parameters (RMR, RML, RMS, RMO, Q10)
 * @see     Afgen.h             for the AFGEN table interpolation utility
 ******************************************************************************/

#ifndef _RESPIRATION_H_
#define _RESPIRATION_H_

#include "../define.h"
#include "CropState.h"                  /* Static crop parameters (RMR, RML, RMS, RMO, Q10, RFSETB table) */
#include "CropDynamicState.h"           /* Dynamic state variables (DVS - development stage, 0=emergence to 2=maturity) */
#include "LeafDynamicsState.h"          /* Leaf biomass state (WLV) */
#include "StemDynamicsState.h"          /* Stem biomass state (WST) */
#include "RootDynamicsState.h"          /* Root biomass state (WRT) */
#include "StorageOrganDynamicsState.h"  /* Storage organ biomass state (WSO) */
#include "RespirationState.h"           /* Output state: PMRES array */
#include "MeteoState.h"                 /* Meteorological driving variables (tmin, tmax for temperature) */
#include "Afgen.h"                      /* AFGEN tabular interpolation utility */
#include <cmath>

/**
 * @class Respiration
 * @brief Computes potential maintenance respiration rates for all crop organs.
 *
 * This class provides a single method, calc_rates(), which is called each
 * simulation day. It operates on Kokkos device arrays via a parallel_for
 * kernel so that respiration is computed simultaneously across all grid cells.
 */
class Respiration
{

public:
    /**
     * @brief Calculate maintenance respiration rates for every computational cell.
     *
     * This method executes a Kokkos parallel kernel that, for each cell:
     *   1. Retrieves live biomass of each organ (roots, leaves, stems, storage organs)
     *   2. Computes the base maintenance respiration as a weighted sum using
     *      organ-specific respiration coefficients (RMR, RML, RMS, RMO)
     *   3. Applies a senescence reduction factor from the RFSETB lookup table,
     *      which reduces respiration as the crop matures past development stage ~1.0
     *   4. Applies a temperature correction using the Q10 model, where Q10 is
     *      the factor by which respiration increases for every 10-degree-Celsius
     *      rise above the reference temperature of 25 degrees Celsius
     *   5. Stores the result in rs.PMRES
     *
     * @param[in,out] rs      RespirationState containing the output PMRES array [kg CH2O/ha/d]
     * @param[in]     p       CropState with static parameters and AFGEN tables
     * @param[in]     s        CropDynamicState with per-cell development stage (DVS) [0..2]
     * @param[in]     m        MeteoState with daily meteorological time series
     * @param[in]     drv_idx  Index into the meteo arrays for the current simulation day
     * @param[in]     lds      LeafDynamicsState with leaf live biomass WLV [kg/ha]
     * @param[in]     sds      StemDynamicsState with stem live biomass WST [kg/ha]
     * @param[in]     rds      RootDynamicsState with root live biomass WRT [kg/ha]
     * @param[in]     sods     StorageOrganDynamicsState with storage organ live biomass WSO [kg/ha]
     */
    void calc_rates(RespirationState &rs, const CropState &p, const CropDynamicState &s,
                    const MeteoState &m, int drv_idx,
                    const LeafDynamicsState &lds,
                    const StemDynamicsState &sds,
                    const RootDynamicsState &rds,
                    const StorageOrganDynamicsState &sods)
    {
        /* ------------------------------------------------------------------- *
         * Retrieve the AFGEN table for the senescence reduction factor (RFSETB).
         * This table maps development stage (DVS) to a dimensionless multiplier
         * in the range [0, 1]. Typical values:
         *   DVS < 1.0  : RFSETB = 1.0  (full respiration, no senescence)
         *   DVS ~ 1.3  : RFSETB ~ 0.5  (partial senescence)
         *   DVS >= 2.0 : RFSETB = 0.0  (complete senescence, no respiration)
         * ------------------------------------------------------------------- */
        auto RFSETB = p.tables.at("RFSETB");

        Kokkos::parallel_for("Respiration_Rates", rs.nCells, KOKKOS_LAMBDA(const int i) {

            /* ================================================================= *
             * Section 1: Retrieve organ biomass state variables
             *
             * Each variable represents the live (non-senesced) dry matter weight
             * of the respective organ. Units: [kg dry matter / ha]
             *
             * Typical ranges for a well-developed cereal crop (e.g., wheat):
             *   WRT (roots):         50 - 500 kg/ha
             *   WLV (leaves):       100 - 3000 kg/ha
             *   WST (stems):        200 - 6000 kg/ha
             *   WSO (storage organs): 0 - 8000 kg/ha
             * ================================================================= */
            real WRT = rds.WRT(i);   /* Live root dry weight [kg/ha] */
            real WLV = lds.WLV(i);   /* Live leaf dry weight [kg/ha] */
            real WST = sds.WST(i);   /* Live stem dry weight [kg/ha] */
            real WSO = sods.WSO(i);  /* Live storage organ dry weight [kg/ha] */
            real DVS = s.DVS(i);     /* Development stage [0..2], unitless */

            /* ================================================================= *
             * Section 2: Retrieve meteorological driving variables
             *
             * Daily mean temperature is computed as the arithmetic mean of
             * daily minimum and maximum air temperature.
             * Units: degrees Celsius
             * Typical range: -10 to +40 deg C depending on climate
             * ================================================================= */
            real TMIN = m.tmin(drv_idx);     /* Daily minimum temperature [deg C] */
            real TMAX = m.tmax(drv_idx);     /* Daily maximum temperature [deg C] */
            real TEMP = (TMIN + TMAX) / 2.0; /* Daily mean temperature [deg C] */

            /* ================================================================= *
             * Section 3: Retrieve crop-specific maintenance respiration coefficients
             *
             * These parameters represent the maintenance respiration rate per unit
             * dry weight of each organ at the reference temperature of 25 deg C.
             * Units: [kg CH2O / (kg dry matter / day)]
             * Typical values (WOFOST defaults for wheat):
             *   RMR (roots):         0.010
             *   RML (leaves):        0.030
             *   RMS (stems):         0.015
             *   RMO (storage organs): 0.010
             *
             * Q10: Temperature sensitivity factor. Each 10-degree-Celsius increase
             *      above the reference temperature (25 deg C) multiplies the
             *      respiration rate by Q10. Typical value: 2.0.
             *      Range: 1.5 to 2.5 depending on crop species.
             * ================================================================= */
            real RMR = p.p.RMR;  /* Maintenance respiration coefficient for roots [kg CH2O kg^-1 d^-1] */
            real RML = p.p.RML;  /* Maintenance respiration coefficient for leaves [kg CH2O kg^-1 d^-1] */
            real RMS = p.p.RMS;  /* Maintenance respiration coefficient for stems [kg CH2O kg^-1 d^-1] */
            real RMO = p.p.RMO;  /* Maintenance respiration coefficient for storage organs [kg CH2O kg^-1 d^-1] */
            real Q10 = p.p.Q10;  /* Q10 temperature sensitivity factor [dimensionless] */

            /* ================================================================= *
             * Section 4: Compute base maintenance respiration (RMRES)
             *
             * The uncorrected maintenance respiration is the weighted sum of
             * organ biomasses multiplied by their respective respiration
             * coefficients. This represents the total carbohydrate cost of
             * maintaining all living tissue at the reference temperature.
             *
             * RMRES = RMR * WRT + RML * WLV + RMS * WST + RMO * WSO
             *
             * Units: [kg CH2O / ha / d]
             * ================================================================= */
            real RMRES = (RMR * WRT) +
                         (RML * WLV) +
                         (RMS * WST) +
                         (RMO * WSO);

            /* ================================================================= *
             * Section 5: Apply senescence reduction factor (RFSETB)
             *
             * As the crop matures (DVS > 1.0), leaves and other organs begin
             * to senesce, reducing the total living biomass and hence the
             * maintenance respiration demand. The RFSETB table provides a
             * dimensionless multiplier in [0, 1] as a function of DVS.
             *
             * Equivalent Python: RMRES *= p.RFSETB(kk["DVS"])
             * ================================================================= */
            real RFSETB_val = Afgen::lookup(RFSETB, DVS);
            RMRES *= RFSETB_val;

            /* ================================================================= *
             * Section 6: Temperature correction using the Q10 model (TEFF)
             *
             * The Q10 model assumes respiration increases exponentially with
             * temperature. The reference temperature is 25 deg C (at which
             * TEFF = 1.0). For every 10-degree-Celsius increase above 25 deg C,
             * respiration multiplies by Q10.
             *
             * TEFF = Q10 ^ ((TEMP - 25.0) / 10.0)
             *
             * Examples with Q10 = 2.0:
             *   TEMP = 5  deg C -> TEFF = 0.25  (75% reduction)
             *   TEMP = 15 deg C -> TEFF = 0.50  (50% reduction)
             *   TEMP = 25 deg C -> TEFF = 1.00  (reference)
             *   TEMP = 35 deg C -> TEFF = 2.00  (doubled)
             *
             * Units: [dimensionless]
             * ================================================================= */
            real TEFF = pow(Q10, (TEMP - 25.0) / 10.0);

            /* ================================================================= *
             * Section 7: Compute and store potential maintenance respiration (PMRES)
             *
             * The final potential maintenance respiration is the product of the
             * senescence-adjusted base rate and the temperature correction factor.
             *
             * PMRES = RMRES * TEFF
             *
             * Units: [kg CH2O / ha / d]
             *
             * This value is later used to partition available assimilates:
             *   - First, PMRES is subtracted from gross assimilation (PGASS)
             *   - The remainder, if positive, is available for growth
             *   - If PMRES exceeds PGASS, the deficit is taken from storage
             * ================================================================= */
            rs.PMRES(i) = RMRES * TEFF;

            /* ================================================================= *
             * Section 8: Debug output (conditional compilation)
             *
             * When DEBUG_CROP_GROWTH_MODEL is defined, detailed per-cell
             * diagnostics are printed for verification and troubleshooting.
             * ================================================================= */
#if DEBUG_CROP_GROWTH_MODEL
            printf("[RESPIRATION-DEBUG] cell=%d, DVS=%.4f, WRT=%.6f, WLV=%.6f, WST=%.6f, WSO=%.6f\n",
                   i, DVS, WRT, WLV, WST, WSO);
            printf("[RESPIRATION-DEBUG] cell=%d, RMR=%.6f, RML=%.6f, RMS=%.6f, RMO=%.6f\n",
                   i, RMR, RML, RMS, RMO);
            printf("[RESPIRATION-DEBUG] cell=%d, RMRES=%.6f, RFSETB=%.6f, TEMP=%.2f, Q10=%.2f, TEFF=%.6f\n",
                   i, RMRES / RFSETB_val, RFSETB_val, TEMP, Q10, TEFF);
            printf("[RESPIRATION-DEBUG] cell=%d, PMRES=%.6f\n", i, rs.PMRES(i));
#endif
        });
    }
};

#endif
