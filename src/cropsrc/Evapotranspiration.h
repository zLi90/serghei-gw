/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _EVAPOTRANSPIRATION_H_
#define _EVAPOTRANSPIRATION_H_

#include "../define.h"
#include "CropState.h"               // Static crop parameters
#include "CropDynamicState.h"        // DVS, LAI, SM (from SERGHEI)
#include "MeteoState.h"              // TMIN, TMAX, IRRAD, VAP, WIND (for ET0 calculation)
#include "EvapotranspirationState.h" // Output TRA, TRAMX, EVS, RFTRA
#include "Afgen.h"                   // For KDIFTB
#include <cmath>
#include <algorithm> // For fmin, fmax

// Helper function from evapotranspiration.py
KOKKOS_INLINE_FUNCTION
real SWEAF(real ET0_val, real DEPNR)
{
    const real A = 0.76;
    const real B = 1.5;
    real sweaf = 1.0 / (A + B * ET0_val) - (5.0 - DEPNR) * 0.10;

    if (DEPNR < 3.0)
    {
        sweaf += (ET0_val - 0.6) / (DEPNR * (DEPNR + 3.0));
    }
    // Limit to [0.10, 0.95]
    if (sweaf < 0.10)
        sweaf = 0.10;
    if (sweaf > 0.95)
        sweaf = 0.95;
    return sweaf;
}

class Evapotranspiration
{

public:
    // Calculates potential and actual evaporation/transpiration rates
    // root_zone_SM: Average volumetric soil moisture content in the root zone (from SERGHEI)
    // The E0, ES0, ET0 values are now calculated in MeteoInit using proper Penman and Penman-Monteith methods
    void calc_rates(EvapotranspirationState &ets, const CropState &p, const CropDynamicState &s,
                    const MeteoState &m, int drv_idx, const realArr &root_zone_SM)
    {

        // Access Tables
        auto KDIFTB = p.tables.at("KDIFTB");

        Kokkos::parallel_for("EvapoTrans_Rates", ets.nCells, KOKKOS_LAMBDA(const int i) {

            
            real DVS = s.DVS(i);
            real LAI = s.LAI(i); // Assuming LAI from CropDynamicState
            // real LAI = 1;
            real SM = root_zone_SM(i); // From SERGHEI

            // 1. Get E0, ES0, ET0 from MeteoState (already calculated in MeteoInit)
            // Units: cm/day (converted from mm/day in MeteoInit)
            ets.E0(i) = m.e0(drv_idx);   // cm/day - Potential evaporation from open water
            ets.ES0(i) = m.es0(drv_idx); // cm/day - Potential evaporation from bare soil
            ets.ET0(i) = m.et0(drv_idx); // cm/day - Potential evapotranspiration from reference crop

            // 2. Crop specific correction on potential transpiration rate
            // CFET is a crop-specific factor for ET0
            real ET0_CROP = fmax(0.0, p.p.CFET * ets.ET0(i));

            // 3. Maximum evaporation and transpiration rates
            // KGLOB: light extinction coefficient for the canopy
            real KGLOB = 0.75 * Afgen::lookup(KDIFTB, DVS);
            real EKL = exp(-KGLOB * LAI);

            // EVS: Maximum evaporation from wet soil below canopy
            // This decreases exponentially with LAI due to canopy shading
            ets.EVS(i) = ets.E0(i) * EKL;
            
            // TRAMX: Maximum transpiration from crop canopy
            // This is the portion of ET0 that reaches the canopy (1 - EKL)
            ets.TRAMX(i) = ET0_CROP * (1.0 - EKL);
#if DEBUG_CROP_GROWTH_MODEL
            printf("Evapotranspiration: i=%d, DVS=%.2f, LAI=%.2f, SM=%.2f, E0=%.6f, ES0=%.6f, ET0=%.6f, EVS=%.6f, TRAMX=%.6f\n",
                   i, DVS, LAI, SM, ets.E0(i), ets.ES0(i), ets.ET0(i), ets.EVS(i), ets.TRAMX(i));
#endif
       
            

            // 4. Critical soil moisture for water stress calculation
            // Note: SMFCF, SMW should come from SoilState, but for now using placeholders
            // In WOFOST, these are soil parameters that define:
            // - SMFCF: Field capacity
            // - SMW: Wilting point
            // - DEPNR: Rooting depth factor
            
            // TODO: Replace these with actual soil parameters from SoilState
            real SMFCF = 0.30; // Placeholder: Field capacity
            real SMW = 0.10;   // Placeholder: Wilting point
            real DEPNR = p.p.DEPNR; // Rooting depth factor from crop parameters

            // Calculate critical soil moisture content for water stress
            real SWDEP = SWEAF(ET0_CROP, DEPNR);
            real SMCR = (1.0 - SWDEP) * (SMFCF - SMW) + SMW;

            // 5. Reduction factor for transpiration in case of water shortage (RFWS)
            // This reduces transpiration when soil moisture is between wilting point and critical point
            real RFWS = (SM - SMW) / (SMCR - SMW);
            if (RFWS < 0.0) RFWS = 0.0;
            if (RFWS > 1.0) RFWS = 1.0;

            // 6. Reduction factor for oxygen shortage (RFOS)
            // Only applied if oxygen stress is enabled (IOX = 1)
            real RFOS = 1.0;
            if (p.p.IOX == 1) {
                // Oxygen stress occurs when soil is saturated
                // This requires SM0 (saturated moisture) and CRAIRC (critical air content)
                // TODO: Get these from SoilState
                real SM0 = 0.45; // Placeholder: Saturated moisture content
                real CRAIRC = 0.06; // Placeholder: Critical air content
                
                // RFOS reduces transpiration when soil is too wet (oxygen deficiency)
                real RFOSMX = (SM0 - SM) / CRAIRC;
                RFOS = std::max(0.0, std::min(1.0, RFOSMX));
            }

            // 7. Total reduction factor for transpiration
            ets.RFTRA(i) = RFOS * RFWS;

            // 8. Actual transpiration rate
            ets.TRA(i) = ets.TRAMX(i) * ets.RFTRA(i); });
    }
};

#endif
