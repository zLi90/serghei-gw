/* -*- mode: c++; c-default-style: "linux" -*- */

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
    static constexpr real PI = 3.14159265358979323846;
    static constexpr real RAD = 0.017453292519943295; // PI / 180.0

public:
    // --- ASTRO Routine ---
    // Calculates astronomic daylength and diurnal radiation characteristics
    KOKKOS_INLINE_FUNCTION
    AstroVars astro(int day_of_year, real latitude, real radiation, bool debug) const
    {
        AstroVars av;

        real LAT = latitude; // In degrees
        real IDAY = (real)day_of_year;
        real AVRAD = radiation; // J/m2/day
        real ANGLE = -4.0;      // degrees

        // Declination and solar constant
        real sin_2345 = sin(23.45 * RAD);
        real cos_time = cos(2.0 * PI * (IDAY + 10.0) / 365.0);
        real DEC = -asin(sin_2345 * cos_time);

        real SC = 1370.0 * (1.0 + 0.033 * cos(2.0 * PI * IDAY / 365.0));

        // Intermediate variables
        real sin_lat = sin(RAD * LAT);
        real cos_lat = cos(RAD * LAT);
        real sin_dec = sin(DEC);
        real cos_dec = cos(DEC);

        av.SINLD = sin_lat * sin_dec;
        av.COSLD = cos_lat * cos_dec;

        real AOB = av.SINLD / av.COSLD;
        real DSINB = 0.0;

        if (std::abs(AOB) <= 1.0)
        {
            av.DAYL = 12.0 * (1.0 + 2.0 * asin(AOB) / PI);
            real sqrt_1_AOB2 = sqrt(1.0 - AOB * AOB);
            DSINB = 3600.0 * (av.DAYL * av.SINLD + 24.0 * av.COSLD * sqrt_1_AOB2 / PI);
            av.DSINBE = 3600.0 * (av.DAYL * (av.SINLD + 0.4 * (av.SINLD * av.SINLD + av.COSLD * av.COSLD * 0.5)) +
                                  12.0 * av.COSLD * (2.0 + 3.0 * 0.4 * av.SINLD) * sqrt_1_AOB2 / PI);
        }
        else
        {
            if (AOB > 1.0)
                av.DAYL = 24.0;
            if (AOB < -1.0)
                av.DAYL = 0.0;
            DSINB = 3600.0 * (av.DAYL * av.SINLD);
            av.DSINBE = 3600.0 * (av.DAYL * (av.SINLD + 0.4 * (av.SINLD * av.SINLD + av.COSLD * av.COSLD * 0.5)));
        }

        real sin_angle = sin(ANGLE * RAD);
        real AOB_CORR = (-sin_angle + av.SINLD) / av.COSLD;

        if (std::abs(AOB_CORR) <= 1.0)
        {
            av.DAYLP = 12.0 * (1.0 + 2.0 * asin(AOB_CORR) / PI);
        }
        else if (AOB_CORR > 1.0)
        {
            av.DAYLP = 24.0;
        }
        else
        {
            av.DAYLP = 0.0;
        }

        av.ANGOT = SC * DSINB;
        if (av.DAYL > 0.0 && av.ANGOT > 0.0)
        {
            av.ATMTR = AVRAD / av.ANGOT;
        }
        else
        {
            av.ATMTR = 0.0;
        }

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

        av.DIFPP = FRDIF * av.ATMTR * 0.5 * SC;

#if DEBUG_CROP_GROWTH_MODEL
        printf("[ASSIM-DEBUG-ASTRO] IDAY=%.1f, SC=%.2f, DEC=%.4f, SINLD=%.4f, COSLD=%.4f, DAYL=%.4f, DSINBE=%.4f, ANGOT=%.4f, ATMTR=%.4f, FRDIF=%.4f, DIFPP=%.4f\n",
               IDAY, SC, DEC, av.SINLD, av.COSLD, av.DAYL, av.DSINBE, av.ANGOT, av.ATMTR, FRDIF, av.DIFPP);
#endif

        return av;
    }

    // --- ASSIM7 Routine ---
    KOKKOS_INLINE_FUNCTION
    real assim7(real AMAX, real EFF, real LAI, real KDIF, real SINB, real PARDIR, real PARDIF, bool debug) const
    {
        const real XGAUSS[3] = {0.1127017, 0.5000000, 0.8872983};
        const real WGAUSS[3] = {0.2777778, 0.4444444, 0.2777778};
        const real SCV = 0.2;

        real sqrt_1_SCV = sqrt(1.0 - SCV);
        real REFH = (1.0 - sqrt_1_SCV) / (1.0 + sqrt_1_SCV);
        real REFS = REFH * 2.0 / (1.0 + 1.6 * SINB);

        real KDIRBL = (0.5 / max(0.001, SINB)) * KDIF / (0.8 * sqrt_1_SCV);
        real KDIRT = KDIRBL * sqrt_1_SCV;

        real FGROS = 0.0;
        for (int i = 0; i < 3; ++i)
        {
            real LAIC = LAI * XGAUSS[i];

            real exp_KDIF = exp(-KDIF * LAIC);
            real exp_KDIRT = exp(-KDIRT * LAIC);
            real exp_KDIRBL = exp(-KDIRBL * LAIC);

            real VISDF = (1.0 - REFS) * PARDIF * KDIF * exp_KDIF;
            real VIST = (1.0 - REFS) * PARDIR * KDIRT * exp_KDIRT;
            real VISD = (1.0 - SCV) * PARDIR * KDIRBL * exp_KDIRBL;

            real VISSHD = VISDF + VIST - VISD;

            real AMAX_eff = max(2.0, AMAX); // Check this logic against Python
            real FGRSH = AMAX * (1.0 - exp(-VISSHD * EFF / AMAX_eff));

            real VISPP = (1.0 - SCV) * PARDIR / max(0.001, SINB);

            real FGRSUN;
            if (VISPP <= 0.0)
            {
                FGRSUN = FGRSH;
            }
            else
            {
                FGRSUN = AMAX * (1.0 - (AMAX - FGRSH) * (1.0 - exp(-VISPP * EFF / AMAX_eff)) / (EFF * VISPP));
            }

            real FSLLA = exp_KDIRBL;
            real FGL = FSLLA * FGRSUN + (1.0 - FSLLA) * FGRSH;

            FGROS += FGL * WGAUSS[i];

#if DEBUG_CROP_GROWTH_MODEL
            printf("[ASSIM-DEBUG-ASSIM7] i=%d, LAIC=%.6f, VISDF=%.6f, VIST=%.6f, VISD=%.6f, VISSHD=%.6f, FGRSH=%.6f, FGRSUN=%.6f, FSLLA=%.6f, FGL=%.6f\n",
                   i, LAIC, VISDF, VIST, VISD, VISSHD, FGRSH, FGRSUN, FSLLA, FGL);
#endif
        }

        FGROS = FGROS * LAI;
        return FGROS;
    }

    // --- TOTASS7 Routine ---
    KOKKOS_INLINE_FUNCTION
    real totass7(real DAYL, real AMAX, real EFF, real LAI, real KDIF, real AVRAD,
                 real DIFPP, real DSINBE, real SINLD, real COSLD, bool debug) const
    {

        const real XGAUSS[3] = {0.1127017, 0.5000000, 0.8872983};
        const real WGAUSS[3] = {0.2777778, 0.4444444, 0.2777778};

        real DTGA = 0.0;

        if (AMAX > 0.0 && LAI > 0.0 && DAYL > 0.0)
        {
            for (int i = 0; i < 3; ++i)
            {
                real HOUR = 12.0 + 0.5 * DAYL * XGAUSS[i];
                real cos_term = cos(2.0 * PI * (HOUR + 12.0) / 24.0);
                real SINB = max(0.0, SINLD + COSLD * cos_term);

                real safe_DSINBE = max(1.0, DSINBE);
                real PAR = 0.5 * AVRAD * SINB * (1.0 + 0.4 * SINB) / safe_DSINBE;

                real PARDIF = min(PAR, SINB * DIFPP);
                real PARDIR = PAR - PARDIF;

#if DEBUG_CROP_GROWTH_MODEL
                printf("[ASSIM-DEBUG-TOTASS7] i=%d, HOUR=%.4f, SINB=%.6f, PAR=%.6f, PARDIF=%.6f, PARDIR=%.6f\n",
                       i, HOUR, SINB, PAR, PARDIF, PARDIR);
#endif

                real FGROS = assim7(AMAX, EFF, LAI, KDIF, SINB, PARDIR, PARDIF, debug);
                DTGA += FGROS * WGAUSS[i];

#if DEBUG_CROP_GROWTH_MODEL
                printf("[ASSIM-DEBUG-TOTASS7] i=%d, FGROS=%.6f, Accumulated DTGA=%.6f\n", i, FGROS, DTGA);
#endif
            }
        }

        DTGA *= DAYL;
        return DTGA;
    }

    // --- Main Calc Rates Function ---
    // --- Main Calc Rates Function ---
    void calc_rates(AssimilateState &as, const CropState &p, const CropDynamicState &s,
                    const MeteoState &m, int drv_idx, int current_day_of_year)
    {

        auto AMAXTB = p.tables.at("AMAXTB");
        auto EFFTB = p.tables.at("EFFTB");
        auto KDIFTB = p.tables.at("KDIFTB");
        auto TMPFTB = p.tables.at("TMPFTB");
        auto TMNFTB = p.tables.at("TMNFTB");

        Kokkos::parallel_for("Assim_Rates", as.nCells, KOKKOS_LAMBDA(const int i) {
            bool debug = (i == 0); // Only debug first cell

            real DVS = s.DVS(i);
            real LAI = s.LAI(i);

            real TMIN = m.tmin(drv_idx);
            real TMAX = m.tmax(drv_idx);

            // [MODIFIED]: Define both Daily Average and Daytime Average variables
            real TEMP = (TMIN + TMAX) / 2.0;          // Daily Average (24h)
            real DTEMP = TMIN + (TMAX - TMIN) * 0.75; // Daytime Average (Photosynthesis active period)

            real IRRAD = m.irrad(drv_idx);
            real LAT = m.lat;

            real TMINRA = TMIN; // Note: Strictly this should be 7-day rolling avg of Tmin, currently simplified

            AstroVars av = astro(current_day_of_year, LAT, IRRAD, debug);

            real AMAX_BASE = Afgen::lookup(AMAXTB, DVS); // Look up base

            // [MODIFIED]: Use DTEMP instead of TEMP for Assimilation parameters
            // Because photosynthesis happens during the day
            real TMPF = Afgen::lookup(TMPFTB, DTEMP);
            real AMAX = AMAX_BASE * TMPF;

            real KDIF = Afgen::lookup(KDIFTB, DVS);

            // [MODIFIED]: Use DTEMP for Light Use Efficiency
            real EFF = Afgen::lookup(EFFTB, DTEMP);

#if DEBUG_CROP_GROWTH_MODEL
            // [MODIFIED]: Updated debug print to show DTEMP
            printf("[ASSIM-DEBUG] cell=%d, DVS=%.4f, LAI=%.6f, TEMP(24h)=%.2f, DTEMP(Day)=%.2f, TMIN=%.2f, IRRAD=%.1f, LAT=%.2f\n",
                   i, DVS, LAI, TEMP, DTEMP, TMIN, IRRAD, LAT);
            printf("[ASSIM-DEBUG] AMAX_BASE=%.4f, TMPF=%.4f, AMAX=%.4f, KDIF=%.4f, EFF=%.4f\n", AMAX_BASE, TMPF, AMAX, KDIF, EFF);
            printf("[ASSIM-DEBUG] DAYL=%.4f, SINLD=%.4f, COSLD=%.4f, DIFPP=%.4f, DSINBE=%.4f\n",
                   av.DAYL, av.SINLD, av.COSLD, av.DIFPP, av.DSINBE);
#endif

            real DTGA = totass7(av.DAYL, AMAX, EFF, LAI, KDIF, IRRAD, av.DIFPP, av.DSINBE, av.SINLD, av.COSLD, debug);

            real TMNF = Afgen::lookup(TMNFTB, TMINRA);
            DTGA *= TMNF;

            as.PGASS(i) = DTGA * 30.0 / 44.0;

#if DEBUG_CROP_GROWTH_MODEL
            printf("[ASSIM-DEBUG] DTGA_FINAL=%.6f, TMNF=%.4f, PGASS=%.6f\n", DTGA, TMNF, as.PGASS(i));
#endif
        });
    }
};

#endif