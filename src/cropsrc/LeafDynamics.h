/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _LEAF_DYNAMICS_H_
#define _LEAF_DYNAMICS_H_

#include "../define.h"
#include "CropState.h"               // Static parameters
#include "CropDynamicState.h"        // Dynamic variables (DVS, ADMI, FL, SAI, PAI, RFTRA)
#include "LeafDynamicsState.h"       // Output variables
#include "EvapotranspirationState.h" // For RFTRA
#include "MeteoState.h"              // Driving variables (Temp)
#include "Afgen.h"                   // Table lookup helper
#include <cmath>
#include <algorithm>
#include <cstdio>

class LeafDynamics
{

public:
    // ... (Helper functions and initialize unchanged) ...
    // Note: initialize prints also adjusted to match if needed, but calc_rates/integrate are key.
    void initialize(LeafDynamicsState &lds, const CropState &p, const CropDynamicState &s,
                    const realArr &FR, const realArr &FL,
                    const realArr &SAI, const realArr &PAI)
    { // Needed for LAI calc

        auto SLATB = p.tables.at("SLATB");
        auto KDIFTB = p.tables.at("KDIFTB");

        // Access Tables
        auto FRTB = p.tables.at("FRTB");
        auto FLTB = p.tables.at("FLTB");
        auto FSTB = p.tables.at("FSTB");
        auto FOTB = p.tables.at("FOTB");

        Kokkos::parallel_for("LeafDynamics_Init", lds.nCells, KOKKOS_LAMBDA(const int i) {
            real TDWI = p.p.TDWI;
            // real fl = FL(i); // 使用传入的 FL
            // real fr = FR(i); // 使用传入的 FR
            real DVS = s.DVS(i);

            s.FR(i) = Afgen::lookup(FRTB, DVS);
            s.FL(i) = Afgen::lookup(FLTB, DVS);
            // s.FS(i) = Afgen::lookup(FSTB, DVS);
            // s.FO(i) = Afgen::lookup(FOTB, DVS);
#if DEBUG_CROP_GROWTH_MODEL
            if (i == 0)
            {
                printf("[LEAF-INIT-KERNEL] TDWI=%.4f, FL=%.4f, FR=%.4f, DVS=%.4f\n", TDWI, s.FL(i), s.FR(i), DVS);
            }
#endif
            // Initial leaf biomass
            real wlv = (TDWI * (1.0 - s.FR(i))) * s.FL(i);

            // --- 强制修正：使用 LAIEM ---
            // WOFOST 标准做法：如果算出来的 LAI 小于 LAIEM，强制使用 LAIEM
            real sla_val = Afgen::lookup(SLATB, DVS);
            real calc_lai = wlv * sla_val;
            real LAIEM = p.p.LAIEM; // 0.07 from input

            if (calc_lai < LAIEM)
            {
                // 如果计算值小于预设出苗 LAI，强制提升 WLV
                wlv = LAIEM / sla_val;
                calc_lai = LAIEM;
            }

            lds.WLV(i) = wlv;
            lds.DWLV(i) = 0.0;
            lds.TWLV(i) = wlv;

            // First leaf class
            lds.num_leaf_classes(i) = 1;
            lds.LV(i, 0) = wlv;
            lds.SLA(i, 0) = sla_val;
            lds.LVAGE(i, 0) = 0.0;

            lds.LAIEM(i) = calc_lai;
            lds.LASUM(i) = calc_lai;
            lds.LAIEXP(i) = calc_lai;
            lds.LAIMAX(i) = calc_lai;

            real sai = SAI(i);
            real pai = PAI(i);
            lds.LAI(i) = lds.LASUM(i) + sai + pai;
#if DEBUG_CROP_GROWTH_MODEL
            printf("[LEAF-INIT-KERNEL] Final Init LAI=%.4f (LASUM=%.4f, SAI=%.4f, PAI=%.4f)\n",
                   lds.LAI(i), lds.LASUM(i), sai, pai);
#endif
        });
    }

    // Calculate rates
    void calc_rates(LeafDynamicsState &lds, const CropState &p, const CropDynamicState &s,
                    const MeteoState &m, int drv_idx,
                    const realArr &ADMI, const EvapotranspirationState &ets)
    {

        auto KDIFTB = p.tables.at("KDIFTB");
        auto SLATB = p.tables.at("SLATB");
        // Access Tables
        auto FRTB = p.tables.at("FRTB");
        auto FLTB = p.tables.at("FLTB");
        auto FSTB = p.tables.at("FSTB");
        auto FOTB = p.tables.at("FOTB");

        Kokkos::parallel_for("LeafDynamics_Rates", lds.nCells, KOKKOS_LAMBDA(const int i) {
            real DVS = s.DVS(i);
            real TMIN = m.tmin(drv_idx);
            real TMAX = m.tmax(drv_idx);
            real TEMP = (TMIN + TMAX) / 2.0;

            // Growth rate leaves
            // -----------------------------------------------------------
            // [FIXED]: Removed hardcoded fl. Now using dynamic FL from CropDynamicState
            // which corresponds to the FL calculated in Partitioning.h based on DVS tables.
            // -----------------------------------------------------------

            s.FR(i) = Afgen::lookup(FRTB, DVS);
            s.FL(i) = Afgen::lookup(FLTB, DVS);
            // s.FS(i) = Afgen::lookup(FSTB, DVS);
            // s.FO(i) = Afgen::lookup(FOTB, DVS);

            lds.GRLV(i) = ADMI(i) * s.FL(i);

            // Death due to water stress
            lds.DSLV1(i) = lds.WLV(i) * (1.0 - ets.RFTRA(i)) * p.p.PERDL;

            // Death due to self shading
            real kdif = Afgen::lookup(KDIFTB, DVS);
            real LAICR = 3.2 / kdif;
            real term = 0.03 * (lds.LAI(i) - LAICR) / LAICR;
            if (term < 0.0)
                term = 0.0;
            if (term > 0.03)
                term = 0.03;
            lds.DSLV2(i) = lds.WLV(i) * term;

            // Death due to frost
            lds.DSLV3(i) = 0.0;

            // Max death rate
            real dslv = fmax(lds.DSLV1(i), fmax(lds.DSLV2(i), lds.DSLV3(i)));
            lds.DSLV(i) = dslv;

            // Death due to aging (DALV)
            real dalv_sum = 0.0;
            real SPAN = p.p.SPAN;
            int n_classes = lds.num_leaf_classes(i);

            for (int c = 0; c < n_classes; ++c)
            {
                if (lds.LVAGE(i, c) > SPAN)
                {
                    dalv_sum += lds.LV(i, c);
                }
            }
            lds.DALV(i) = dalv_sum;

            // Total death rate
            lds.DRLV(i) = fmax(dslv, dalv_sum);

            // Physiologic ageing
            real TBASE = p.p.TBASE;
            real fysage = (TEMP - TBASE) / (35.0 - TBASE);
            if (fysage < 0.0)
                fysage = 0.0;
            lds.FYSAGE(i) = fysage;

            // Specific leaf area
            lds.SLAT(i) = Afgen::lookup(SLATB, DVS);

            // Leaf area expansion
            real gla = 0.0; // For debug print
            if (lds.LAIEXP(i) < 6.0)
            {
                real dteff = TEMP - TBASE;
                if (dteff < 0.0)
                    dteff = 0.0;

                lds.GLAIEX(i) = lds.LAIEXP(i) * p.p.RGRLAI * dteff; // Potential
                lds.GLASOL(i) = lds.GRLV(i) * lds.SLAT(i);          // Source Limited

                // Actual growth rate is minimum of potential and source-limited
                gla = fmin(lds.GLAIEX(i), lds.GLASOL(i));

                // [FIX]: Update GLAIEX to be the ACTUAL rate, so LAIEXP updates correctly in integrate()
                lds.GLAIEX(i) = gla;

                if (lds.GRLV(i) > 0.0)
                {
                    // Adjust specific leaf area for new leaves
                    lds.SLAT(i) = gla / lds.GRLV(i);
                }

// DEBUG: Leaf area expansion (Match PY-LEAF-RATES-EXPANSION)
#if DEBUG_CROP_GROWTH_MODEL
                if (i == 0)
                {
                    printf("[CPP-LEAF-RATES-EXPANSION] cell=%d, LAIEXP=%.6f, GLAIEX=%.6f, GLASOL=%.6f, GLA=%.6f, SLAT=%.6f\n",
                           i, lds.LAIEXP(i), lds.GLAIEX(i), lds.GLASOL(i), gla, lds.SLAT(i));
                }
#endif
            }
            else
            {
                lds.GLAIEX(i) = 0.0;
                lds.GLASOL(i) = 0.0;
            }

// DEBUG: General Rates (Match Python Prints)
#if DEBUG_CROP_GROWTH_MODEL
            if (i == 0)
            {
                printf("[CPP-LEAF-RATES-DETAIL] cell=%d, FL=%.6f, ADMI=%.6f, GRLV=%.6f, WLV=%.6f\n",
                       i, s.FL(i), ADMI(i), lds.GRLV(i), lds.WLV(i));
                printf("[CPP-LEAF-RATES-DETAIL] cell=%d, DVS=%.4f, RFTRA=%.6f\n",
                       i, s.DVS(i), ets.RFTRA(i));
                printf("[CPP-LEAF-RATES-DEATH] cell=%d, DSLV1=%.6f, DSLV2=%.6f, DALV=%.6f, DRLV=%.6f\n",
                       i, lds.DSLV1(i), lds.DSLV2(i), lds.DALV(i), lds.DRLV(i));
                printf("[CPP-LEAF-RATES-AGE] cell=%d, TEMP=%.2f, TBASE=%.2f, FYSAGE=%.6f\n",
                       i, TEMP, TBASE, lds.FYSAGE(i));
            }
#endif
        });
    }

    // Integrate states
    void integrate(LeafDynamicsState &lds, const CropDynamicState &s, real delt = 1.0)
    {

        Kokkos::parallel_for("LeafDynamics_Integrate", lds.nCells, KOKKOS_LAMBDA(const int i) {

// DEBUG: Before integration (Match PY-LEAF-INTEGRATE-START)
#if DEBUG_CROP_GROWTH_MODEL
            if (i == 0)
            {
                printf("[CPP-LEAF-INTEGRATE-START] cell=%d, WLV=%.6f, LASUM=%.6f, LAI=%.6f, GRLV=%.6f, DRLV=%.6f\n",
                       i, lds.WLV(i), lds.LASUM(i), lds.LAI(i), lds.GRLV(i), lds.DRLV(i));
            }
#endif

            real tDRLV = lds.DRLV(i) * delt;
            int n_classes = lds.num_leaf_classes(i);

// DEBUG: Death amount (Match PY-LEAF-INTEGRATE-DEATH)
#if DEBUG_CROP_GROWTH_MODEL
            if (i == 0)
            {
                printf("[CPP-LEAF-INTEGRATE-DEATH] cell=%d, tDRLV=%.6f, n_classes=%d\n",
                       i, tDRLV, n_classes);
            }
#endif

            // Iterate from Oldest to Newest
            while (n_classes > 0 && tDRLV > 0.0)
            {
                int oldest_idx = n_classes - 1;
                real lv_weight = lds.LV(i, oldest_idx);

                if (tDRLV >= lv_weight)
                {
                    tDRLV -= lv_weight;
                    n_classes--;
                }
                else
                {
                    lds.LV(i, oldest_idx) -= tDRLV;
                    tDRLV = 0.0;
                }
            }
            lds.num_leaf_classes(i) = n_classes;

            // Physiological Ageing
            real fysage = lds.FYSAGE(i) * delt;
            for (int c = 0; c < n_classes; ++c)
            {
                lds.LVAGE(i, c) += fysage;
            }

            // Leaf Growth
            if (n_classes < LeafDynamicsState::MAX_LEAF_CLASSES)
            {
                for (int c = n_classes; c > 0; --c)
                {
                    lds.LV(i, c) = lds.LV(i, c - 1);
                    lds.SLA(i, c) = lds.SLA(i, c - 1);
                    lds.LVAGE(i, c) = lds.LVAGE(i, c - 1);
                }
                lds.LV(i, 0) = lds.GRLV(i) * delt;
                lds.SLA(i, 0) = lds.SLAT(i);
                lds.LVAGE(i, 0) = 0.0;
                lds.num_leaf_classes(i)++;
            }

            // Update Aggregates
            real sum_lv = 0.0;
            real sum_lasum = 0.0;
            int current_n = lds.num_leaf_classes(i);

            for (int c = 0; c < current_n; ++c)
            {
                sum_lv += lds.LV(i, c);
                sum_lasum += lds.LV(i, c) * lds.SLA(i, c);
            }

            lds.WLV(i) = sum_lv;
            lds.LASUM(i) = sum_lasum;

            lds.DWLV(i) += lds.DRLV(i) * delt;
            lds.TWLV(i) = lds.WLV(i) + lds.DWLV(i);

            real sai = s.SAI(i);
            real pai = s.PAI(i);
            real old_lai = lds.LAI(i);
            lds.LAI(i) = lds.LASUM(i) + sai + pai;

            if (lds.LAI(i) > lds.LAIMAX(i))
            {
                lds.LAIMAX(i) = lds.LAI(i);
            }

            // Update LAIEXP with the ACTUAL growth rate (GLAIEX was updated in calc_rates)
            lds.LAIEXP(i) += lds.GLAIEX(i) * delt;

// DEBUG: After integration (Match PY-LEAF-INTEGRATE-END)
#if DEBUG_CROP_GROWTH_MODEL
            if (i == 0)
            {
                printf("[CPP-LEAF-INTEGRATE-END] cell=%d, WLV=%.6f, LASUM=%.6f, LAI=%.6f (old=%.6f), LAIEXP=%.6f\n",
                       i, lds.WLV(i), lds.LASUM(i), lds.LAI(i), old_lai, lds.LAIEXP(i));
                printf("[CPP-LEAF-INTEGRATE-END] cell=%d, SAI=%.6f, PAI=%.6f, LAIMAX=%.6f\n",
                       i, sai, pai, lds.LAIMAX(i));
            }
#endif
        });
    }
};

#endif