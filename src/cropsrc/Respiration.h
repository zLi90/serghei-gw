/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _RESPIRATION_H_
#define _RESPIRATION_H_

#include "../define.h"
#include "CropState.h"        // Static parameters
#include "CropDynamicState.h" // Dynamic variables (Biomass, DVS)
#include "RespirationState.h" // Output variables
#include "MeteoState.h"       // Driving variables (Temp)
#include "Afgen.h"            // Table lookup helper
#include <cmath>

class Respiration
{

public:
       // Calculates maintenance respiration rate
       // drv_idx: Index in the meteo data for the current day
       void calc_rates(RespirationState &rs, const CropState &p, const CropDynamicState &s,
                       const MeteoState &m, int drv_idx)
       {

              // Access Table: Reduction factor for senescence
              // 对应 Python: p.RFSETB(kk["DVS"])
              auto RFSETB = p.tables.at("RFSETB");

              Kokkos::parallel_for("Respiration_Rates", rs.nCells, KOKKOS_LAMBDA(const int i) {
                     // 1. 获取状态变量 (从 CropDynamicState)
                     // 假设 s 中包含各器官干重 [kg ha-1] 和发育阶段
                     real WRT = s.WRT(i); // Dry weight of living roots
                     real WLV = s.WLV(i); // Dry weight of living leaves
                     real WST = s.WST(i); // Dry weight of living stems
                     real WSO = s.WSO(i); // Dry weight of living storage organs
                     real DVS = s.DVS(i); // Development stage

                     // 2. 获取气象变量
                     real TMIN = m.tmin(drv_idx);
                     real TMAX = m.tmax(drv_idx);
                     real TEMP = (TMIN + TMAX) / 2.0; // Daily mean temperature

                     // 3. 获取参数 (从 CropState)
                     // 维护呼吸系数 [kg CH2O kg-1 d-1]
                     real RMR = p.p.RMR; // Roots
                     real RML = p.p.RML; // Leaves
                     real RMS = p.p.RMS; // Stems
                     real RMO = p.p.RMO; // Storage Organs
                     real Q10 = p.p.Q10; // Q10 factor

                     // 4. 计算相对维持呼吸 (Relative Maintenance Respiration)
                     // RMRES = (p.RMR * kk["WRT"] + p.RML * kk["WLV"] + ...)
                     real RMRES = (RMR * WRT) +
                                  (RML * WLV) +
                                  (RMS * WST) +
                                  (RMO * WSO);

                     // 5. 针对衰老进行修正 (Correction for senescence)
                     // RMRES *= p.RFSETB(kk["DVS"])
                     real RFSETB_val = Afgen::lookup(RFSETB, DVS);
                     RMRES *= RFSETB_val;

                     // 6. 温度效应修正 (Temperature effect)
                     // TEFF = p.Q10**((drv.TEMP-25.)/10.)
                     // 对应 C++: pow(base, exponent)
                     real TEFF = pow(Q10, (TEMP - 25.0) / 10.0);

                     // 7. 计算潜在维持呼吸并存储 (Potential maintenance respiration)
                     // self.rates.PMRES = RMRES * TEFF
                     rs.PMRES(i) = RMRES * TEFF;

// 8. 添加调试输出
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
