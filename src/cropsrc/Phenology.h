#ifndef _PHENOLOGY_H_
#define _PHENOLOGY_H_

#include "../define.h"
#include "CropState.h"        // 静态参数（如温度阈值、发育阶段阈值）
#include "CropDynamicState.h" // 动态状态（如当前发育阶段、DVS值、累积积温）
#include "MeteoState.h"       // 气象数据（日最低/最高温）
#include "Afgen.h"            // 查表工具（用于温度-发育速率映射）
#include <cmath>

class Phenology {
public:
    // 初始化作物发育状态（播种或出苗起始）
    void initialize(CropDynamicState &s, const CropState &p) {
        // 【关键逻辑】根据初始发育状态（DVSI）决定起始阶段
        // DVSI < 0: 从播种（DVS=-0.1）开始；DVSI >=0: 从出苗（DVS=0.0）开始
        Kokkos::parallel_for("Pheno_Init", s.nCells, KOKKOS_LAMBDA(const int i) {
            if (p.p.DVSI < 0.0) {
                s.STAGE(i) = CropStage::EMERGING;  // 起始阶段：出苗（EMERGING）
                s.DVS(i) = p.p.DVSI;               // 例如：-0.1（播种时DVS）
            } else {
                s.STAGE(i) = CropStage::VEGETATIVE; // 起始阶段：营养生长期
                s.DVS(i) = p.p.DVSI;               // 例如：0.0（已出苗）
            }

            // 初始化累积积温（TSUM）和春化量（VERN）
            s.TSUM(i) = 0.0;   // 用于营养生长期的累积积温
            s.TSUME(i) = 0.0;  // 用于出苗阶段的累积积温（代码中未使用，可能为冗余）
            s.VERN(i) = 0.0;   // 春化量（Vernalisation accumulation）
            s.ISVERNALISED(i) = 0; // 春化完成标志（0=未完成，1=已完成）
        });
    }

    // 计算每日发育速率（DVR）和春化速率（VERNR）
    void calc_rates(CropDynamicState &s, const CropState &p, const MeteoState &m, int drv_idx, int current_day_of_year) {
        // 获取发育速率查表（DTSMTB）和春化速率查表（VERNRTB）
        auto DTSMTB = p.tables.at("DTSMTB");  // 温度-发育速率映射表
        auto VERNRTB = p.hasTable("VERNRTB") ? p.tables.at("VERNRTB") : realArr(); // 春化速率表（若不存在则用空表）

        Kokkos::parallel_for("Pheno_Rates", s.nCells, KOKKOS_LAMBDA(const int i) {
            // 1. 获取日均温
            real TMIN = m.tmin(drv_idx);  // 日最低温
            real TMAX = m.tmax(drv_idx);  // 日最高温
            real TEMP = (TMIN + TMAX) / 2.0; // 日均温（℃）

            // 2. 【春化逻辑】仅当作物类型支持春化（IDSL>=2）且处于营养生长期
            real VERNFAC = 1.0; // 春化速率修正因子（默认1.0）
            if (p.p.IDSL >= 2) { // IDSL=2表示需春化（如冬小麦）
                if (s.STAGE(i) == CropStage::VEGETATIVE && s.ISVERNALISED(i) == 0) {
                    // 【关键】春化条件：DVS < 临界值（VERNDVS=0.3）
                    real VERNDVS = 0.3; // 临界DVS（春化启动阈值，需从CropParams读取）
                    
                    if (s.DVS(i) < VERNDVS) {
                        // 从春化速率表查表：VERNRTB(TEMP) → 获取春化速率
                        real rate = Afgen::lookup(VERNRTB, TEMP);
                        s.VERNR(i) = rate; // 每日春化速率

                        // 【春化修正因子】根据累积春化量调整速率
                        real req_base = p.p.VERNBASE; // 春化启动阈值（需补充参数）
                        real req_sat = p.p.VERNSAT;   // 春化饱和阈值（需补充参数）
                        real r = (s.VERN(i) - req_base) / (req_sat - req_base);
                        r = std::max(0.0, std::min(1.0, r)); // 限制在[0,1]区间
                        VERNFAC = r; // 修正因子（累积春化量不足时速率降低）
                    } else {
                        // 【关键】达到临界DVS，强制春化完成
                        s.VERNR(i) = 0.0; // 停止春化
                        VERNFAC = 1.0;
                    }
                } else {
                    s.VERNR(i) = 0.0; // 非春化阶段或已完成春化
                    VERNFAC = 1.0;
                }
            }

            // 3. 【光周期敏感逻辑】仅当IDSL>=1（需光周期敏感）
            real DVRED = 1.0; // 光周期修正因子（默认1.0）
            if (p.p.IDSL >= 1) {
                // 调用 calculate_daylength 计算日长
                // 纬度从 MeteoState 中获取
                real DAYLP = this->calculate_daylength(current_day_of_year, m.lat);
                
                // 计算光周期对发育速率的折减因子
                // (DAYLP - DLC) / (DLO - DLC)
                // DLO: Optimum daylength
                // DLC: Critical daylength
                real DLO = p.p.DLO;
                real DLC = p.p.DLC;
                
                real ratio = (DAYLP - DLC) / (DLO - DLC);
                
                // Limit to [0, 1]
                if (ratio < 0.0) ratio = 0.0;
                if (ratio > 1.0) ratio = 1.0;
                
                DVRED = ratio;
            }

            // 4. 【核心】根据当前阶段计算发育速率（DVR）
            if (s.STAGE(i) == CropStage::EMERGING) {
                // 出苗阶段：DVS从-0.1→0.0（需累积TSUMEM）
                real dtsume = TEMP - p.p.TBASEM; // 有效温度（TBASEM=出苗温度基点）
                if (dtsume < 0.0) dtsume = 0.0; // 低于基点无发育
                // 【关键】DVS变化量 = 0.1 * (有效温度 / TSUMEM)
                // （因DVS从-0.1→0.0需0.1增量，TSUMEM为达到0.0所需总积温）
                s.DVR(i) = 0.1 * dtsume / p.p.TSUMEM; 
            } 
            else if (s.STAGE(i) == CropStage::VEGETATIVE) {
                // 营养生长期：DVS从0.0→1.0
                s.DTSUM(i) = Afgen::lookup(DTSMTB, TEMP) * VERNFAC * DVRED; // 有效积温（含春化/光周期修正）
                s.DVR(i) = s.DTSUM(i) / p.p.TSUM1; // 发育速率 = 有效积温 / 达到DVS=1所需积温
            }
            else if (s.STAGE(i) == CropStage::REPRODUCTIVE) {
                // 生殖生长期：DVS从1.0→DVSEND（通常=2.0）
                s.DTSUM(i) = Afgen::lookup(DTSMTB, TEMP); // 有效积温
                s.DVR(i) = s.DTSUM(i) / p.p.TSUM2; // 发育速率
            }
            else if (s.STAGE(i) == CropStage::MATURE) {
                s.DTSUM(i) = 0.0; // 成熟期停止发育
                s.DVR(i) = 0.0;
            }
            //报错机制Problem: no stage defined
            else{
                printf("Error: no stage defined in phenology calculation.\n");
            }

        });
    }

    // 更新作物状态（积分发育过程）
    void integrate(CropDynamicState &s, const CropState &p, real delt = 1.0) {
        Kokkos::parallel_for("Pheno_Integrate", s.nCells, KOKKOS_LAMBDA(const int i) {
            // 1. 【春化积分】仅当IDSL>=2且处于VEGETATIVE阶段
            
            if (p.p.IDSL >= 2) {
                if (s.STAGE(i) == CropStage::VEGETATIVE) {
                    s.VERN(i) += s.VERNR(i) * delt; // 累积春化量
                    // 【关键】春化完成条件
                    if (s.VERN(i) >= p.p.VERNSAT) { // 达到饱和阈值
                        s.ISVERNALISED(i) = 1; // 标记春化完成
                    } else if (s.DVS(i) >= 0.3 && s.ISVERNALISED(i) == 0) {
                        // 【备用逻辑】若DVS达到临界值（0.3）且未完成春化，强制完成
                        s.ISVERNALISED(i) = 1;
                    }
                }
            }

            // 2. 【发育阶段积分】
            if (s.STAGE(i) == CropStage::EMERGING) {
                s.DVS(i) += s.DVR(i) * delt; // 更新DVS（从-0.1→0.0）
                if (s.DVS(i) >= 0.0) { // 达到出苗阈值
                    s.STAGE(i) = CropStage::VEGETATIVE; // 转入营养生长期
                    s.DVS(i) = 0.0; // 重置DVS（从0.0开始）
                    s.TSUM(i) = 0.0; // 重置累积积温
                }
            }
            else if (s.STAGE(i) == CropStage::VEGETATIVE) {
                s.TSUM(i) += s.DTSUM(i) * delt; // 累积营养生长期积温
                s.DVS(i) += s.DVR(i) * delt;    // 更新DVS（从0.0→1.0）
                if (s.DVS(i) >= 1.0) { // 达到营养生长期终点
                    s.STAGE(i) = CropStage::REPRODUCTIVE; // 转入生殖生长期
                    // 【注】标准WOFOST中DVS不会重置为1.0，此处按Python代码逻辑强制赋值
                    s.DVS(i) = 1.0; 
                }
            }
            else if (s.STAGE(i) == CropStage::REPRODUCTIVE) {
                s.TSUM(i) += s.DTSUM(i) * delt; // 累积生殖生长期积温
                s.DVS(i) += s.DVR(i) * delt;    // 更新DVS（从1.0→DVSEND）
                if (s.DVS(i) >= p.p.DVSEND) { // 达到成熟阈值（通常=2.0）
                    s.STAGE(i) = CropStage::MATURE; // 转入成熟期
                    s.DVS(i) = p.p.DVSEND; // 重置DVS为终点值
                }
            }
            else if (s.STAGE(i) == CropStage::MATURE) {
                // 成熟期不再变化

            }

            else{
                printf("Error: no stage defined in phenology integration.\n");
            }
        });
    }


    // --- 移植的 daylength 函数 ---
    // 计算日长
    // day_of_year: 年积日 (1-365)
    // latitude: 纬度 (degrees)
    // angle: 太阳角度阈值 (degrees)，默认为 -4.0
    KOKKOS_INLINE_FUNCTION
    real calculate_daylength(int day_of_year, real latitude, real angle = -4.0) const {
        const real PI = 3.14159265358979323846;
        const real RAD = 0.017453292519943295; // PI / 180.0

        // 纬度范围检查 (在高性能代码中通常为了速度省略报错，或者仅在Debug模式检查)
        // if (std::abs(latitude) > 90.0) { ... }

        // calculate daylength logic derived from WOFOST ASTRO.FOR
        // DEC calculation
        real sin_obliquity = sin(23.45 * RAD); 
        real cos_pos = cos(2.0 * PI * (real(day_of_year) + 10.0) / 365.0);
        real DEC = -asin(sin_obliquity * cos_pos);

        real SINLD = sin(RAD * latitude) * sin(DEC);
        real COSLD = cos(RAD * latitude) * cos(DEC);
        
        // AOB calculation
        // AOB = (-sin(ANGLE*RAD)+SINLD)/COSLD
        real AOB = (-sin(angle * RAD) + SINLD) / COSLD;

        real DAYLP = 0.0;
        if (std::abs(AOB) <= 1.0) {
            DAYLP = 12.0 * (1.0 + 2.0 * asin(AOB) / PI);
        } else if (AOB > 1.0) {
            DAYLP = 24.0;
        } else {
            DAYLP = 0.0;
        }

        return DAYLP;
    }



};

#endif