/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _SOURCESINK_CROP_H_
#define _SOURCESINK_CROP_H_

#include "define.h"
#include "GwDomain.h"
#include "GwState.h"
#include "SourceSink.h"
#include "cropsrc/CropDynamicState.h" // Needed to access RD and LAI

// 扩展的 SourceSink 处理类，专门用于处理与作物模型耦合时的源汇项
class SourceSinkCrop {

public:
    // 计算每个地表单元对应的根区平均体积含水率
    // 输入：
    //   gw: 地下水状态 (包含 wc - 含水率)
    //   gdom: 地下水域 (包含网格几何信息 dz, depth)
    //   RD: 作物根深数组 (来自 WOFOST, dimension: nCellSwMem)
    // 输出：
    //   root_zone_SM: 根区平均含水率数组 (dimension: nCellSwMem)
    static void compute_root_zone_moisture(const GwState& gw, const GwDomain& gdom, 
                                           const realArr& RD, realArr& root_zone_SM) {
        printf("Computing root zone moisture...\n");
        
        // 遍历所有地表单元 (Surface Cells)
        Kokkos::parallel_for("CalcRootSM", gdom.nCellSw, KOKKOS_LAMBDA(const int i_surf) {
            
            // 1. 获取当前单元的根深 (cm -> m)
            real rd = RD(i_surf) * 0.01; 

            // 处理根深为0或极小的情况
            if (rd <= 1e-6) {
                // 映射到地表单元对应的顶层地下网格单元 (k=0, top layer)
                int j = i_surf / gdom.nx;
                int i = i_surf % gdom.nx;
                int iGlob_top = (hc + 0) * gdom.nxhc * gdom.nyhc + 
                                (hc + j) * gdom.nxhc + 
                                (hc + i);
                
                // 边界检查
                if (iGlob_top >= 0 && iGlob_top < gw.wc.extent(0)) {
                    root_zone_SM(i_surf) = gw.wc(iGlob_top, 1); // Use current step (index 1)
                } else {
                    root_zone_SM(i_surf) = 0.0; // Fallback
                }
                return;
            }

            real total_water_depth = 0.0; // [m]
            real total_eff_depth = 0.0;   // [m]
            
            // 2. 遍历该地表单元对应的垂直土柱
            int j = i_surf / gdom.nx;
            int i = i_surf % gdom.nx;
            
            // 从顶层 (k=0) 开始向下遍历
            for (int k = 0; k < gdom.nz; k++) {
                int iGlob = (hc + k) * gdom.nxhc * gdom.nyhc + 
                            (hc + j) * gdom.nxhc + 
                            (hc + i);
                
                // 边界检查
                if (iGlob < 0 || iGlob >= gw.wc.extent(0)) {
                    break;
                }
                
                real dz = gdom.dz(iGlob);     // Layer thickness [m]
                real depth_center = gdom.depth(iGlob); // Depth of cell center [m]
                real depth_top = depth_center - 0.5 * dz;
                real depth_bot = depth_center + 0.5 * dz;
                
                // 计算根区与当前层的重叠部分
                real overlap_top = fmax(0.0, depth_top);
                real overlap_bot = fmin(rd, depth_bot);
                
                real eff_dz = overlap_bot - overlap_top;
                
                if (eff_dz <= 0.0) {
                    // 当前层完全在根区以下，由于深度随k增加，后续层也都在根区以下
                    break; 
                }
                
                // 累加水量
                real wc_val = gw.wc(iGlob, 1); // Use current step (index 1)
                total_water_depth += wc_val * eff_dz; 
                total_eff_depth += eff_dz;
                
                // 防止无限循环的安全检查
                if (total_eff_depth > rd * 10.0) {
                    break; // 已经累积了过多的深度，防止死循环
                }
            }
            
            if (total_eff_depth > 0.0) {
                root_zone_SM(i_surf) = total_water_depth / total_eff_depth;
            } else {
                // 如果没有重叠，使用表层含水率
                int iGlob_top = (hc + 0) * gdom.nxhc * gdom.nyhc + 
                                (hc + j) * gdom.nxhc + 
                                (hc + i);
                if (iGlob_top >= 0 && iGlob_top < gw.wc.extent(0)) {
                    root_zone_SM(i_surf) = gw.wc(iGlob_top, 1);
                } else {
                    root_zone_SM(i_surf) = 0.0;
                }
            }

        });

    }

    // 更新源汇项参数（基于 WOFOST 输出）
    // 输入：
    //   gdom: 地下水域
    //   ss: 源汇项对象 (将被修改)
    //   LAI: 叶面积指数数组 (来自 WOFOST)
    //   RD: 根深数组 (来自 WOFOST)
    //   current_time: 当前模拟时间
    static void update_crop_parameters(const GwDomain& gdom, SourceSink& ss, 
                                       const realArr& LAI, const realArr& RD, 
                                       real current_time) {
        
        // 遍历所有定义的地下水源汇项
        // 通常蒸散发 (ET) 是作为一个特定的 SourceSink (type 0) 定义的
        for (auto& gwss_item : ss.gwss) {
            if (gwss_item.sstype == 0) { // Type 0 is ET
                
                // 设置实时数据标志
                gwss_item.use_realtime_data = true;
                
                // 分配实时数据数组（如果尚未分配）
                if (gwss_item.lai_realtime.data() == nullptr) {
                    gwss_item.lai_realtime = realArr("lai_realtime", gdom.nCellSwMem);
                }
                if (gwss_item.rd_realtime.data() == nullptr) {
                    gwss_item.rd_realtime = realArr("rd_realtime", gdom.nCellSwMem);
                }
                
                // 复制 LAI 和 RD 数据到 GwSS 对象
                // LAI 和 RD 来自 WOFOST，大小为 nCellSwMem（表面网格，包含halo）
                Kokkos::deep_copy(gwss_item.lai_realtime, LAI);
                Kokkos::deep_copy(gwss_item.rd_realtime, RD);
            }
        }
    }
    
    // 这是一个辅助函数，用于在 Host 端准备数据，或者作为上述 update_crop_parameters 的替代
    // 它计算 T_split = (1 - exp(-k*LAI)) 和 E_split = exp(-k*LAI)
    // 供后续 Kernel 使用
    static void compute_evap_partitioning(const GwDomain& gdom, 
                                          const realArr& LAI, 
                                          realArr& T_split, realArr& E_split,
                                          real k_ext = 0.5) {
        
        Kokkos::parallel_for("CalcETSplit", gdom.nCellSwMem, KOKKOS_LAMBDA(const int i_surf) {
            real lai = LAI(i_surf);
            real exp_term = exp(-k_ext * lai);
            
            T_split(i_surf) = 1.0 - exp_term;
            E_split(i_surf) = exp_term;
        });
    }
};

#endif
