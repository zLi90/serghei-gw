/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _CROP_STATE_H_
#define _CROP_STATE_H_

// --- 恢复 SERGHEI 原生头文件 ---
#include "../const.h"
#include "../define.h"
#include "../SArray.h"
#include <map>
#include <string>
#include <iostream>

struct CropParams {
    // 标量参数
    real TBASEM, TEFFMX, TSUMEM;
    int IDSL;
    real DLO, DLC, TSUM1, TSUM2, DVSI, DVSEND;
    real TDWI, LAIEM, RGRLAI;
    real SPAN, TBASE, SPA;
    real CVL, CVO, CVR, CVS;
    real Q10, RML, RMO, RMR, RMS;
    real PERDL;
    real CFET, DEPNR;
    int IAIRDU, IOX;
    real RDI, RRI, RDMCR;

    // Vernalization parameters
    real VERNBASE, VERNSAT;
    
    // Biomass parameters
    real NMINSO, NMINVE, NMAXSO, NMAXVE;
    real PMINSO, PMINVE, PMAXSO, PMAXVE;
    real KMINSO, KMINVE, KMAXSO, KMAXVE;
    real YZERO, NFIX;
    // ... 其他参数如果需要 ...
};

class CropState {
public:
    CropParams p; 
    std::string cropName;
    
    // 表格参数 (Kokkos View)
    std::map<std::string, realArr> tables;

    // 分配表格内存
    // 修改 int size 为 size_t size，匹配 Kokkos 的维度定义
// 分配表格内存
    // 1. 参数类型改为 size_t
    void allocateTable(std::string name, size_t size) {
        // 2. 直接传递 std::string 类型的 name，Kokkos 原生支持
        tables[name] = realArr(name, size);
    }

    // 检查表格是否存在
    bool hasTable(std::string name) const {
        return tables.find(name) != tables.end();
    }

    // 获取表格 (Kokkos View)
    realArr getTable(std::string name) {
        if (tables.find(name) != tables.end()) {
            return tables[name];
        } else {
            // 在 SERGHEI 中通常使用特定的错误处理，这里简单打印
            // 注意：在 Device 代码中调用此函数可能会有问题（std::map 不支持 Device），
            // 建议在初始化阶段将需要的 View 提取出来传递给 Kernel。
            std::cerr << "Error: Table " << name << " not found." << std::endl;
            return realArr(); // 返回空 View
        }
    }
};
#endif
