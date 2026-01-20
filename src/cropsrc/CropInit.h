/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _CROP_INIT_H_
#define _CROP_INIT_H_


#include "../define.h"
#include "../Parallel.h"
#include "CropState.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <cstring>

class CropInit {
    // 内部类：解析行
    class PsLn {
    public:
        std::string line;
        std::string key;
        std::stringstream value;
        void lowercase() { std::for_each(line.begin(), line.end(), [](char &c) { c = ::tolower(c); }); }
        void parse() {
            key.clear();
            value.clear();
            if (!line.empty() && line.find("//", 0) != 0) {
                size_t splitloc = line.find(':', 0);
                if (splitloc != std::string::npos) {
                    key = line.substr(0, splitloc);
                    // 移除 Key 中的空格和制表符
                    key.erase(std::remove(key.begin(), key.end(), ' '), key.end());
                    key.erase(std::remove(key.begin(), key.end(), '\t'), key.end());
                    
                    std::string val = line.substr(splitloc + 1, line.length() - splitloc);
                    // 移除注释
                    size_t splitter = val.find("//", 0);
                    std::string strloc;
                    if (splitter != std::string::npos) { strloc = val.substr(0, splitter); }
                    else { strloc = val; }
                    value.clear();
                    value.str(strloc);
                }
            }
        }
    };

public:
    int readCropParameters(std::string fNameIn, CropState &cs, Parallel &par) {
        std::ifstream fInStream(fNameIn);
        std::string line;
        PsLn pline;
        
        if (!fInStream.is_open()) {
            if (par.masterproc) std::cerr << RERROR "Unable to open crop parameter file: " << fNameIn << "\n";
            return 0;
        }

        while (std::getline(fInStream, line)) {
            pline.line = line;
            // pline.lowercase(); // WOFOST keys are often uppercase, keep original case or enforce convention
            pline.parse();

            if (!pline.key.empty()) {
                if (!strcmp("CRPNAM", pline.key.c_str())) { pline.value >> cs.cropName; }
                // --- Scalars ---
                else if (!strcmp("TBASEM", pline.key.c_str())) { pline.value >> cs.p.TBASEM; }
                else if (!strcmp("TEFFMX", pline.key.c_str())) { pline.value >> cs.p.TEFFMX; }
                else if (!strcmp("TSUMEM", pline.key.c_str())) { pline.value >> cs.p.TSUMEM; }
                else if (!strcmp("IDSL", pline.key.c_str()))   { pline.value >> cs.p.IDSL; }
                else if (!strcmp("DLO", pline.key.c_str()))    { pline.value >> cs.p.DLO; }
                else if (!strcmp("DLC", pline.key.c_str()))    { pline.value >> cs.p.DLC; }
                else if (!strcmp("TSUM1", pline.key.c_str()))  { pline.value >> cs.p.TSUM1; }
                else if (!strcmp("TSUM2", pline.key.c_str()))  { pline.value >> cs.p.TSUM2; }
                else if (!strcmp("DVSI", pline.key.c_str()))   { pline.value >> cs.p.DVSI; }
                else if (!strcmp("DVSEND", pline.key.c_str())) { pline.value >> cs.p.DVSEND; }
                else if (!strcmp("TDWI", pline.key.c_str()))   { pline.value >> cs.p.TDWI; }
                else if (!strcmp("LAIEM", pline.key.c_str()))  { pline.value >> cs.p.LAIEM; }
                else if (!strcmp("RGRLAI", pline.key.c_str())) { pline.value >> cs.p.RGRLAI; }
                else if (!strcmp("SPAN", pline.key.c_str()))   { pline.value >> cs.p.SPAN; }
                else if (!strcmp("TBASE", pline.key.c_str()))  { pline.value >> cs.p.TBASE; }
                else if (!strcmp("SPA", pline.key.c_str()))    { pline.value >> cs.p.SPA; }
                else if (!strcmp("CVL", pline.key.c_str()))    { pline.value >> cs.p.CVL; }
                else if (!strcmp("CVO", pline.key.c_str()))    { pline.value >> cs.p.CVO; }
                else if (!strcmp("CVR", pline.key.c_str()))    { pline.value >> cs.p.CVR; }
                else if (!strcmp("CVS", pline.key.c_str()))    { pline.value >> cs.p.CVS; }
                else if (!strcmp("Q10", pline.key.c_str()))    { pline.value >> cs.p.Q10; }
                else if (!strcmp("RML", pline.key.c_str()))    { pline.value >> cs.p.RML; }
                else if (!strcmp("RMO", pline.key.c_str()))    { pline.value >> cs.p.RMO; }
                else if (!strcmp("RMR", pline.key.c_str()))    { pline.value >> cs.p.RMR; }
                else if (!strcmp("RMS", pline.key.c_str()))    { pline.value >> cs.p.RMS; }
                else if (!strcmp("PERDL", pline.key.c_str()))  { pline.value >> cs.p.PERDL; }
                else if (!strcmp("CFET", pline.key.c_str()))   { pline.value >> cs.p.CFET; }
                else if (!strcmp("DEPNR", pline.key.c_str()))  { pline.value >> cs.p.DEPNR; }
                else if (!strcmp("IAIRDU", pline.key.c_str())) { pline.value >> cs.p.IAIRDU; }
                else if (!strcmp("IOX", pline.key.c_str()))    { pline.value >> cs.p.IOX; }
                else if (!strcmp("RDI", pline.key.c_str()))    { pline.value >> cs.p.RDI; }
                else if (!strcmp("RRI", pline.key.c_str()))    { pline.value >> cs.p.RRI; }
                else if (!strcmp("RDMCR", pline.key.c_str()))  { pline.value >> cs.p.RDMCR; }
                // Nutrients
                else if (!strcmp("NMINSO", pline.key.c_str())) { pline.value >> cs.p.NMINSO; }
                else if (!strcmp("NMINVE", pline.key.c_str())) { pline.value >> cs.p.NMINVE; }
                else if (!strcmp("NMAXSO", pline.key.c_str())) { pline.value >> cs.p.NMAXSO; }
                else if (!strcmp("NMAXVE", pline.key.c_str())) { pline.value >> cs.p.NMAXVE; }
                else if (!strcmp("PMINSO", pline.key.c_str())) { pline.value >> cs.p.PMINSO; }
                else if (!strcmp("PMINVE", pline.key.c_str())) { pline.value >> cs.p.PMINVE; }
                else if (!strcmp("PMAXSO", pline.key.c_str())) { pline.value >> cs.p.PMAXSO; }
                else if (!strcmp("PMAXVE", pline.key.c_str())) { pline.value >> cs.p.PMAXVE; }
                else if (!strcmp("KMINSO", pline.key.c_str())) { pline.value >> cs.p.KMINSO; }
                else if (!strcmp("KMINVE", pline.key.c_str())) { pline.value >> cs.p.KMINVE; }
                else if (!strcmp("KMAXSO", pline.key.c_str())) { pline.value >> cs.p.KMAXSO; }
                else if (!strcmp("KMAXVE", pline.key.c_str())) { pline.value >> cs.p.KMAXVE; }
                else if (!strcmp("YZERO", pline.key.c_str()))  { pline.value >> cs.p.YZERO; }
                else if (!strcmp("NFIX", pline.key.c_str()))   { pline.value >> cs.p.NFIX; }
                
                // --- Tables (Arrays) ---
                else if (isTableKey(pline.key)) {
                    std::vector<real> values;
                    real val;
                    // 读取剩余的值到 vector
                    while (pline.value >> val) {
                        values.push_back(val);
                    }
                    
                    if (values.size() % 2 != 0) {
                        if (par.masterproc) std::cerr << YEXC "Warning: Odd number of elements in table " << pline.key << ". Check input file.\n";
                    }
                    
                    // 1. 在 CropState 中分配 View
                    cs.allocateTable(pline.key, values.size());
                    
                    // 2. 获取 View 引用
                    realArr tableRef = cs.getTable(pline.key);
                    
                    // 3. 填充数据 (使用 Kokkos 标准流程)
                    // 创建一个 HostMirror (主机端副本)
                    realArr::HostMirror h_table = Kokkos::create_mirror_view(tableRef);
                    
                    // 将 vector 数据填入 HostMirror
                    for (size_t i = 0; i < values.size(); ++i) {
                        h_table(i) = values[i];
                    }
                    
                    // 将数据从 HostMirror 深拷贝到 Device View
                    Kokkos::deep_copy(tableRef, h_table);
                }
            }
        }

        // 在 WofostInit.h 中
            if (par.masterproc) {
                printf("Crop Params Check: TDWI=%.4f, LAIEM=%.4f\n", cs.p.TDWI, cs.p.LAIEM);
            }
        
        if (par.masterproc) std::cout << GOK "Crop parameters read from " << fNameIn << "\n";
        return 1;
    }

private:
    bool isTableKey(const std::string& key) {
        static const std::vector<std::string> tableKeys = {
            "DTSMTB", "SLATB", "SSATB", "KDIFTB", "EFFTB", "AMAXTB", 
            "TMPFTB", "TMNFTB", "RFSETB", "FRTB", "FLTB", "FSTB", 
            "FOTB", "RDRRTB", "RDRSTB"
        };
        return std::find(tableKeys.begin(), tableKeys.end(), key) != tableKeys.end();
    }
};
#endif