/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _WOFOST_INIT_H_
#define _WOFOST_INIT_H_
#include "../define.h"
#include "Wofost72.h"
#include "CropState.h"
#include "MeteoState.h"
#include "CropInit.h"
#include "MeteoInit.h"
#include "../GwState.h"
#include "../GwDomain.h"
#include "../Parallel.h"
#include <string>
#include <iostream>
#include <ctime>

class WofostInit
{

public:
    // Main initialization function
    // Returns 1 on success, 0 on failure
    int initialize_wofost(Wofost72 &wofost, CropState &cropParam, MeteoState &meteo,
                          CropInit &cropInit, MeteoInit &meteoInit,
                          GwState &gw, GwDomain &gdom, Domain &dom, // dom needed for start time
                          Parallel &par, FileIO &io, std::string inFolder, std::string outFolder)
    {

        if (par.masterproc)
            std::cout << GOK "Initializing WOFOST Crop Model..." << std::endl;

        // 1. Read Crop Parameters
        std::string cropFile = inFolder + "crop.input";
        if (!cropInit.readCropParameters(cropFile, cropParam, par))
        {
            if (par.masterproc)
                std::cerr << RERROR "Failed to read " << cropFile << std::endl;
            return 0;
        }

        // 2. Read Meteo Data
        std::string meteoFile = inFolder + "meteo.input";
        if (!meteoInit.readMeteoFile(meteoFile, meteo, par))
        {
            if (par.masterproc)
                std::cerr << RERROR "Failed to read " << meteoFile << std::endl;
            return 0;
        }

        // 3. Allocate Wofost Memory based on subsurface grid (with Halo)
        // int n_crop_cells = gdom.nCellSw;
        int n_crop_cells = dom.nCell;

        // Ensure Wofost memory is allocated
        wofost.allocate(n_crop_cells);

        // 4. Calculate Initial Root Zone Soil Moisture
        // We need a temporary view to hold initial SM for initialization
        realArr init_root_zone_SM("InitRootSM", n_crop_cells);

        // Assuming we can compute this from gw.wc (volumetric water content)
        // We need to know the initial rooting depth (RDI) to average over.
        real RDI = cropParam.p.RDI;

        // 计算初始根区含水量方法一
        // wofost初始根区含水率由serghei中gw读取
        // 但是这里没有根除初始根长计算根区含水量，而是简化为取地表层（第一层）含水量作为初始根区含水量
        // Kokkos::parallel_for("WofostInit_CalcSM", n_crop_cells, KOKKOS_LAMBDA(int iGlob) {
        //     int ii, jj;
        //     // nxhc and nyhc are dimensions with halo
        //     unpackIndicesUniformGrid(iGlob, gdom.ny + 2 * hc, gdom.nx + 2 * hc, jj, ii);
        //     real total_water = 0.0;
        //     real total_depth = 0.0;
        //     int k_top = 0; // Assuming k=0 is top physical layer
        //     int iGlobGWTop = (hc+k_top)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
        //     real top_wc = gw.wc(iGlobGWTop, 0); // Use index 0 for init
        //     init_root_zone_SM(iGlob) = top_wc;
        // });

        // 计算初始含水量方法二：根据初始根长EDI计算根区含水量
        //  Cache halo cells to avoid macro expansion issues

        const int nx = gdom.nx;
        const int nxhc = gdom.nxhc;
        const int nyhc = gdom.nyhc;
        const int nz = gdom.nz;

        // Convert RDI from cm to meters
        real rd_m = RDI * 0.01;

        Kokkos::parallel_for("CalcRootSM", gdom.nCellSw, KOKKOS_LAMBDA(const int i_surf) {
            if (rd_m <= 1e-6)
            {
                // 如果根深极小，取表层含水率
                // 映射到地表单元对应的顶层地下网格单元 (k=0, top layer)
                int j = i_surf / nx;
                int i = i_surf % nx;
                int iGlob_top = (hc + 0) * nxhc * nyhc +
                                (hc + j) * nxhc +
                                (hc + i);
                init_root_zone_SM(i_surf) = gw.wc(iGlob_top, 0); // Use init step (index 0)
                return;
            }

            real total_water_depth = 0.0; // [m]
            real total_eff_depth = 0.0;   // [m]

            // 2. 遍历该地表单元对应的垂直土柱
            int i, j;

            // Manual unpack for uniform grid (row-major)
            j = i_surf / nx;
            i = i_surf % nx;

            // In SERGHEI, k=0 is TOP LAYER (smallest depth), k=nz-1 is BOTTOM (largest depth)
            // So we iterate from top (k=0) to bottom (k=nz-1)

            for (int k = 0; k < nz; k++)
            {
                int iGlob = (hc + k) * nxhc * nyhc +
                            (hc + j) * nxhc +
                            (hc + i);

                real dz = gdom.dz(iGlob);              // Layer thickness [m]
                real depth_center = gdom.depth(iGlob); // Depth of cell center [m]
                real depth_top = depth_center - 0.5 * dz;
                real depth_bot = depth_center + 0.5 * dz;

                // Intersection of layer [depth_top, depth_bot] with root zone [0, rd_m]
                real overlap_top = fmax(0.0, depth_top);
                real overlap_bot = fmin(rd_m, depth_bot);

                real eff_dz = overlap_bot - overlap_top;

                if (eff_dz <= 0.0)
                {
                    // Layer is completely below root zone
                    break;
                }

                // Accumulate water
                // gw.wc is volumetric water content [-]
                total_water_depth += gw.wc(iGlob, 0) * eff_dz; // Use init step (index 0)
                total_eff_depth += eff_dz;
            }

            if (total_eff_depth > 0.0)
            {
                init_root_zone_SM(i_surf) = total_water_depth / total_eff_depth;
            }
            else
            {
                init_root_zone_SM(i_surf) = 0.0; // Should not happen if rd_m > 0
            }
#if DEBUG_CROP_GROWTH_MODEL
            printf("i_surf: %d, rd: %f, root_zone_SM: %f\n", i_surf, rd_m, init_root_zone_SM(i_surf));

#endif
        });

        // 5. Calculate start DOY (Day of Year) - REPLACED LOGIC
        // Use the start DOY detected from the Meteo file
        int start_doy = meteo.start_doy;

        if (par.masterproc)
        {
            printf("Simulation Start Configuration:\n");
            printf("  > Start Year: %d\n", meteo.start_year);
            printf("  > Start DOY : %d\n", start_doy);
        }

        // 6. Initialize WOFOST Internal States
        wofost.initialize(cropParam, meteo, init_root_zone_SM, start_doy);

        // 初始化输出结果文件
        io.outputIniCrop(wofost, dom, par, outFolder);

        // Debug: Check initial values
        if (par.masterproc)
        {
            printf("[INIT-DEBUG] After WOFOST initialize:\n");
            printf("  s.LAI(0) = %.4f\n", wofost.s.LAI(0));
            printf("  s.WLV(0) = %.4f\n", wofost.s.WLV(0));
            printf("  s.DVS(0) = %.4f\n", wofost.s.DVS(0));
            printf("  lds.LAI(0) = %.4f\n", wofost.lds.LAI(0));
            printf("  lds.WLV(0) = %.4f\n", wofost.lds.WLV(0));
            printf("  lds.LASUM(0) = %.4f\n", wofost.lds.LASUM(0));
        }

        if (par.masterproc)
            std::cout << GOK "WOFOST Initialized." << std::endl;

        return 1;
    }
};

#endif
