/**
 * @file WofostInit.h
 * @brief Top-level initialization routine for the WOFOST 7.2 crop growth model.
 *
 * Orchestrates the complete initialization sequence:
 *   1. Read crop parameters from the input file (via CropInit).
 *   2. Read meteorological driving data from the input file (via MeteoInit).
 *   3. Allocate memory for all WOFOST dynamic state arrays.
 *   4. Compute the initial root-zone soil moisture by integrating the
 *      volumetric water content from the groundwater model over the
 *      initial rooting depth (RDI).
 *   5. Determine the simulation start day-of-year (DOY) from meteo data.
 *   6. Initialize all internal WOFOST states (phenology, partitioning,
 *      organ biomass, LAI, etc.).
 *   7. Initialize the crop output file.
 *
 * This class bridges the SERGHEI hydrological framework (Domain, GwState,
 * GwDomain) with the standalone WOFOST crop model (Wofost72).
 *
 * @see Wofost72.h    for the main model class being initialized.
 * @see CropInit.h    for the crop parameter file reader.
 * @see MeteoInit.h   for the meteorological data file reader.
 */

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
    /**
     * @brief Main initialization entry point for the WOFOST crop model.
     *
     * Reads all input data, allocates memory, computes initial conditions
     * (particularly root-zone soil moisture from the groundwater state),
     * and initializes all WOFOST sub-modules.
     *
     * @param wofost       WOFOST model container to initialize.
     * @param cropParam    CropState to populate with static parameters.
     * @param meteo        MeteoState to populate with weather data.
     * @param cropInit     Crop parameter reader instance.
     * @param meteoInit    Meteorological data reader instance.
     * @param gw           Groundwater state (provides volumetric water content).
     * @param gdom         Groundwater domain grid (provides cell geometry).
     * @param dom          Surface domain (provides nCell and start time).
     * @param par          Parallel controller (for masterproc output).
     * @param io           File I/O controller (for output file initialization).
     * @param inFolder     Input directory path (trailing slash expected).
     * @param outFolder    Output directory path (trailing slash expected).
     * @return 1 on success, 0 on failure.
     */
    int initialize_wofost(Wofost72 &wofost, CropState &cropParam, MeteoState &meteo,
                          CropInit &cropInit, MeteoInit &meteoInit,
                          GwState &gw, GwDomain &gdom, Domain &dom,
                          Parallel &par, FileIO &io, std::string inFolder, std::string outFolder)
    {

        if (par.masterproc)
            std::cout << GOK "Initializing WOFOST Crop Model..." << std::endl;

        /* ================================================================ */
        /* Step 1: Read crop parameters from input file                     */
        /* ================================================================ */
        std::string cropFile = inFolder + "cropparameter.input";
        if (!cropInit.readCropParameters(cropFile, cropParam, par))
        {
            if (par.masterproc)
                std::cerr << RERROR "Failed to read " << cropFile << std::endl;
            return 0;
        }

        /* ================================================================ */
        /* Step 2: Read meteorological driving data                         */
        /* ================================================================ */
        std::string meteoFile = inFolder + "cropmeteo.input";
        if (!meteoInit.readMeteoFile(meteoFile, meteo, par))
        {
            if (par.masterproc)
                std::cerr << RERROR "Failed to read " << meteoFile << std::endl;
            return 0;
        }

        /* ================================================================ */
        /* Step 3: Allocate all WOFOST dynamic state arrays                 */
        /*                                                                  */
        /* Uses the number of surface grid cells (without halo) so that     */
        /* each cell has its own crop simulation state.                     */
        /* ================================================================ */
        int n_crop_cells = dom.nCell;
        wofost.allocate(n_crop_cells);

        /* ================================================================ */
        /* Step 4: Compute initial root-zone soil moisture                  */
        /*                                                                  */
        /* The initial root-zone soil moisture is calculated by integrating */
        /* the volumetric water content (gw.wc) over the vertical soil      */
        /* column from the surface down to the initial rooting depth (RDI). */
        /*                                                                  */
        /* For each surface cell:                                           */
        /*   - RDI is converted from cm to metres.                          */
        /*   - The vertical soil column is traversed top to bottom.         */
        /*   - Each layer's contribution is weighted by its overlap with    */
        /*     the root zone [0, RDI_metres].                               */
        /*   - The result is a depth-averaged volumetric water content.     */
        /*                                                                  */
        /* Edge case: if RDI is effectively zero, the top-layer water       */
        /* content is used directly.                                        */
        /* ================================================================ */
        realArr init_root_zone_SM("InitRootSM", n_crop_cells);

        /* Initial rooting depth from crop parameters */
        real RDI = cropParam.p.RDI;

        /* Cache halo and grid dimensions for index computation */
        const int nx   = gdom.nx;
        const int nxhc = gdom.nxhc;
        const int nyhc = gdom.nyhc;
        const int nz   = gdom.nz;

        /* Convert RDI from [cm] to [m] */
        real rd_m = RDI * 0.01;

        Kokkos::parallel_for("CalcRootSM", n_crop_cells, KOKKOS_LAMBDA(const int i_surf) {

            /* Edge case: negligible rooting depth -> use top-layer water content */
            if (rd_m <= 1e-6)
            {
                /*
                 * Map the surface cell index to the corresponding top-layer
                 * groundwater cell (k=0 is the top layer in SERGHEI).
                 */
                int j = i_surf / nx;
                int i = i_surf % nx;
                int iGlob_top = (gdom.hc + 0) * nxhc * nyhc +
                                (gdom.hc + j) * nxhc +
                                (gdom.hc + i);
                init_root_zone_SM(i_surf) = gw.wc(iGlob_top, 0);
                return;
            }

            /* Accumulate depth-weighted water content over the root zone */
            real total_water_depth = 0.0;   /* [m] depth-integrated water */
            real total_eff_depth  = 0.0;    /* [m] total effective depth   */

            /* Unpack surface cell index to (i, j) on the uniform grid */
            int i, j;
            j = i_surf / nx;
            i = i_surf % nx;

            /*
             * Iterate over the vertical soil column.
             * In SERGHEI: k=0 is the TOP layer, k=nz-1 is the BOTTOM layer.
             * Each layer's overlap with the root zone [0, rd_m] is computed
             * and the water content is weighted accordingly.
             */
            for (int k = 0; k < nz; k++)
            {
                int iGlob = (gdom.hc + k) * nxhc * nyhc +
                            (gdom.hc + j) * nxhc +
                            (gdom.hc + i);

                real dz           = gdom.dz(iGlob);           /* Layer thickness [m]      */
                real depth_center = gdom.depth(iGlob);        /* Depth of cell centre [m] */
                real depth_top    = depth_center - 0.5 * dz;  /* Top of layer [m]         */
                real depth_bot    = depth_center + 0.5 * dz;  /* Bottom of layer [m]      */

                /* Compute intersection of this layer with the root zone [0, rd_m] */
                real overlap_top = fmax(0.0, depth_top);
                real overlap_bot = fmin(rd_m, depth_bot);
                real eff_dz = overlap_bot - overlap_top;

                if (eff_dz <= 0.0)
                {
                    /* This layer is entirely below the root zone; stop iterating */
                    break;
                }

                /*
                 * Accumulate water: gw.wc is volumetric water content [m3/m3].
                 * Multiplying by effective thickness gives water depth [m].
                 */
                total_water_depth += gw.wc(iGlob, 0) * eff_dz;
                total_eff_depth  += eff_dz;
            }

            /* Depth-averaged volumetric water content in the root zone */
            if (total_eff_depth > 0.0)
            {
                init_root_zone_SM(i_surf) = total_water_depth / total_eff_depth;
            }
            else
            {
                init_root_zone_SM(i_surf) = 0.0;  /* Should not occur if rd_m > 0 */
            }
        });

        /* ================================================================ */
        /* Step 5: Determine simulation start day-of-year (DOY)             */
        /*                                                                  */
        /* The start DOY is extracted from the meteorological input file,   */
        /* which records the date of the first weather record.              */
        /* ================================================================ */
        int start_doy = meteo.start_doy;

        if (par.masterproc)
        {
            printf("Simulation Start Configuration:\n");
            printf("  > Start Year: %d\n", meteo.start_year);
            printf("  > Start DOY : %d\n", start_doy);
        }

        /* ================================================================ */
        /* Step 6: Initialize all WOFOST internal states                   */
        /*                                                                  */
        /* This sets phenology (DVS, TSUM, STAGE), partitioning fractions, */
        /* organ biomass from TDWI, LAI from LAIEM, rooting depth from RDI,*/
        /* and zeroes all cumulative accumulators.                          */
        /* ================================================================ */
        wofost.initialize(cropParam, meteo, init_root_zone_SM, start_doy);

        /* ================================================================ */
        /* Step 7: Initialize the crop output results file                  */
        /* ================================================================ */
        io.outputIniCrop(wofost, dom, par, outFolder);

        /* --- Debug: verify initial values on master process --- */
        if (par.masterproc)
        {
            printf("[INIT-DEBUG] After WOFOST initialize:\n");
            printf("  s.LAI(0) = %.4f\n", wofost.s.LAI(0));
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

#endif /* _WOFOST_INIT_H_ */
