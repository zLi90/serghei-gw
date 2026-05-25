/* -*- mode: c++; c-default-style: "linux" -*- */

/*******************************************************************************
 * @file    SourceSinkCrop.h
 * @brief   Crop-hydrology coupling layer: interfaces between the WOFOST crop
 *          growth model and the SERGHEI groundwater/soil moisture model.
 *
 * @details This module provides the bridge between two coupled simulation systems:
 *
 *   1. **WOFOST crop model** - simulates crop phenology, biomass accumulation,
 *      root growth, and leaf area index (LAI) as functions of weather and
 *      soil moisture stress.
 *
 *   2. **SERGHEI groundwater model** - simulates subsurface flow and soil
 *      moisture distribution on a 3D finite-difference grid.
 *
 * The coupling is two-way:
 *   - **WOFOST -> SERGHEI**: Crop parameters (LAI, root depth RD) are passed
 *     to the groundwater model to partition potential evapotranspiration into
 *     soil evaporation and transpiration, and to determine the extraction depth.
 *   - **SERGHEI -> WOFOST**: Soil moisture from the groundwater model is
 *     aggregated over the root zone and passed back to WOFOST to compute
 *     water stress factors that regulate crop growth.
 *
 * Key functions:
 *   - compute_root_zone_moisture(): Spatially averages soil moisture over
 *     the root zone for each surface cell.
 *   - update_crop_parameters(): Passes LAI and root depth to the
 *     groundwater source/sink terms.
 *   - compute_evap_partitioning(): Splits potential ET into soil evaporation
 *     and plant transpiration using an exponential extinction model.
 *
 * @see     SourceSink.h          for the groundwater source/sink term definitions
 * @see     CropDynamicState.h    for WOFOST dynamic state (LAI, RD)
 ******************************************************************************/

#ifndef _SOURCESINK_CROP_H_
#define _SOURCESINK_CROP_H_

#include "define.h"
#include "GwDomain.h"
#include "GwState.h"
#include "SourceSink.h"
#include "cropsrc/CropDynamicState.h"  /* Needed to access RD (root depth) and LAI (leaf area index) */

/**
 * @class SourceSinkCrop
 * @brief Extended SourceSink handler for coupled crop-hydrology simulations.
 *
 * This class extends the standard SERGHEI source/sink framework to support
 * crop model coupling. It provides static methods for computing root zone
 * soil moisture from the groundwater model state, updating crop parameters
 * in the source/sink objects, and partitioning evapotranspiration.
 */
class SourceSinkCrop
{

public:
    /* ======================================================================= *
     * Root Zone Soil Moisture Computation
     * ======================================================================= */

    /**
     * @brief Compute root zone average volumetric moisture content for each
     *        surface cell.
     *
     * This function integrates the WOFOST root depth (RD) with the 3D soil
     * moisture field from the SERGHEI groundwater model. For each surface
     * cell, it traverses the vertical soil column and computes a depth-
     * weighted average of volumetric water content over the root zone.
     *
     * The algorithm:
     *   1. Convert WOFOST root depth from [cm] to [m]
     *   2. For each vertical layer in the soil column:
     *      a. Compute the overlap between the layer and the root zone [0, RD]
     *      b. Accumulate water content * overlap thickness
     *      c. Accumulate overlap thickness
     *   3. Divide total water depth by total overlap thickness to get
     *      the depth-averaged root zone moisture content
     *
     * The result is a dimensionless volumetric water content [m^3/m^3]
     * (or equivalently [cm^3/cm^3]) representing the average moisture
     * available to the crop roots in each surface cell.
     *
     * Typical volumetric moisture content ranges:
     *   Wilting point:     0.05 - 0.15 (crop-dependent)
     *   Field capacity:    0.20 - 0.40 (soil-texture-dependent)
     *   Saturation:        0.35 - 0.50 (soil-texture-dependent)
     *
     * @param[in]  gw            GwState containing the 3D soil moisture field
     *                           (gw.wc is indexed by [iGlob, time_step])
     * @param[in]  gdom          GwDomain containing grid geometry (dz, depth,
     *                           cell counts, halo padding)
     * @param[in]  RD            Root depth array from WOFOST, dimension: nCellSw [cm]
     *                           Typical range: 0 - 150 cm depending on crop and growth stage
     * @param[out] root_zone_SM  Root zone average moisture content array,
     *                           dimension: nCellSw [m^3/m^3, volumetric]
     */
    static void compute_root_zone_moisture(const GwState &gw, const GwDomain &gdom,
                                           const realArr &RD, realArr &root_zone_SM)
    {
        printf("Computing root zone moisture...\n");

        /* ----------------------------------------------------------------- *
         * Parallel loop over all surface cells.
         *
         * nCellSw is the number of surface water cells (2D footprint of
         * the groundwater domain, excluding halo cells in SERGHEI's
         * internal indexing but mapped through the halo offset 'hc').
         * ----------------------------------------------------------------- */
        Kokkos::parallel_for("CalcRootSM", gdom.nCellSw, KOKKOS_LAMBDA(const int i_surf) {

            /* ------------------------------------------------------------- *
             * Step 1: Retrieve root depth and convert from cm to m.
             *
             * Root depth (RD) comes from WOFOST and represents the maximum
             * depth to which crop roots extend below the soil surface.
             * In WOFOST, RD increases during vegetative growth and may
             * decrease slightly during senescence.
             * Units: input [cm], converted to [m] for consistency with
             * the groundwater model grid.
             * ------------------------------------------------------------- */
            real rd = RD(i_surf) * 0.01;  /* [cm] -> [m] */

            /* ------------------------------------------------------------- *
             * Edge case: root depth is zero or negligible (e.g., before crop
             * emergence or after complete senescence).
             *
             * In this case, use the volumetric moisture content of the
             * topmost soil layer (k=0) as a fallback.
             * ------------------------------------------------------------- */
            if (rd <= 1e-6)
            {
                /* Map the surface cell index to 2D (i,j) coordinates.
                 * i_surf is a linear index over the nCellSw surface cells.
                 * j = row index (y-direction), i = column index (x-direction). */
                int j = i_surf / gdom.nx;
                int i = i_surf % gdom.nx;

                /* Compute the global 3D index for the top layer (k=0).
                 * The halo offset 'hc' pads the domain on all sides to
                 * handle boundary conditions in the finite-difference scheme. */
                int iGlob_top = (gdom.hc + 0) * gdom.nxhc * gdom.nyhc +
                                (gdom.hc + j) * gdom.nxhc +
                                (gdom.hc + i);

                /* Boundary safety check before accessing the moisture array. */
                if (iGlob_top >= 0 && iGlob_top < gw.wc.extent(0))
                {
                    root_zone_SM(i_surf) = gw.wc(iGlob_top, 1); /* Current time step (index 1) */
                }
                else
                {
                    root_zone_SM(i_surf) = 0.0;  /* Fallback for out-of-bounds */
                }
                return;
            }

            /* Accumulators for depth-weighted moisture averaging. */
            real total_water_depth = 0.0;  /* Sum of (wc * effective_layer_thickness) [m] */
            real total_eff_depth = 0.0;    /* Sum of effective layer thicknesses within root zone [m] */

            /* ------------------------------------------------------------- *
             * Step 2: Decompose surface cell index into 2D (i,j) coordinates
             * for the vertical column traversal.
             * ------------------------------------------------------------- */
            int j = i_surf / gdom.nx;
            int i = i_surf % gdom.nx;

            /* ------------------------------------------------------------- *
             * Step 3: Traverse the vertical soil column from top (k=0) downward.
             *
             * For each layer, compute the overlap between the layer extent
             * and the root zone [0, rd], then accumulate the moisture
             * contribution.
             *
             * The loop breaks early when a layer is entirely below the
             * root zone (no overlap), since deeper layers will also have
             * no overlap (depth is monotonically increasing with k).
             * ------------------------------------------------------------- */
            for (int k = 0; k < gdom.nz; k++)
            {
                /* Global 3D index for layer k at surface position (i, j). */
                int iGlob = (gdom.hc + k) * gdom.nxhc * gdom.nyhc +
                            (gdom.hc + j) * gdom.nxhc +
                            (gdom.hc + i);

                /* Boundary check for the global index. */
                if (iGlob < 0 || iGlob >= gw.wc.extent(0))
                {
                    break;
                }

                /* Layer geometry:
                 * dz          = layer thickness [m]
                 * depth_center = depth of cell center below surface [m]
                 * depth_top    = depth of cell top face [m]
                 * depth_bot    = depth of cell bottom face [m] */
                real dz = gdom.dz(iGlob);              /* [m] */
                real depth_center = gdom.depth(iGlob); /* [m] */
                real depth_top = depth_center - 0.5 * dz;
                real depth_bot = depth_center + 0.5 * dz;

                /* --------------------------------------------------------- *
                 * Compute overlap between the root zone [0, rd] and the
                 * current soil layer [depth_top, depth_bot].
                 *
                 * overlap_top = max(0, depth_top)  : root zone starts at 0
                 * overlap_bot = min(rd, depth_bot) : root zone ends at rd
                 * eff_dz = overlap_bot - overlap_top : effective thickness
                 * --------------------------------------------------------- */
                real overlap_top = fmax(0.0, depth_top);
                real overlap_bot = fmin(rd, depth_bot);

                real eff_dz = overlap_bot - overlap_top;

                if (eff_dz <= 0.0)
                {
                    /* This layer is entirely below the root zone.
                     * Since depth increases with k, all subsequent layers
                     * are also below the root zone -- terminate the loop. */
                    break;
                }

                /* Accumulate the depth-weighted moisture contribution.
                 * wc_val is volumetric water content [m^3/m^3].
                 * Multiplying by eff_dz gives a water depth [m]. */
                real wc_val = gw.wc(iGlob, 1); /* Current time step (index 1) */
                total_water_depth += wc_val * eff_dz;
                total_eff_depth += eff_dz;

                /* Safety check: prevent runaway accumulation due to
                 * unexpected grid geometry (should not normally trigger). */
                if (total_eff_depth > rd * 10.0)
                {
                    break;
                }
            }

            /* ------------------------------------------------------------- *
             * Step 4: Compute the depth-averaged root zone moisture content.
             *
             * root_zone_SM = total_water_depth / total_eff_depth
             *
             * This gives the effective volumetric water content experienced
             * by the crop roots, used by WOFOST to compute water stress.
             * ------------------------------------------------------------- */
            if (total_eff_depth > 0.0)
            {
                root_zone_SM(i_surf) = total_water_depth / total_eff_depth;
            }
            else
            {
                /* No overlap between root zone and any soil layer.
                 * Fall back to the top-layer moisture content. */
                int iGlob_top = (gdom.hc + 0) * gdom.nxhc * gdom.nyhc +
                                (gdom.hc + j) * gdom.nxhc +
                                (gdom.hc + i);
                if (iGlob_top >= 0 && iGlob_top < gw.wc.extent(0))
                {
                    root_zone_SM(i_surf) = gw.wc(iGlob_top, 1);
                }
                else
                {
                    root_zone_SM(i_surf) = 0.0;
                }
            }

            /* ------------------------------------------------------------- *
             * Debug/testing override (disabled by default).
             * Can be uncommented to force zero root zone moisture at a
             * specific column (i==1) to test water stress effects on
             * crop growth in a controlled scenario.
             * ------------------------------------------------------------- */
            // if(i == 1) {
            //     root_zone_SM(i_surf) = 0.0;
            // }
        });
    }

    /* ======================================================================= *
     * Crop Parameter Update for Groundwater Source/Sink Terms
     * ======================================================================= */

    /**
     * @brief Update the groundwater source/sink terms with current WOFOST
     *        crop parameters (LAI and root depth).
     *
     * This function transfers the spatially-distributed LAI and root depth
     * arrays from the WOFOST crop model to the SERGHEI groundwater model's
     * source/sink term objects. These parameters allow the groundwater model
     * to:
     *   - Partition total potential ET into soil evaporation and transpiration
     *     based on LAI (via an extinction coefficient model)
     *   - Apply root water extraction to the correct vertical layers
     *     based on root depth (RD)
     *
     * @param[in]     gdom           Groundwater domain with grid dimensions
     * @param[in,out] ss             SourceSink object containing groundwater
     *                               source/sink terms (will be modified)
     * @param[in]     LAI            Leaf area index array from WOFOST [m^2/m^2]
     *                               Typical range: 0 - 7, dimension: nCellSwSurface
     * @param[in]     RD             Root depth array from WOFOST [cm]
     *                               Typical range: 0 - 150 cm, dimension: nCellSwSurface
     * @param[in]     current_time   Current simulation time (unused in current
     *                               implementation but available for time-varying logic)
     *
     * @note The source/sink type 0 (sstype == 0) corresponds to evapotranspiration.
     *       Other types may include wells, drains, or boundary fluxes.
     */
    static void update_crop_parameters(const GwDomain &gdom, SourceSink &ss,
                                       const realArr &LAI, const realArr &RD,
                                       real current_time)
    {
        /* Iterate over all defined groundwater source/sink terms.
         * Typically, evapotranspiration (ET) is registered as a specific
         * source/sink type (type 0). */
        for (auto &gwss_item : ss.gwss)
        {
            if (gwss_item.sstype == 0)
            { /* Type 0 = Evapotranspiration */

                /* Mark this source/sink as using real-time WOFOST data,
                 * so the groundwater solver knows to apply LAI-based
                 * partitioning and root-depth-based extraction. */
                gwss_item.use_realtime_data = true;

                /* Total number of surface grid cells (including halo).
                 * nCellSwSurface = nCellSw + halo cells on all sides. */
                int nCellSwSurface = gdom.nxhc * gdom.nyhc;

                /* Allocate Kokkos device arrays for LAI and RD on the
                 * source/sink object if they have not been allocated yet.
                 * These persist across time steps. */
                if (gwss_item.wofost_lai.data() == nullptr)
                {
                    gwss_item.wofost_lai = realArr("wofost_lai", nCellSwSurface);
                }
                if (gwss_item.wofost_rd.data() == nullptr)
                {
                    gwss_item.wofost_rd = realArr("wofost_rd", nCellSwSurface);
                }

                /* Deep copy LAI and RD from the WOFOST arrays to the
                 * source/sink object's device arrays.
                 *
                 * Both arrays are sized for the full surface grid including
                 * halo cells. Halo values will typically be zero or
                 * extrapolated, depending on the boundary treatment. */
                Kokkos::deep_copy(gwss_item.wofost_lai, LAI);
                Kokkos::deep_copy(gwss_item.wofost_rd, RD);
            }
        }
    }

    /* ======================================================================= *
     * Evapotranspiration Partitioning (Soil Evaporation vs. Transpiration)
     * ======================================================================= */

    /**
     * @brief Partition potential evapotranspiration into soil evaporation and
     *        plant transpiration fractions based on leaf area index.
     *
     * Uses a simple exponential extinction model (Beer-Lambert law analogy)
     * to estimate how much of the net radiation reaches the soil surface
     * versus being intercepted by the canopy:
     *
     *   T_split = 1 - exp(-k_ext * LAI)   (transpiration fraction)
     *   E_split = exp(-k_ext * LAI)        (soil evaporation fraction)
     *
     * where k_ext is the extinction coefficient for global radiation.
     *
     * Physical interpretation:
     *   - At LAI = 0 (bare soil):  T_split = 0, E_split = 1 (all evaporation)
     *   - At LAI = 3:              T_split ~ 0.78, E_split ~ 0.22 (with k=0.5)
     *   - At LAI = 5+:             T_split ~ 0.92, E_split ~ 0.08 (mostly transpiration)
     *
     * @param[in]  gdom     Groundwater domain (for surface grid dimensions)
     * @param[in]  LAI      Leaf area index array [m^2/m^2], dimension: nCellSwSurface
     * @param[out] T_split  Transpiration fraction array [dimensionless, 0-1]
     * @param[out] E_split  Soil evaporation fraction array [dimensionless, 0-1]
     * @param[in]  k_ext    Extinction coefficient for global radiation [dimensionless].
     *                       Default: 0.5 (typical for many C3 crops).
     *                       Range: 0.3 (erectophile canopies) to 0.8 (planophile canopies).
     */
    static void compute_evap_partitioning(const GwDomain &gdom,
                                          const realArr &LAI,
                                          realArr &T_split, realArr &E_split,
                                          real k_ext = 0.5)
    {
        /* Total number of surface grid cells (including halo region). */
        int nCellSwSurface = gdom.nxhc * gdom.nyhc;

        Kokkos::parallel_for("CalcETSplit", nCellSwSurface, KOKKOS_LAMBDA(const int i_surf) {
            real lai = LAI(i_surf);
            real exp_term = exp(-k_ext * lai);

            /* Transpiration fraction: increases with LAI as the canopy
             * intercepts more radiation. */
            T_split(i_surf) = 1.0 - exp_term;

            /* Soil evaporation fraction: decreases with LAI as the canopy
             * shades the soil surface. */
            E_split(i_surf) = exp_term;
        });
    }
};

#endif
