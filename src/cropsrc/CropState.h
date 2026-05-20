/**
 * @file CropState.h
 * @brief Static crop parameters and parameter-table storage for the WOFOST model.
 *
 * This file defines:
 *   - CropParams: a struct holding all scalar crop parameters (constants for a
 *     given crop cultivar that do not change during simulation). Parameters are
 *     grouped by function: phenology, assimilation, respiration, root/soil,
 *     and nutrient (N/P/K) thresholds.
 *   - CropState: a wrapper class that owns the CropParams struct, the crop name,
 *     and a map of tabulated parameter functions (e.g., partitioning fractions
 *     as functions of development stage). Tables are stored as Kokkos Views
 *     (realArr) for GPU compatibility.
 *
 * Parameter units and typical ranges are documented inline below.
 * Reference: WOFOST 7.1 User Guide, Boogaard et al. (1998).
 *
 * @see Afgen.h     for the interpolation routine used on these tables.
 * @see CropInit.h  for the file reader that populates CropParams and tables.
 * @see Wofost72.h  for the main model class that consumes these parameters.
 */

/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _CROP_STATE_H_
#define _CROP_STATE_H_

#include "../const.h"
#include "../define.h"
#include "../SArray.h"
#include <map>
#include <string>
#include <iostream>

/**
 * @struct CropParams
 * @brief All scalar (non-tabulated) crop parameters for a single cultivar.
 *
 * These values are read once from the crop parameter input file and remain
 * constant throughout the simulation.  They define the physiological
 * characteristics of the crop being modelled.
 *
 * Conventions:
 *   [deg.C]   = degrees Celsius
 *   [deg.C.d] = degree-days (temperature sum)
 *   [d-1]     = per day
 *   [kg/ha]   = kilograms per hectare
 *   [m]       = metres
 *   [cm]      = centimetres
 *   [-]       = dimensionless
 */
struct CropParams {

    /* ================================================================== */
    /* Phenology: Emergence                                               */
    /* ================================================================== */

    real TBASEM;  /**< Lower temperature threshold for emergence [deg.C].
                        Typical range: 0.0 -- 10.0. Below this temperature,
                        emergence does not progress. */
    real TEFFMX;  /**< Maximum effective temperature for emergence [deg.C].
                        Typical range: 20.0 -- 40.0. Above this, the effective
                        temperature for emergence is capped. */
    real TSUMEM;  /**< Temperature sum required for emergence [deg.C.d].
                        Typical range: 50 -- 300. The crop emerges once the
                        accumulated effective temperature reaches this value. */

    /* ================================================================== */
    /* Phenology: Development Rate                                        */
    /* ================================================================== */

    int  IDSL;    /**< Development stage mode selector [-].
                        0 = not temperature- or daylength-driven (e.g., permanent pasture).
                        1 = temperature-driven only.
                        2 = both temperature- and daylength-driven. */
    real DLO;     /**< Optimal daylength for development [h].
                        Typical range: 10.0 -- 18.0. Only used when IDSL >= 2. */
    real DLC;     /**< Critical daylength for development [h].
                        Typical range: 8.0 -- 16.0. Only used when IDSL >= 2. */
    real TSUM1;   /**< Temperature sum from emergence to anthesis [deg.C.d].
                        Typical range: 500 -- 2500. */
    real TSUM2;   /**< Temperature sum from anthesis to maturity [deg.C.d].
                        Typical range: 300 -- 2000. */
    real DVSI;    /**< Initial development stage at emergence [-].
                        Typically 0.0. */
    real DVSEND;  /**< Development stage at which the crop is considered mature [-].
                        Typically 2.0 (WOFOST convention: 0--1 vegetative,
                        1--2 reproductive). */

    /* ================================================================== */
    /* Initial Conditions                                                 */
    /* ================================================================== */

    real TDWI;    /**< Initial total dry weight of the crop [kg/ha].
                        Typical range: 10 -- 300. Distributed among organs
                        via partitioning fractions at initialization. */
    real LAIEM;   /**< Leaf area index at emergence [ha/ha].
                        Typical range: 0.01 -- 0.5. */
    real RGRLAI;  /**< Maximum relative increase in LAI [ha/ha/d].
                        Typical range: 0.005 -- 0.05. Governs exponential
                        LAI growth in the juvenile stage. */

    /* ================================================================== */
    /* Leaf Dynamics                                                      */
    /* ================================================================== */

    real SPAN;    /**< Life span of leaves under optimal conditions [d].
                        Typical range: 200 -- 1000. After this age, leaves
                        senesce at the maximum rate. */
    real TBASE;   /**< Lower threshold temperature for ageing of leaves [deg.C].
                        Typical range: 0.0 -- 10.0. */
    real SPA;     /**< Specific pod area [m2/kg]. Used when PAI > 0.
                        Typical range: 0.0 -- 5.0. */

    /* ================================================================== */
    /* Conversion Efficiency (Biomass Conversion Factors)                 */
    /* ================================================================== */

    real CVL;     /**< Efficiency of conversion of assimilates to leaf
                        biomass [kg/kg]. Typical range: 0.60 -- 0.75.
                        (= 1 / cost in CH2O per kg leaf dry matter). */
    real CVO;     /**< Efficiency of conversion of assimilates to storage
                        organ biomass [kg/kg]. Typical range: 0.65 -- 0.85. */
    real CVR;     /**< Efficiency of conversion of assimilates to root
                        biomass [kg/kg]. Typical range: 0.60 -- 0.75. */
    real CVS;     /**< Efficiency of conversion of assimilates to stem
                        biomass [kg/kg]. Typical range: 0.60 -- 0.75. */

    /* ================================================================== */
    /* Maintenance Respiration                                            */
    /* ================================================================== */

    real Q10;     /**< Relatively increase in respiration rate per 10 deg.C
                        temperature increase [-]. Typical value: 2.0. */
    real RML;     /**< Relative maintenance respiration rate of leaves
                        [kg CH2O/kg DM/d]. Typical range: 0.010 -- 0.040. */
    real RMO;     /**< Relative maintenance respiration rate of storage organs
                        [kg CH2O/kg DM/d]. Typical range: 0.001 -- 0.010. */
    real RMR;     /**< Relative maintenance respiration rate of roots
                        [kg CH2O/kg DM/d]. Typical range: 0.005 -- 0.025. */
    real RMS;     /**< Relative maintenance respiration rate of stems
                        [kg CH2O/kg DM/d]. Typical range: 0.005 -- 0.025. */

    /* ================================================================== */
    /* Senescence & Mortality                                             */
    /* ================================================================== */

    real PERDL;   /**< Maximum relative death rate of leaves due to water
                        stress [d-1]. Typical range: 0.01 -- 0.10. */

    /* ================================================================== */
    /* Evapotranspiration                                                 */
    /* ================================================================== */

    real CFET;    /**< Crop factor for potential evapotranspiration [-].
                        Typical range: 0.8 -- 1.2. Scales potential ET. */
    real DEPNR;   /**< Dependency of crop water use on soil moisture
                        (critical soil moisture deficit) [-].
                        Typical range: 2.5 -- 6.0. Higher values indicate
                        greater sensitivity to drought. */

    /* ================================================================== */
    /* Root Dynamics                                                      */
    /* ================================================================== */

    int  IAIRDU;  /**< Switch for air ducts in roots (aeration) [-].
                        0 = no air ducts, 1 = air ducts present. */
    int  IOX;     /**< Switch for oxygen stress calculation [-].
                        0 = no oxygen stress, 1 = oxygen stress applied. */
    real RDI;     /**< Initial rooting depth [cm]. Typical range: 5 -- 30. */
    real RRI;     /**< Maximum daily increase in rooting depth [cm/d].
                        Typical range: 1.0 -- 5.0. */
    real RDMCR;   /**< Maximum rooting depth of the crop [cm].
                        Typical range: 50 -- 200. */

    /* ================================================================== */
    /* Soil Physical Characteristics                                      */
    /* ================================================================== */

    real SMW;     /**< Soil moisture content at wilting point [cm3/cm3].
                        Typical range: 0.05 -- 0.20. */
    real SMFCF;   /**< Soil moisture content at field capacity [cm3/cm3].
                        Typical range: 0.15 -- 0.45. */
    real SM0;     /**< Soil porosity (saturation moisture content) [cm3/cm3].
                        Typical range: 0.35 -- 0.60. */
    real CRAIRC;  /**< Critical air content for aeration stress [cm3/cm3].
                        Typical range: 0.02 -- 0.10. If air-filled porosity
                        drops below this, root functioning is impaired. */

    /* ================================================================== */
    /* Vernalization                                                      */
    /* ================================================================== */

    real VERNBASE; /**< Base vernalization temperature [deg.C].
                         Typical range: -5.0 -- 5.0. */
    real VERNSAT;  /**< Saturated (full) vernalization temperature [deg.C].
                         Typical range: 5.0 -- 15.0. */

    /* ================================================================== */
    /* Nutrient Parameters (N, P, K)                                      */
    /* ================================================================== */

    /* --- Nitrogen [kg N / kg DM] --- */
    real NMINSO;  /**< Minimum N concentration in storage organs
                        [kg N/kg DM]. Typical range: 0.005 -- 0.015. */
    real NMINVE;  /**< Minimum N concentration in vegetative organs
                        [kg N/kg DM]. Typical range: 0.005 -- 0.020. */
    real NMAXSO;  /**< Maximum N concentration in storage organs
                        [kg N/kg DM]. Typical range: 0.010 -- 0.040. */
    real NMAXVE;  /**< Maximum N concentration in vegetative organs
                        [kg N/kg DM]. Typical range: 0.020 -- 0.060. */

    /* --- Phosphorus [kg P / kg DM] --- */
    real PMINSO;  /**< Minimum P concentration in storage organs
                        [kg P/kg DM]. Typical range: 0.001 -- 0.004. */
    real PMINVE;  /**< Minimum P concentration in vegetative organs
                        [kg P/kg DM]. Typical range: 0.001 -- 0.004. */
    real PMAXSO;  /**< Maximum P concentration in storage organs
                        [kg P/kg DM]. Typical range: 0.002 -- 0.008. */
    real PMAXVE;  /**< Maximum P concentration in vegetative organs
                        [kg P/kg DM]. Typical range: 0.002 -- 0.008. */

    /* --- Potassium [kg K / kg DM] --- */
    real KMINSO;  /**< Minimum K concentration in storage organs
                        [kg K/kg DM]. Typical range: 0.003 -- 0.010. */
    real KMINVE;  /**< Minimum K concentration in vegetative organs
                        [kg K/kg DM]. Typical range: 0.003 -- 0.015. */
    real KMAXSO;  /**< Maximum K concentration in storage organs
                        [kg K/kg DM]. Typical range: 0.005 -- 0.020. */
    real KMAXVE;  /**< Maximum K concentration in vegetative organs
                        [kg K/kg DM]. Typical range: 0.008 -- 0.030. */

    /* --- General Nutrient --- */
    real YZERO;   /**< Residual (zero-yield) biomass fraction [-].
                        Typical range: 0.0 -- 0.1. */
    real NFIX;    /**< Biological nitrogen fixation rate [kg N/ha/d].
                        Typical range: 0.0 -- 3.0 (significant for legumes). */
};

/**
 * @class CropState
 * @brief Container for all static crop data: scalar parameters, name, and
 *        tabulated parameter functions.
 *
 * The tabulated functions are stored in a std::map<string, realArr> where each
 * realArr is a Kokkos View holding flat [x0, y0, x1, y1, ...] data. These
 * tables are populated during initialization (CropInit::readCropParameters)
 * and are read (via Afgen::lookup) during the simulation kernel.
 *
 * IMPORTANT: The std::map is host-only. Device kernels must receive the
 * individual Kokkos Views (extracted via getTable) as capture-by-value
 * arguments; they cannot access the map directly.
 */
class CropState {
public:
    CropParams p;                    /**< Scalar crop parameters (see struct above). */
    std::string cropName;            /**< Crop / cultivar name from input file. */

    /** @brief Map of tabulated parameter functions, keyed by WOFOST parameter name.
     *  Example keys: "DTSMTB", "FLTB", "FRTB", "SLATB", "AMAXTB", etc. */
    std::map<std::string, realArr> tables;

    /**
     * @brief Allocate a Kokkos View for a named parameter table.
     *
     * @param name  Parameter table name (e.g., "DTSMTB").
     * @param size  Number of elements in the flat array (must be even).
     *
     * @note The size parameter is size_t to match Kokkos dimension types.
     */
    void allocateTable(std::string name, size_t size) {
        tables[name] = realArr(name, size);
    }

    /**
     * @brief Check whether a table with the given name has been loaded.
     * @param name  Parameter table name.
     * @return true if the table exists in the map.
     */
    bool hasTable(std::string name) const {
        return tables.find(name) != tables.end();
    }

    /**
     * @brief Retrieve a Kokkos View for a named parameter table.
     *
     * @param name  Parameter table name.
     * @return      The Kokkos View containing the tabulated data.
     *
     * @warning This method uses std::map, which is NOT accessible on the GPU
     *          device. Only call this from host code (e.g., during
     *          initialization or when extracting Views to pass to kernels).
     */
    realArr getTable(std::string name) {
        if (tables.find(name) != tables.end()) {
            return tables[name];
        } else {
            std::cerr << "Error: Table " << name << " not found." << std::endl;
            return realArr();  /* Return an empty (unallocated) View */
        }
    }
};

#endif /* _CROP_STATE_H_ */
