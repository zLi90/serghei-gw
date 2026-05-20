/**
 * @file CropInit.h
 * @brief Crop parameter file reader for the WOFOST model.
 *
 * Reads a colon-delimited text file containing scalar crop parameters and
 * tabulated parameter functions. Each line has the format:
 *
 *   KEY : value1 value2 ...  // optional comment
 *
 * Scalar parameters are stored directly in CropParams. Tabulated parameters
 * (identified by key name) are stored as flat [x0, y0, x1, y1, ...] arrays
 * in CropState::tables via Kokkos Views, with host-to-device deep copy.
 *
 * Supported table keys:
 *   DTSMTB  - Effective temperature as f(mean temperature) [deg.C vs deg.C]
 *   SLATB   - Specific leaf area as f(DVS) [m2/kg vs -]
 *   SSATB   - Specific stem area as f(DVS) [m2/kg vs -]
 *   KDIFTB  - Extinction coefficient for diffuse visible light as f(DVS) [- vs -]
 *   EFFTB   - Light-use efficiency as f(mean temperature) [kg HA/(MJ) vs deg.C]
 *   AMAXTB  - Maximum assimilation rate as f(DVS) [kg/ha/h vs -]
 *   TMPFTB  - Temperature reduction factor for AMAX as f(mean temp) [- vs deg.C]
 *   TMNFTB  - Reduction factor for low minimum temperature as f(Tmin) [- vs deg.C]
 *   RFSETB  - Reduction factor for transpiration as f(SM) [- vs cm3/cm3]
 *   FRTB    - Fraction of biomass to roots as f(DVS) [- vs -]
 *   FLTB    - Fraction of above-ground biomass to leaves as f(DVS) [- vs -]
 *   FSTB    - Fraction of above-ground biomass to stems as f(DVS) [- vs -]
 *   FOTB    - Fraction of above-ground biomass to storage organs as f(DVS) [- vs -]
 *   RDRRTB  - Relative death rate of roots as f(DVS) [d-1 vs -]
 *   RDRSTB  - Relative death rate of stems as f(DVS) [d-1 vs -]
 *
 * @see CropState.h  for the data structures populated by this reader.
 * @see Afgen.h      for the interpolation function used on the tables.
 */

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

class CropInit
{
    /* ================================================================== */
    /* Internal line parser class                                         */
    /* ================================================================== */

    /**
     * @class PsLn
     * @brief Parses a single line from the crop parameter file.
     *
     * Splits a line into a key (before the colon) and a value stream
     * (after the colon). Comment lines (starting with "//") produce
     * empty keys and are silently skipped.
     */
    class PsLn
    {
    public:
        std::string line;       /**< Raw input line. */
        std::string key;        /**< Extracted parameter key (trimmed, no whitespace). */
        std::stringstream value;/**< Stream of values after the colon, comments stripped. */

        /** @brief Convert the entire line to lowercase (currently unused
         *         because WOFOST keys are typically uppercase). */
        void lowercase()
        {
            std::for_each(line.begin(), line.end(), [](char &c)
                          { c = ::tolower(c); });
        }

        /**
         * @brief Parse the line into key and value.
         *
         * Lines starting with "//" or empty lines are ignored.
         * The key is the substring before the first colon, with all spaces
         * and tabs removed. The value is the substring after the colon,
         * with any trailing "//..." comment removed.
         */
        void parse()
        {
            key.clear();
            value.clear();

            /* Skip empty lines and comment-only lines */
            if (!line.empty() && line.find("//", 0) != 0)
            {
                /* Find the colon separator */
                size_t splitloc = line.find(':', 0);
                if (splitloc != std::string::npos)
                {
                    /* Extract key: everything before the colon, strip whitespace */
                    key = line.substr(0, splitloc);
                    key.erase(std::remove(key.begin(), key.end(), ' '), key.end());
                    key.erase(std::remove(key.begin(), key.end(), '\t'), key.end());

                    /* Extract value: everything after the colon, strip inline comment */
                    std::string val = line.substr(splitloc + 1, line.length() - splitloc);
                    size_t splitter = val.find("//", 0);
                    std::string strloc;
                    if (splitter != std::string::npos)
                    {
                        strloc = val.substr(0, splitter);
                    }
                    else
                    {
                        strloc = val;
                    }
                    value.clear();
                    value.str(strloc);
                }
            }
        }
    };

public:
    /**
     * @brief Read all crop parameters from an input file.
     *
     * Opens the specified file, parses each line, and populates the scalar
     * parameters in CropState::p and the tabulated functions in
     * CropState::tables. Table data is deep-copied from host to device
     * via Kokkos::create_mirror_view / Kokkos::deep_copy.
     *
     * @param fNameIn  Path to the crop parameter input file.
     * @param cs       CropState object to populate.
     * @param par      Parallel controller (for masterproc error reporting).
     * @return 1 on success, 0 on failure (file not found or parse error).
     */
    int readCropParameters(std::string fNameIn, CropState &cs, Parallel &par)
    {
        std::ifstream fInStream(fNameIn);
        std::string line;
        PsLn pline;

        /* --- Open the parameter file --- */
        if (!fInStream.is_open())
        {
            if (par.masterproc)
                std::cerr << RERROR "Unable to open crop parameter file: " << fNameIn << "\n";
            return 0;
        }

        /* --- Parse each line --- */
        while (std::getline(fInStream, line))
        {
            pline.line = line;
            pline.parse();

            if (!pline.key.empty())
            {
                /* ====================================================== */
                /* Crop name                                              */
                /* ====================================================== */
                if (!strcmp("CROPNAM", pline.key.c_str()))
                {
                    pline.value >> cs.cropName;
                }

                /* ====================================================== */
                /* Phenology: Emergence parameters                        */
                /* ====================================================== */
                else if (!strcmp("TBASEM", pline.key.c_str()))
                {
                    pline.value >> cs.p.TBASEM;
                }
                else if (!strcmp("TEFFMX", pline.key.c_str()))
                {
                    pline.value >> cs.p.TEFFMX;
                }
                else if (!strcmp("TSUMEM", pline.key.c_str()))
                {
                    pline.value >> cs.p.TSUMEM;
                }

                /* ====================================================== */
                /* Phenology: Development parameters                      */
                /* ====================================================== */
                else if (!strcmp("IDSL", pline.key.c_str()))
                {
                    pline.value >> cs.p.IDSL;
                }
                else if (!strcmp("DLO", pline.key.c_str()))
                {
                    pline.value >> cs.p.DLO;
                }
                else if (!strcmp("DLC", pline.key.c_str()))
                {
                    pline.value >> cs.p.DLC;
                }
                else if (!strcmp("TSUM1", pline.key.c_str()))
                {
                    pline.value >> cs.p.TSUM1;
                }
                else if (!strcmp("TSUM2", pline.key.c_str()))
                {
                    pline.value >> cs.p.TSUM2;
                }
                else if (!strcmp("DVSI", pline.key.c_str()))
                {
                    pline.value >> cs.p.DVSI;
                }
                else if (!strcmp("DVSEND", pline.key.c_str()))
                {
                    pline.value >> cs.p.DVSEND;
                }

                /* ====================================================== */
                /* Initial conditions                                     */
                /* ====================================================== */
                else if (!strcmp("TDWI", pline.key.c_str()))
                {
                    pline.value >> cs.p.TDWI;
                }
                else if (!strcmp("LAIEM", pline.key.c_str()))
                {
                    pline.value >> cs.p.LAIEM;
                }
                else if (!strcmp("RGRLAI", pline.key.c_str()))
                {
                    pline.value >> cs.p.RGRLAI;
                }

                /* ====================================================== */
                /* Leaf dynamics                                          */
                /* ====================================================== */
                else if (!strcmp("SPAN", pline.key.c_str()))
                {
                    pline.value >> cs.p.SPAN;
                }
                else if (!strcmp("TBASE", pline.key.c_str()))
                {
                    pline.value >> cs.p.TBASE;
                }
                else if (!strcmp("SPA", pline.key.c_str()))
                {
                    pline.value >> cs.p.SPA;
                }

                /* ====================================================== */
                /* Conversion efficiency (biomass conversion factors)     */
                /* ====================================================== */
                else if (!strcmp("CVL", pline.key.c_str()))
                {
                    pline.value >> cs.p.CVL;
                }
                else if (!strcmp("CVO", pline.key.c_str()))
                {
                    pline.value >> cs.p.CVO;
                }
                else if (!strcmp("CVR", pline.key.c_str()))
                {
                    pline.value >> cs.p.CVR;
                }
                else if (!strcmp("CVS", pline.key.c_str()))
                {
                    pline.value >> cs.p.CVS;
                }

                /* ====================================================== */
                /* Maintenance respiration                                */
                /* ====================================================== */
                else if (!strcmp("Q10", pline.key.c_str()))
                {
                    pline.value >> cs.p.Q10;
                }
                else if (!strcmp("RML", pline.key.c_str()))
                {
                    pline.value >> cs.p.RML;
                }
                else if (!strcmp("RMO", pline.key.c_str()))
                {
                    pline.value >> cs.p.RMO;
                }
                else if (!strcmp("RMR", pline.key.c_str()))
                {
                    pline.value >> cs.p.RMR;
                }
                else if (!strcmp("RMS", pline.key.c_str()))
                {
                    pline.value >> cs.p.RMS;
                }

                /* ====================================================== */
                /* Senescence & mortality                                 */
                /* ====================================================== */
                else if (!strcmp("PERDL", pline.key.c_str()))
                {
                    pline.value >> cs.p.PERDL;
                }

                /* ====================================================== */
                /* Evapotranspiration                                     */
                /* ====================================================== */
                else if (!strcmp("CFET", pline.key.c_str()))
                {
                    pline.value >> cs.p.CFET;
                }
                else if (!strcmp("DEPNR", pline.key.c_str()))
                {
                    pline.value >> cs.p.DEPNR;
                }

                /* ====================================================== */
                /* Root dynamics                                          */
                /* ====================================================== */
                else if (!strcmp("IAIRDU", pline.key.c_str()))
                {
                    pline.value >> cs.p.IAIRDU;
                }
                else if (!strcmp("IOX", pline.key.c_str()))
                {
                    pline.value >> cs.p.IOX;
                }
                else if (!strcmp("RDI", pline.key.c_str()))
                {
                    pline.value >> cs.p.RDI;
                }
                else if (!strcmp("RRI", pline.key.c_str()))
                {
                    pline.value >> cs.p.RRI;
                }
                else if (!strcmp("RDMCR", pline.key.c_str()))
                {
                    pline.value >> cs.p.RDMCR;
                }

                /* ====================================================== */
                /* Nutrient parameters (N, P, K)                          */
                /* ====================================================== */
                else if (!strcmp("NMINSO", pline.key.c_str()))
                {
                    pline.value >> cs.p.NMINSO;
                }
                else if (!strcmp("NMINVE", pline.key.c_str()))
                {
                    pline.value >> cs.p.NMINVE;
                }
                else if (!strcmp("NMAXSO", pline.key.c_str()))
                {
                    pline.value >> cs.p.NMAXSO;
                }
                else if (!strcmp("NMAXVE", pline.key.c_str()))
                {
                    pline.value >> cs.p.NMAXVE;
                }
                else if (!strcmp("PMINSO", pline.key.c_str()))
                {
                    pline.value >> cs.p.PMINSO;
                }
                else if (!strcmp("PMINVE", pline.key.c_str()))
                {
                    pline.value >> cs.p.PMINVE;
                }
                else if (!strcmp("PMAXSO", pline.key.c_str()))
                {
                    pline.value >> cs.p.PMAXSO;
                }
                else if (!strcmp("PMAXVE", pline.key.c_str()))
                {
                    pline.value >> cs.p.PMAXVE;
                }
                else if (!strcmp("KMINSO", pline.key.c_str()))
                {
                    pline.value >> cs.p.KMINSO;
                }
                else if (!strcmp("KMINVE", pline.key.c_str()))
                {
                    pline.value >> cs.p.KMINVE;
                }
                else if (!strcmp("KMAXSO", pline.key.c_str()))
                {
                    pline.value >> cs.p.KMAXSO;
                }
                else if (!strcmp("KMAXVE", pline.key.c_str()))
                {
                    pline.value >> cs.p.KMAXVE;
                }
                else if (!strcmp("YZERO", pline.key.c_str()))
                {
                    pline.value >> cs.p.YZERO;
                }
                else if (!strcmp("NFIX", pline.key.c_str()))
                {
                    pline.value >> cs.p.NFIX;
                }

                /* ====================================================== */
                /* Soil physical characteristics                          */
                /* ====================================================== */
                else if (!strcmp("SMW", pline.key.c_str()))
                {
                    pline.value >> cs.p.SMW;
                }
                else if (!strcmp("SMFCF", pline.key.c_str()))
                {
                    pline.value >> cs.p.SMFCF;
                }
                else if (!strcmp("SM0", pline.key.c_str()))
                {
                    pline.value >> cs.p.SM0;
                }
                else if (!strcmp("CRAIRC", pline.key.c_str()))
                {
                    pline.value >> cs.p.CRAIRC;
                }

                /* ====================================================== */
                /* Tabulated parameter functions                          */
                /*                                                       */
                /* Values are read into a temporary host vector, then     */
                /* copied into a Kokkos View and deep-copied to device.   */
                /* Table format: x0 y0 x1 y1 ... (even number of values) */
                /* ====================================================== */
                else if (isTableKey(pline.key))
                {
                    std::vector<real> values;
                    real val;

                    /* Read all numeric values after the colon */
                    while (pline.value >> val)
                    {
                        values.push_back(val);
                    }

                    /* Tables must have an even number of elements (x,y pairs) */
                    if (values.size() % 2 != 0)
                    {
                        if (par.masterproc)
                            std::cerr << YEXC "Warning: Odd number of elements in table " << pline.key << ". Check input file.\n";
                    }

                    /* Step 1: Allocate the Kokkos View in CropState */
                    cs.allocateTable(pline.key, values.size());

                    /* Step 2: Get a reference to the allocated View */
                    realArr tableRef = cs.getTable(pline.key);

                    /* Step 3: Create a host mirror and fill with parsed data */
                    realArr::HostMirror h_table = Kokkos::create_mirror_view(tableRef);

                    for (size_t i = 0; i < values.size(); ++i)
                    {
                        h_table(i) = values[i];
                    }

                    /* Step 4: Deep copy from host mirror to device View */
                    Kokkos::deep_copy(tableRef, h_table);
                }
            }
        }

        /* --- Report success --- */
        if (par.masterproc)
            std::cout << GOK "Crop parameters read from " << fNameIn << "\n";
        return 1;
    }

private:
    /**
     * @brief Check whether a key corresponds to a tabulated parameter function.
     *
     * Tabulated parameters are stored as flat arrays of (x, y) pairs rather
     * than single scalar values. This method checks the key against the
     * known list of WOFOST table parameter names.
     *
     * @param key  Parameter name from the input file.
     * @return     true if key is a table parameter name.
     */
    bool isTableKey(const std::string &key)
    {
        static const std::vector<std::string> tableKeys = {
            "DTSMTB",   /**< Effective temperature as f(mean temperature) */
            "SLATB",    /**< Specific leaf area as f(DVS) */
            "SSATB",    /**< Specific stem area as f(DVS) */
            "KDIFTB",   /**< Diffuse light extinction coefficient as f(DVS) */
            "EFFTB",    /**< Light-use efficiency as f(mean temperature) */
            "AMAXTB",   /**< Maximum CO2 assimilation rate as f(DVS) */
            "TMPFTB",   /**< Temperature reduction factor for AMAX as f(Tmean) */
            "TMNFTB",   /**< Low-temperature reduction factor as f(Tmin) */
            "RFSETB",   /**< Transpiration reduction factor as f(soil moisture) */
            "FRTB",     /**< Root partitioning fraction as f(DVS) */
            "FLTB",     /**< Leaf partitioning fraction as f(DVS) */
            "FSTB",     /**< Stem partitioning fraction as f(DVS) */
            "FOTB",     /**< Storage organ partitioning fraction as f(DVS) */
            "RDRRTB",   /**< Root death rate as f(DVS) */
            "RDRSTB"    /**< Stem death rate as f(DVS) */
        };
        return std::find(tableKeys.begin(), tableKeys.end(), key) != tableKeys.end();
    }
};

#endif /* _CROP_INIT_H_ */
