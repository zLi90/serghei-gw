/* -*- mode: c++; c-default-style: "linux" -*- */
/*******************************************************************************
 * @file    RTInitSW.h
 * @brief   Initialization routines for Surface Water Reactive Transport Model (SW-RTM)
 *
 * This header defines the RTInitSW class, which extends the Initializer base class
 * to handle all initialization steps for the surface water reactive transport module.
 * The module supports multiple dissolved species (solutes) and provides:
 *
 *   - Parsing of key:value input files (via the PsLn helper class)
 *   - Reading transport parameters (diffusion, dispersion, advection scheme)
 *   - Reading boundary conditions (constant, zero-gradient, time-varying)
 *   - Reading source/sink terms (with polygon-based spatial selection)
 *   - Reading reaction parameters (decay, nitrification, denitrification,
 *     aerobic respiration, mineralization, and surface-water-specific processes)
 *   - Setting initial solute concentrations from constants or spatial data files
 *
 * Supported input file formats:
 *   - Simplified ASCII grid (3-line header: ncols, nrows, nodata)
 *   - ESRI ASCII Grid (6-line header: NCOLS, NROWS, XLLCORNER, YLLCORNER,
 *     CELLSIZE, NODATA_value)
 *
 * @note    Surface water RTM does not consider a solid phase (unlike groundwater
 *          RTM), so only aqueous-phase concentrations and parameters are managed.
 *
 * @see     RTInitGW.h   - Groundwater RTM initialization (reference implementation)
 * @see     RTStateSW.h  - Surface water reactive transport state variables
 ******************************************************************************/

#ifndef _RT_INIT_SW_H_
#define _RT_INIT_SW_H_

#include "define.h"
#include "Domain.h"
#include "BC.h"
#include "FileIO.h"
#include "State.h"
#include "Parallel.h"
#include "Parser.h"
#include "SourceSink.h"
#include "RTInitGW.h"
#include "RTStateSW.h"

/* ========================================================================== */
/*                          CLASS: RTInitSW                                    */
/* ========================================================================== */

/**
 * @class RTInitSW
 * @brief Initializes the Surface Water Reactive Transport Model (SW-RTM).
 *
 * Extends the Initializer base class. Orchestrates the full initialization
 * pipeline for surface water reactive transport: reading parameter files,
 * setting initial conditions, applying boundary conditions, and writing
 * initial output. Supports multiple solute species with independent
 * transport and reaction parameters per species.
 */
class RTInitSW : public Initializer
{

public:
    /* ====================================================================== */
    /*                     NESTED CLASS: PsLn (Parser Line)                    */
    /* ====================================================================== */

    /**
     * @class PsLn
     * @brief Lightweight line parser for colon-separated key:value input files.
     *
     * Splits an input line into a key (before the colon) and a value (after
     * the colon). Handles comment lines (starting with "//"), strips whitespace
     * and tabs from keys, replaces semicolons with spaces in values for
     * downstream stringstream parsing, and trims inline comments ("//...").
     *
     * Example input line:
     *   n_mass: 5  // number of species
     * Parses to: key="n_mass", value="5"
     */
    class PsLn
    {
    public:
        std::string line;       ///< The raw input line to parse
        std::string key;        ///< Extracted key (left of the colon, whitespace stripped)
        std::stringstream value;///< Extracted value (right of the colon, comments removed)

        /**
         * @brief Convert the raw line to lowercase in-place.
         *
         * Useful for case-insensitive key matching. Applies ::tolower to
         * every character in the line string.
         */
        void lowercase()
        {
            std::for_each(line.begin(), line.end(), [](char &c)
                          { c = ::tolower(c); });
        }

        /**
         * @brief Parse the raw line into key and value components.
         *
         * Processing steps:
         *   1. Clear any previous key/value (safe for reuse of PsLn objects).
         *   2. Skip empty lines and comment lines (starting with "//").
         *   3. Split on the first colon (':') to separate key from value.
         *   4. Strip spaces and tabs from the key.
         *   5. Trim inline comments ("//...") from the value portion.
         *   6. Replace semicolons (';') with spaces to support multi-value
         *      parsing via stringstream (e.g., arrays of parameters).
         */
        void parse()
        {
            // Clear previous results to allow reuse of this PsLn object
            key.clear();
            value.clear();

            // Only process non-empty, non-comment lines
            if (!line.empty() && line.find("//", 0) != 0)
            {
                // Split at the first colon
                uint splitloc = line.find(':', 0);
                key = line.substr(0, splitloc);

                // Remove all spaces and tabs from the key
                key.erase(std::remove(key.begin(), key.end(), ' '), key.end());
                key.erase(std::remove(key.begin(), key.end(), '\t'), key.end());

                // Extract the value portion (everything after the colon)
                std::string val = line.substr(splitloc + 1, line.length() - splitloc);

                // Trim inline comments ("// ...") from the value
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

                // Replace semicolons with spaces for proper stringstream parsing
                std::replace(strloc.begin(), strloc.end(), ';', ' ');
                value.str(strloc);
            }
        }

        /**
         * @brief Print the parsed line, key, and value to stdout (for debugging).
         */
        void print() { std::cout << "line: " << line << "\tkey: " << key << "\tvalue: " << value.str() << std::endl; }
    };

public:
    /* ====================================================================== */
    /*              MAIN INITIALIZATION: initialize_rt()                       */
    /* ====================================================================== */

    /**
     * @brief Main entry point for surface water reactive transport initialization.
     *
     * Orchestrates the full initialization pipeline in the correct order:
     *   1. Read transport parameters from rttransportsw.input (determines n_mass)
     *   2. Allocate state arrays for all species (requires n_mass)
     *   3. Re-apply transport parameters to Kokkos views (after allocation)
     *   4. Read boundary conditions from rtswbc.input
     *   5. Read source/sink terms from rtswss.input
     *   6. Read reaction parameters from reactionsw.input (if ReactionModule enabled)
     *   7. Set initial solute concentrations from constants or files
     *   8. Apply boundary conditions to all species
     *   9. Copy c(1) -> c(0) so both time levels have valid initial data
     *  10. Write initial output
     *
     * @param rtsw       Surface water reactive transport state [modifed]
     * @param state      Simulation state (flow variables) [read-only]
     * @param dom        Computational domain (grid, indices) [read-only]
     * @param ss         Source/sink data container [modified]
     * @param ebc        External boundary conditions [modified]
     * @param io         File I/O handler for output [read-only]
     * @param par        Parallel decomposition info [read-only]
     * @param inFolder   Input directory path (trailing slash expected)
     * @param outFolder  Output directory path (trailing slash expected)
     *
     * @return 1 on success, 0 on failure
     */
    int initialize_rt(RTStateSW &rtsw, State &state, Domain &dom, SourceSink &ss, ExternalBoundaries &ebc, FileIO &io, Parallel &par, std::string inFolder, std::string outFolder)
    {

        // --- Step 1: Read transport parameters (must be first to determine n_mass) ---
        std::string fNameIn = inFolder + "rttransportsw.input";

        if (!readRTFileSW(fNameIn, rtsw, par))
        {
            std::cerr << GOK << " SWRTM: Reading in SW reactive transport  input file failed." << std::endl;
            return 0;
        }

        // --- Step 2: Allocate surface water state arrays (requires n_mass) ---
        rtsw.allocate(dom);

        // --- Step 3: Re-apply transport parameters after allocation ---
        // Transport parameters were read into temporary std::vector storage before
        // allocation. Now copy them into the Kokkos views that were created by
        // allocate(). This two-step approach is necessary because Kokkos views
        // cannot be sized until n_mass is known and allocate() is called.
        for (int i = 0; i < rtsw.n_mass; i++)
        {
            rtsw.diffusion_molecular(i) = rtsw.temp_diffusion_molecular[i];  ///< Molecular diffusion coefficient [m^2/s]
            rtsw.alpha_L(i) = rtsw.temp_alpha_L[i];                          ///< Longitudinal dispersivity [m]
            rtsw.alpha_T(i) = rtsw.temp_alpha_T[i];                          ///< Transverse dispersivity [m]
            rtsw.RainCon_arr(i) = rtsw.RainCon[i];                           ///< Rainfall concentration [mg/L or kg/m^3]
        }

        // --- Step 4: Read boundary conditions ---
        fNameIn = inFolder + "rtswbc.input";
        if (!readRTBCFileSW(fNameIn, dom, ebc, par, state, rtsw))
        {

            if (par.masterproc)
            {
                std::cerr << RERROR "SWRTM: Unable to read rtsw BC SW File " << std::endl;
                return 0;
            }
        }

        // --- Step 5: Read source/sink terms ---
        fNameIn = inFolder + "rtswss.input";
        if (!readRTSSFileSW(fNameIn, dom, ss.swss, par, state, rtsw))
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "SWRTM: Error encountered reading rtswss.input file " << std::endl;
                return 0;
            }
        }

        // --- Step 6: Read reaction parameters (only if reaction module is enabled) ---
        if (rtsw.ReactionModule == 1)
        {
            fNameIn = inFolder + "reactionsw.input";
            if (!readRTReactionFile(fNameIn, rtsw, par))
            {
                if (par.masterproc)
                {
                    std::cerr << RERROR "SWRTM: Unable to read reactionsw.input file" << std::endl;
                    return 0;
                }
            }
        }

        // --- Step 7: Set initial solute concentrations ---
        fNameIn = inFolder;
        if (!setRtStateSW(fNameIn, rtsw, dom, par))
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "SWRTM: Unable to set SW solute IC " << std::endl;
                return 0;
            }
        }

        // --- Step 8: Apply boundary conditions to all species ---
        // Each boundary's applyrtbc() internally iterates over all species
        // (iSpec = 0 to n_mass-1) and applies species-specific boundary conditions.
        for (int k = 0; k < ebc.rtbc.size(); k++)
        {
            ebc.rtbc[k].applyrtbc(rtsw, dom, state);
        }
        if (par.masterproc)
        {
            std::cerr << GOK "SWRTM: SW reactive transport boundary conditions applied successfully" << std::endl;
        }

        // --- Step 9: Copy time level 1 concentrations to time level 0 ---
        // Initialize c(0) from c(1) so that both the current and previous time
        // levels contain valid concentration data at the start of the simulation.
        Kokkos::parallel_for("update", dom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
            for (int iSpec = 0; iSpec < rtsw.n_mass; iSpec++)
            {
                rtsw.c(iSpec, iGlob, 0) = rtsw.c(iSpec, iGlob, 1);
            } });

        // --- Step 10: Write initial output ---
        io.outputIniRTSW(rtsw, state, dom, par, outFolder);

        return 1;
    }

    /* ====================================================================== */
    /*              SET INITIAL STATE: setRtStateSW()                          */
    /* ====================================================================== */

    /**
     * @brief Set initial aqueous-phase solute concentrations for all species.
     *
     * For each species (iSpec = 0 .. n_mass-1), the initial concentration
     * field is set according to the species' initialization mode:
     *   - Mode 0 (constant): All cells receive the same concentration value.
     *   - Mode 1 (file):     Concentrations are read from a spatial data file
     *                         (simplified or ESRI ASCII Grid format).
     *
     * Both time levels c(iSpec, *, 0) and c(iSpec, *, 1) are initialized.
     *
     * @param inFolder  Input directory path containing IC files
     * @param rtsw      Surface water RT state [modified]
     * @param dom       Computational domain [read-only]
     * @param par       Parallel decomposition info [read-only]
     *
     * @return 1 on success; exits with error code on failure
     */
    int setRtStateSW(std::string inFolder, RTStateSW &rtsw, Domain &dom, Parallel &par)
    {
        int ii, jj, kk, idx, ivg, iGlob, iGlobSW;
        std::string tempStr;

        // Zero-initialize all concentration arrays for all species and time levels
        for (iGlob = 0; iGlob < dom.nCellMem; iGlob++)
        {
            for (int iSpec = 0; iSpec < rtsw.n_mass; iSpec++)
            {
                rtsw.c(iSpec, iGlob, 0) = 0.0;
                rtsw.c(iSpec, iGlob, 1) = 0.0;
            }
        }

        // Set initial conditions independently for each species
        for (int iSpec = 0; iSpec < rtsw.n_mass; iSpec++)
        {
            if (rtsw.RT_Init_Aq_Mode[iSpec] == 0)
            {
                // Mode 0 (constant): Assign a uniform concentration to every cell
                real const_val = rtsw.RT_Init_Aq_Const[iSpec];
                for (iGlob = 0; iGlob < dom.nCellMem; iGlob++)
                {
                    rtsw.c(iSpec, iGlob, 0) = const_val;
                    rtsw.c(iSpec, iGlob, 1) = const_val;
                }
                if (par.masterproc)
                    std::cerr << GOK "SWRTM: Species " << iSpec << ": aq_mode=0(const), set to aq_val=" << const_val << std::endl;
            }
            else if (rtsw.RT_Init_Aq_Mode[iSpec] == 1)
            {
                // Mode 1 (file): Read spatially distributed concentrations from file
                tempStr = rtsw.RT_Init_Aq_File[iSpec];
                if (par.masterproc)
                    std::cerr << GOK "SWRTM: Species " << iSpec << ": aq_mode=1(file), reading from aq_file=" << tempStr << std::endl;
                if (!readRtICFileSW(tempStr, inFolder, rtsw, dom, par, iSpec))
                {
                    std::cerr << RERROR "SWRTM: Failed to read initial condition file for species " << iSpec << std::endl;
                    std::cerr << "      File: " << tempStr << std::endl;
                    exit(-1);
                }
            }
            else
            {
                std::cerr << RERROR "SWRTM: RT_Init_Aq_" << iSpec << ": mode must be 0(const) or 1(file)" << std::endl;
                exit(-1);
            }
        }

        if (par.masterproc)
        {
            std::cerr << GOK "SWRTM: Surface water reactive transport initial condition set" << std::endl;
        }
        return 1;
    }

    /* ====================================================================== */
    /*                    HELPER: toUpperString()                              */
    /* ====================================================================== */

    /**
     * @brief Convert a string to uppercase.
     *
     * Used for case-insensitive comparison of header field names in input
     * files (e.g., "ncols", "NCOLS", "Ncols" are all treated as equivalent).
     *
     * @param s  Input string
     * @return   Uppercase copy of the input string
     */
    std::string toUpperString(const std::string &s)
    {
        std::string result = s;
        std::for_each(result.begin(), result.end(), [](char &c)
                      { c = ::toupper(c); });
        return result;
    }

    /* ====================================================================== */
    /*                    HELPER: trimString()                                 */
    /* ====================================================================== */

    /**
     * @brief Trim leading and trailing whitespace from a string.
     *
     * Removes spaces, tabs, carriage returns, and newlines from both ends
     * of the input string. Returns an empty string if the input consists
     * entirely of whitespace.
     *
     * @param s  Input string (possibly with leading/trailing whitespace)
     * @return   Trimmed string
     */
    std::string trimString(const std::string &s)
    {
        size_t start = s.find_first_not_of(" \t\r\n");
        if (start == std::string::npos)
            return "";
        size_t end = s.find_last_not_of(" \t\r\n");
        return s.substr(start, end - start + 1);
    }

    /* ====================================================================== */
    /*           READ INITIAL CONDITION FILE: readRtICFileSW()                 */
    /* ====================================================================== */

    /**
     * @brief Read spatially distributed initial concentrations from an ASCII grid file.
     *
     * Supports two input file formats:
     *
     * **Simplified format (3-line header):**
     *   <ncols>
     *   <nrows>
     *   <nodata_value>
     *   <data...>  (row-major, nrows x ncols values)
     *
     * **ESRI ASCII Grid format (6-line header):**
     *   NCOLS <ncols>
     *   NROWS <nrows>
     *   XLLCORNER <x>       (ignored, logged)
     *   YLLCORNER <y>       (ignored, logged)
     *   CELLSIZE <size>     (ignored, logged)
     *   NODATA_value <value>
     *   <data...>  (row-major, nrows x ncols values)
     *
     * The file dimensions (NCOLS x NROWS) must exactly match the global domain
     * dimensions (dom.nx_glob x dom.ny_glob). Nodata values in the data are
     * currently passed through (not masked).
     *
     * @param fNameIn   IC filename (without directory prefix)
     * @param fDirIn    Directory prefix for the IC file
     * @param rtsw      Surface water RT state [modified - concentration arrays]
     * @param dom       Computational domain [read-only]
     * @param par       Parallel decomposition info [read-only]
     * @param iSpec     Species index (0-based) to load data for
     *
     * @return 1 on success, 0 on failure
     */
    int readRtICFileSW(std::string fNameIn, std::string fDirIn, RTStateSW &rtsw, Domain &dom, Parallel &par, int iSpec = -1)
    {
        std::string fname = fDirIn + fNameIn;
        std::ifstream fInStream(fname);
        std::string line;
        int tnx = -999, tny = -999, iGlob, iGlobSW, idx, ivg, ii, jj, kk, ii2;
        real tmp, nodata_value = -9999;
        int ndata = dom.ny_glob * dom.nx_glob;
        realArr tmpVar = realArr("var", ndata);
        std::string str;

        // Header parsing state
        bool hasNCOLS = false, hasNROWS = false, hasNODATA = false;
        int headerLinesRead = 0;
        const int MAX_HEADER_LINES = 10; ///< Maximum number of header lines to scan

        if (fInStream.is_open())
        {
            // --- Phase 1: Parse header lines until all required fields are found ---
            while (headerLinesRead < MAX_HEADER_LINES && !fInStream.eof())
            {
                std::getline(fInStream, line);
                headerLinesRead++;

                // Skip blank lines
                if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos)
                    continue;

                // Split line into key and value at the first whitespace
                size_t spacePos = line.find_first_of(" \t");
                std::string key, value;

                if (spacePos != std::string::npos)
                {
                    key = toUpperString(trimString(line.substr(0, spacePos)));
                    value = trimString(line.substr(spacePos + 1));
                }
                else
                {
                    // No whitespace found: treat the entire trimmed line as a value
                    // (simplified format where each header line is just a number)
                    key = "";
                    value = trimString(line);
                }

                // --- Handle ESRI ASCII Grid header fields ---
                if (key == "NCOLS")
                {
                    try
                    {
                        tnx = std::stoi(value);
                        hasNCOLS = true;
                    }
                    catch (...)
                    {
                        if (par.masterproc)
                        {
                            std::cerr << RERROR "SWRTM: Invalid NCOLS value in IC file: " << value << "\n";
                            std::cerr << "      File: " << fname << ", Line: " << headerLinesRead << "\n";
                        }
                        return 0;
                    }
                }
                else if (key == "NROWS")
                {
                    try
                    {
                        tny = std::stoi(value);
                        hasNROWS = true;
                    }
                    catch (...)
                    {
                        if (par.masterproc)
                        {
                            std::cerr << RERROR "SWRTM: Invalid NROWS value in IC file: " << value << "\n";
                            std::cerr << "      File: " << fname << ", Line: " << headerLinesRead << "\n";
                        }
                        return 0;
                    }
                }
                else if (key == "NODATA_VALUE" || key == "NODATAVALUE" || key == "NODATA")
                {
                    try
                    {
                        nodata_value = std::stod(value);
                        hasNODATA = true;
                    }
                    catch (...)
                    {
                        if (par.masterproc)
                        {
                            std::cerr << RERROR "SWRTM: Invalid NODATA_value in IC file: " << value << "\n";
                            std::cerr << "      File: " << fname << ", Line: " << headerLinesRead << "\n";
                        }
                        return 0;
                    }
                }
                else if (key == "XLLCORNER" || key == "YLLCORNER" || key == "CELLSIZE")
                {
                    // ESRI ASCII Grid georeferencing fields -- not used by the model,
                    // but recognized and skipped without error
                    if (par.masterproc)
                    {
                        std::cout << BDASH << "SWRTM: Skipping ESRI ASCII Grid header: " << key << " = " << value << "\n";
                    }
                }
                else if (key.empty())
                {
                    // Simplified format: header lines with no key, just values
                    // Line 1 = ncols, Line 2 = nrows, Line 3 = nodata_value
                    if (!hasNCOLS && tnx == -999)
                    {
                        try
                        {
                            tnx = std::stoi(value);
                            hasNCOLS = true;
                        }
                        catch (...)
                        {
                            if (par.masterproc)
                            {
                                std::cerr << RERROR "SWRTM: Failed to parse NCOLS from line " << headerLinesRead << ": " << value << "\n";
                            }
                            return 0;
                        }
                    }
                    else if (!hasNROWS && tny == -999)
                    {
                        try
                        {
                            tny = std::stoi(value);
                            hasNROWS = true;
                        }
                        catch (...)
                        {
                            if (par.masterproc)
                            {
                                std::cerr << RERROR "SWRTM: Failed to parse NROWS from line " << headerLinesRead << ": " << value << "\n";
                            }
                            return 0;
                        }
                    }
                    else if (!hasNODATA)
                    {
                        try
                        {
                            nodata_value = std::stod(value);
                            hasNODATA = true;
                        }
                        catch (...)
                        {
                            if (par.masterproc)
                            {
                                std::cerr << RERROR "SWRTM: Failed to parse NODATA_value from line " << headerLinesRead << ": " << value << "\n";
                            }
                            return 0;
                        }
                    }
                }
                else
                {
                    // Unrecognized header field -- silently ignored
                    // (Could enable warning output here for debugging)
                }

                // Exit header loop once all three required fields are found
                if (hasNCOLS && hasNROWS && hasNODATA)
                {
                    break;
                }
            }

            // --- Phase 2: Validate required header fields ---
            if (!hasNCOLS)
            {
                if (par.masterproc)
                {
                    std::cerr << RERROR "SWRTM: Missing required header field 'NCOLS' in IC file\n";
                    std::cerr << "      File: " << fname << "\n";
                    std::cerr << "      Expected format:\n";
                    std::cerr << "        NCOLS <ncols>\n";
                    std::cerr << "        NROWS <nrows>\n";
                    std::cerr << "        NODATA_value <value>\n";
                    std::cerr << "        <data...>\n";
                }
                return 0;
            }

            if (!hasNROWS)
            {
                if (par.masterproc)
                {
                    std::cerr << RERROR "SWRTM: Missing required header field 'NROWS' in IC file\n";
                    std::cerr << "      File: " << fname << "\n";
                    std::cerr << "      Expected format:\n";
                    std::cerr << "        NCOLS <ncols>\n";
                    std::cerr << "        NROWS <nrows>\n";
                    std::cerr << "        NODATA_value <value>\n";
                    std::cerr << "        <data...>\n";
                }
                return 0;
            }

            if (!hasNODATA)
            {
                if (par.masterproc)
                {
                    std::cerr << RERROR "SWRTM: Missing required header field 'NODATA_value' in IC file\n";
                    std::cerr << "      File: " << fname << "\n";
                    std::cerr << "      Expected format:\n";
                    std::cerr << "        NCOLS <ncols>\n";
                    std::cerr << "        NROWS <nrows>\n";
                    std::cerr << "        NODATA_value <value>\n";
                    std::cerr << "        <data...>\n";
                }
                return 0;
            }

            // Validate that dimensions are positive
            if (tnx <= 0 || tny <= 0)
            {
                if (par.masterproc)
                {
                    std::cerr << RERROR "SWRTM: Invalid dimensions in IC file: NCOLS=" << tnx << ", NROWS=" << tny << "\n";
                    std::cerr << "      Dimensions must be positive integers.\n";
                    std::cerr << "      File: " << fname << "\n";
                }
                return 0;
            }

            // --- Phase 3: Verify dimensions match the DEM / domain ---
            if (dom.ny_glob != tny || dom.nx_glob != tnx)
            {
                if (par.masterproc)
                {
                    std::cerr << RERROR "SWRTM: RTM SW IC file dimensions don't match DEM parameters.\n";
                    std::cerr << "      IC file: NCOLS=" << tnx << ", NROWS=" << tny << "\n";
                    std::cerr << "      DEM:     nx_glob=" << dom.nx_glob << ", ny_glob=" << dom.ny_glob << "\n";
                    std::cerr << "      File: " << fname << "\n";
                }
                return 0;
            }

            if (par.masterproc)
            {
                std::cerr << GOK << "SWRTM: Successfully read IC file header: NCOLS=" << tnx
                          << ", NROWS=" << tny << ", NODATA_value=" << nodata_value << "\n";
            }

            // --- Phase 4: Read the data grid into a temporary host-side array ---
            for (int ii = 0; ii < ndata; ii++)
            {
                if (!fInStream.fail() && !fInStream.eof())
                {
                    fInStream >> tmp;
                    tmpVar(ii) = tmp;
                }
                else
                {
                    if (par.masterproc)
                    {
                        std::cerr << RERROR "SWRTM: Error reading RTM SW IC file. Not enough data.\n";
                        std::cerr << "      Expected " << ndata << " values (" << tny << " rows x " << tnx << " cols)\n";
                        std::cerr << "      Only read " << ii << " values before EOF.\n";
                        std::cerr << "      File: " << fname << "\n";
                        std::cerr << "      Please check that the data section contains exactly " << ndata << " values.\n";
                    }
                    return 0;
                }
            }
            fInStream.close();
        }
        else if (par.masterproc)
        {
            std::cerr << RERROR "SWRTM: Unable to open IC file: " << fname << "\n";
            std::cerr << "      Please check that the file exists and is readable.\n";
            exit(-1);
        }

        // --- Phase 5: Map global grid data into local domain cells ---
        // Surface water uses nCellMem, which includes halo cells.
        // Only interior cells (dom.nCell) are populated from the global data.
        for (idx = 0; idx < dom.nCell; idx++)
        {
            dom.unpackIndices(idx, jj, ii);

            // Compute the local memory index (including halo offset)
            iGlob = (dom.hc + jj) * (dom.nx + 2 * dom.hc) + (dom.hc + ii);

            // Compute the corresponding global data index for this MPI partition
            ii2 = (par.j_beg + jj) * (dom.nx_glob) + par.i_beg + ii;
            if (ii2 < ndata)
            {
                rtsw.c(iSpec, iGlob, 1) = tmpVar(ii2);
                rtsw.c(iSpec, iGlob, 0) = tmpVar(ii2);
            }
        }
        return 1;
    }

    /* ====================================================================== */
    /*         READ TRANSPORT PARAMETERS: readRTFileSW()                       */
    /* ====================================================================== */

    /**
     * @brief Read surface water reactive transport parameters from rttransportsw.input.
     *
     * Parses a colon-separated key:value file containing transport configuration.
     * Parameters are organized per-species using an "id" key to select the current
     * species, followed by species-specific parameters.
     *
     * **Required global keys:**
     *   - n_mass:           Number of solute species (must be the first species-related key)
     *   - ReactionModule:   Flag to enable/disable reaction module (0 or 1)
     *   - Advection_Scheme: Numerical scheme for advective transport
     *
     * **Per-species keys (preceded by "id: <species_index>"):**
     *   - aq_mode:              Aqueous IC mode (0=constant, 1=file)
     *   - aq_val:               Constant aqueous IC value [mg/L or kg/m^3]
     *   - aq_file:              Filename for spatially distributed aqueous IC
     *   - diffusion_molecular:  Molecular diffusion coefficient [m^2/s]
     *   - alpha_L:              Longitudinal dispersivity [m]
     *   - alpha_T:              Transverse dispersivity [m]
     *   - RainCon:              Rainfall concentration for this species [mg/L or kg/m^3]
     *
     * @param fNameIn  Full path to the rttransportsw.input file
     * @param rtsw     Surface water RT state [modified]
     * @param par      Parallel decomposition info [read-only]
     *
     * @return 1 on success, 0 on failure; exits on missing required parameters
     */
    int readRTFileSW(std::string fNameIn, RTStateSW &rtsw, Parallel &par)
    {
        // Initialize all read-in values to -999 (sentinel for "not set")
        rtsw.RtInitialModeSW = -999;
        rtsw.ReactionModule = -999;
        rtsw.Advection_Scheme = -999;
        rtsw.c0 = -999;
        rtsw.n_mass = -999;
        rtsw.InitialCon = -999;

        rtsw.ReactModel_Type = -999;
        rtsw.Nitrogen_Cycle_Simulation = -999;

        // Clear per-species initialization arrays (resized after n_mass is read)
        rtsw.RT_Init_Aq_Mode.clear();
        rtsw.RT_Init_Aq_Const.clear();
        rtsw.RT_Init_Aq_File.clear();
        rtsw.RainCon.clear();

        // Track which species block is currently being parsed
        int current_species_id = -1;

        // Temporary storage for transport parameters (moved to Kokkos views after allocate())
        rtsw.temp_diffusion_molecular.clear();
        rtsw.temp_alpha_L.clear();
        rtsw.temp_alpha_T.clear();

        std::string strAux;

        // Read the colon-separated key:value file line by line
        std::ifstream fInStream(fNameIn);
        std::string line;
        PsLn pline;
        if (fInStream.is_open())
        {
            while (std::getline(fInStream, line))
            {
                pline.line = line;
                pline.parse();
                if (!pline.key.empty())
                {
                    if (!strcmp("n_mass", pline.key.c_str()))
                    {
                        pline.value >> rtsw.n_mass;
                        std::cout << BDASH "SWRTM: rtsw.n_mass: " << rtsw.n_mass << std::endl;

                        // Resize per-species arrays now that the number of species is known
                        rtsw.RT_Init_Aq_Mode.resize(rtsw.n_mass, -999);
                        rtsw.RT_Init_Aq_Const.resize(rtsw.n_mass, 0.0);
                        rtsw.RT_Init_Aq_File.resize(rtsw.n_mass, "");

                        // Rainfall concentration per species [mg/L or kg/m^3]
                        rtsw.RainCon.resize(rtsw.n_mass, 0.0);

                        // Pre-allocate temporary transport parameter storage
                        rtsw.temp_diffusion_molecular.resize(rtsw.n_mass, 0.0);
                        rtsw.temp_alpha_L.resize(rtsw.n_mass, 0.0);
                        rtsw.temp_alpha_T.resize(rtsw.n_mass, 0.0);
                    }
                    else if (!strcmp("id", pline.key.c_str()))
                    {
                        // Select the current species block for subsequent per-species keys
                        pline.value >> current_species_id;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: Reading species " << current_species_id << " parameters" << std::endl;
                    }
                    else if (!strcmp("aq_mode", pline.key.c_str()))
                    {
                        // Aqueous-phase IC mode: 0 = constant value, 1 = read from file
                        if (current_species_id >= 0 && current_species_id < rtsw.n_mass)
                        {
                            pline.value >> rtsw.RT_Init_Aq_Mode[current_species_id];
                        }
                    }
                    else if (!strcmp("aq_val", pline.key.c_str()))
                    {
                        // Constant aqueous IC concentration [mg/L or kg/m^3]
                        if (current_species_id >= 0 && current_species_id < rtsw.n_mass)
                        {
                            pline.value >> rtsw.RT_Init_Aq_Const[current_species_id];
                            rtsw.RT_Init_Aq_File[current_species_id] = "";
                        }
                    }
                    else if (!strcmp("aq_file", pline.key.c_str()))
                    {
                        // Filename for spatially distributed aqueous IC data
                        if (current_species_id >= 0 && current_species_id < rtsw.n_mass)
                        {
                            pline.value >> rtsw.RT_Init_Aq_File[current_species_id];
                            rtsw.RT_Init_Aq_Const[current_species_id] = 0.0;
                        }
                    }
                    else if (!strcmp("diffusion_molecular", pline.key.c_str()))
                    {
                        // Molecular diffusion coefficient for the current species [m^2/s]
                        if (current_species_id >= 0 && current_species_id < rtsw.n_mass)
                        {
                            pline.value >> rtsw.temp_diffusion_molecular[current_species_id];
                        }
                    }
                    else if (!strcmp("alpha_L", pline.key.c_str()))
                    {
                        // Longitudinal dispersivity for the current species [m]
                        if (current_species_id >= 0 && current_species_id < rtsw.n_mass)
                        {
                            pline.value >> rtsw.temp_alpha_L[current_species_id];
                        }
                    }
                    else if (!strcmp("alpha_T", pline.key.c_str()))
                    {
                        // Transverse dispersivity for the current species [m]
                        if (current_species_id >= 0 && current_species_id < rtsw.n_mass)
                        {
                            pline.value >> rtsw.temp_alpha_T[current_species_id];
                        }
                    }
                    else if (!strcmp("RainCon", pline.key.c_str()))
                    {
                        // Rainfall concentration for the current species [mg/L or kg/m^3]
                        if (current_species_id >= 0 && current_species_id < rtsw.n_mass)
                        {
                            pline.value >> rtsw.RainCon[current_species_id];
                        }
                    }

                    else if (!strcmp("ReactionModule", pline.key.c_str()))
                    {
                        pline.value >> rtsw.ReactionModule;
                    }
                    else if (!strcmp("Advection_Scheme", pline.key.c_str()))
                    {
                        pline.value >> rtsw.Advection_Scheme;
                    }
                }
            }
        }
        else
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "SWRTM: Unable to open " << fNameIn << "\n";
                return 0;
            }
        }

        // --- Validate that all required global parameters were set ---
        if (rtsw.n_mass == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "SWRTM: key n_mass not set." << std::endl;
            exit(-1);
        }

        if (rtsw.ReactionModule == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "SWRTM: key ReactionModule not set." << std::endl;
            exit(-1);
        }
        if (rtsw.Advection_Scheme == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "SWRTM: key Advection_Scheme not set." << std::endl;
            exit(-1);
        }

        // --- Validate that per-species aqueous IC modes are all set ---
        for (int i = 0; i < rtsw.n_mass; i++)
        {
            if (rtsw.RT_Init_Aq_Mode[i] == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: aq_mode for species " << i << " not set." << std::endl;
                exit(-1);
            }
        }

        // --- Print summary of initial conditions (aqueous phase only; no solid phase in SW) ---
        if (par.masterproc)
        {
            std::cerr << GOK "SWRTM: Initial conditions:" << std::endl;
            for (int i = 0; i < rtsw.n_mass; i++)
            {
                std::cerr << "  Species " << i << ": ";
                if (rtsw.RT_Init_Aq_Mode[i] == 0)
                    std::cerr << "aq_mode=0(const), aq_val=" << rtsw.RT_Init_Aq_Const[i];
                else if (rtsw.RT_Init_Aq_Mode[i] == 1)
                    std::cerr << "aq_mode=1(file), aq_file=" << rtsw.RT_Init_Aq_File[i];
                std::cerr << std::endl;
            }
        }

        // --- Print transport parameters (read from temp storage, applied after allocate()) ---
        // At this point the Kokkos views for transport parameters have not yet been
        // allocated. The values stored in temp_* vectors will be copied to the Kokkos
        // views in initialize_rt() after rtsw.allocate() is called.
        if (par.masterproc)
        {
            std::cerr << GOK "SWRTM: Transport properties (read from file):" << std::endl;
            for (int i = 0; i < rtsw.n_mass; i++)
            {
                std::cerr << "  Species " << i << ": diffusion_molecular=" << rtsw.temp_diffusion_molecular[i]
                          << ", alpha_L=" << rtsw.temp_alpha_L[i]
                          << ", alpha_T=" << rtsw.temp_alpha_T[i]
                          << ", RainCon=" << rtsw.RainCon[i] << std::endl;
            }
            std::cerr << GOK "SWRTM: Transport parameters read (will be applied after allocation)\n";
        }
        return 1;
    }

    /* ====================================================================== */
    /*       READ BOUNDARY CONDITIONS: readRTBCFileSW()                        */
    /* ====================================================================== */

    /**
     * @brief Read surface water reactive transport boundary conditions from rtswbc.input.
     *
     * Parses a colon-separated key:value file defining boundary conditions for
     * each external boundary. The file supports per-species boundary types and
     * values, as well as time-varying boundary conditions via external time-series
     * files.
     *
     * **Boundary types (rtbctype / spec_X_bctype):**
     *   - 1 (RTSW_BC_NONE):      No reactive transport BC applied
     *   - 2 (RTSW_BC_CON_CONST): Constant concentration BC [mg/L or kg/m^3]
     *   - 3 (RTSW_BC_ZEROGRAD):  Zero-gradient (Neumann) BC
     *   - 4 (RTSW_BC_CON_T):     Time-varying concentration BC (requires timeseries file)
     *
     * **Per-species keys (within each boundary "id" block):**
     *   - spec_X_bctype:  Boundary type for species X
     *   - spec_X_bcval:   Constant BC value for species X (required for type 2)
     *   - spec_X_bcfile:  Timeseries filename for species X (required for type 4)
     *
     * @param fNameIn  Full path to the rtswbc.input file
     * @param dom      Computational domain [read-only]
     * @param ebc      External boundary conditions [modified]
     * @param par      Parallel decomposition info [read-only]
     * @param state    Simulation state [read-only]
     * @param rtsw     Surface water RT state [read-only]
     *
     * @return 1 on success, 0 on failure
     */
    int readRTBCFileSW(std::string fNameIn, Domain &dom, ExternalBoundaries &ebc, Parallel &par, State &state, RTStateSW &rtsw)
    {
        std::ifstream fInStream(fNameIn);
        std::string dir;
        std::vector<std::string> timeseriesFile;
        std::string line;
        PsLn pline;
        // Extract directory from filename (remove "rtswbc.input" = 12 chars)
        dir = fNameIn.substr(0, fNameIn.length() - 12);

        int ibc = -1;
        int bccount = 0;
        int bccountFound = 0;
        if (fInStream.is_open())
        {
            int bccount = ebc.extbc.size();
            ebc.rtbc.resize(bccount);
            timeseriesFile.resize(bccount);

            // Copy cell index arrays from the existing external BC (flow BC) structure
            for (int k = 0; k < ebc.extbc.size(); k++)
            {
                ebc.rtbc[k].ncellsBC = ebc.extbc[k].ncellsBC;
                ebc.rtbc[k].bcells = intArr("bcells", ebc.rtbc[k].ncellsBC);
                Kokkos::deep_copy(ebc.rtbc[k].bcells, ebc.extbc[k].bcells);
            }

            // Parse the boundary condition file line by line
            while (std::getline(fInStream, line))
            {
                pline.line = line;
                pline.parse();
#if SERGHEI_DEBUG_BOUNDARY
                pline.print();
#endif
                if (!pline.key.empty())
                {
                    if (!strcmp("id", pline.key.c_str()))
                    {
                        // Start a new boundary condition block
                        ibc++;

                        // Initialize per-species data structures for this boundary
                        if (bccount > 0 && ibc >= 0)
                        {
                            ebc.rtbc[ibc].spec_bctype.resize(rtsw.n_mass, -999);

                            // Create independent realArr storage for each species'
                            // spatially distributed boundary values
                            ebc.rtbc[ibc].spec_bcvals.clear();
                            for (int iSpec = 0; iSpec < rtsw.n_mass; iSpec++)
                            {
                                ebc.rtbc[ibc].spec_bcvals.push_back(realArr());
                            }
                            ebc.rtbc[ibc].spec_bcfile.resize(rtsw.n_mass, "");
                            ebc.rtbc[ibc].spec_ts.resize(rtsw.n_mass);
                            ebc.rtbc[ibc].spec_bcval_const.resize(rtsw.n_mass, 0.0);
                            ebc.rtbc[ibc].has_bcfile.resize(rtsw.n_mass, 0);
                        }
                    }
                    else if (!strcmp("rtbctype", pline.key.c_str()) && ibc >= 0)
                    {
                        // Global/default boundary type for this boundary
                        if (bccount > 0)
                        {
                            pline.value >> ebc.rtbc[ibc].bctype;
                            // Validate: must be one of the recognized BC types
                            if (ebc.rtbc[ibc].bctype != RTSW_BC_NONE &&
                                ebc.rtbc[ibc].bctype != RTSW_BC_CON_CONST &&
                                ebc.rtbc[ibc].bctype != RTSW_BC_ZEROGRAD &&
                                ebc.rtbc[ibc].bctype != RTSW_BC_CON_T)
                            {
                                if (par.masterproc)
                                {
                                    std::cerr << RERROR "SWRTM: Invalid rtbctype=" << ebc.rtbc[ibc].bctype
                                              << " for boundary id=" << ibc << ". "
                                              << "Valid values are: 1(NONE), 2(CON_CONST), 3(ZEROGRAD), 4(CON_T)." << std::endl;
                                }
                                return 0;
                            }
                        }
                    }
                    else if (!strcmp("direction", pline.key.c_str()) && ibc >= 0)
                    {
                        // Boundary normal direction vector (x, y components)
                        if (bccount > 0)
                            pline.value >> ebc.rtbc[ibc].normalx >> ebc.rtbc[ibc].normaly;
                    }
                    // --- Per-species boundary condition keys (format: spec_X_<param>) ---
                    else if (pline.key.length() > 5 && pline.key.substr(0, 5) == "spec_" && ibc >= 0)
                    {
                        // Parse species ID and parameter type from the key
                        // Example: "spec_0_bctype" -> species_id=0, param_type="bctype"
                        std::string spec_str = pline.key.substr(5); // Remove "spec_" prefix
                        size_t underscore_pos = spec_str.find('_');
                        if (underscore_pos != std::string::npos)
                        {
                            int spec_id = std::stoi(spec_str.substr(0, underscore_pos));
                            std::string param_type = spec_str.substr(underscore_pos + 1);

                            if (spec_id >= 0 && spec_id < rtsw.n_mass)
                            {
                                if (param_type == "bctype")
                                {
                                    // Boundary type for this species
                                    pline.value >> ebc.rtbc[ibc].spec_bctype[spec_id];
                                    // Validate boundary type
                                    if (ebc.rtbc[ibc].spec_bctype[spec_id] != RTSW_BC_NONE &&
                                        ebc.rtbc[ibc].spec_bctype[spec_id] != RTSW_BC_CON_CONST &&
                                        ebc.rtbc[ibc].spec_bctype[spec_id] != RTSW_BC_ZEROGRAD &&
                                        ebc.rtbc[ibc].spec_bctype[spec_id] != RTSW_BC_CON_T)
                                    {
                                        if (par.masterproc)
                                        {
                                            std::cerr << RERROR "SWRTM: Invalid spec_" << spec_id << "_bctype=" << ebc.rtbc[ibc].spec_bctype[spec_id]
                                                      << " for boundary id=" << ibc << ", species " << spec_id << ". "
                                                      << "Valid values are: 1(NONE), 2(CON_CONST), 3(ZEROGRAD), 4(CON_T)." << std::endl;
                                        }
                                        return 0;
                                    }
                                    if (par.masterproc)
                                        std::cerr << GOK "SWRTM:   Boundary id=" << ibc << ", Species " << spec_id << ", bctype=" << ebc.rtbc[ibc].spec_bctype[spec_id] << std::endl;
                                }
                                else if (param_type == "bcval")
                                {
                                    // Constant BC value for this species [mg/L or kg/m^3]
                                    // Validate that the value is numeric (not a filename)
                                    std::string val_str = pline.value.str();
                                    if (val_str.find('.') != std::string::npos &&
                                        (val_str.find(".input") != std::string::npos))
                                    {
                                        if (par.masterproc)
                                        {
                                            std::cerr << RERROR "SWRTM:   Boundary id=" << ibc << ", Species " << spec_id
                                                      << ": bcval should be a numeric value, not a filename (" << val_str << ")." << std::endl;
                                        }
                                        return 0;
                                    }

                                    pline.value >> ebc.rtbc[ibc].spec_bcval_const[spec_id];
                                    if (par.masterproc)
                                        std::cerr << GOK "SWRTM:   Boundary id=" << ibc << ", Species " << spec_id << ", bcval=" << ebc.rtbc[ibc].spec_bcval_const[spec_id] << std::endl;

                                    // Cross-check: time-varying BCs (type 4) should use bcfile, not bcval
                                    if (ebc.rtbc[ibc].spec_bctype[spec_id] == RTSW_BC_CON_T)
                                    {
                                        if (par.masterproc)
                                        {
                                            std::cerr << RERROR "SWRTM:   Boundary id=" << ibc << ", Species " << spec_id << ": bcval specified for T boundary type (type="
                                                      << ebc.rtbc[ibc].spec_bctype[spec_id] << "). Use bcfile instead." << std::endl;
                                        }
                                        return 0;
                                    }
                                }
                                else if (param_type == "bcfile" || param_type == "file")
                                {
                                    // Timeseries filename for this species (used with time-varying BCs)
                                    pline.value >> ebc.rtbc[ibc].spec_bcfile[spec_id];
                                    ebc.rtbc[ibc].has_bcfile[spec_id] = 1;
                                    if (par.masterproc)
                                        std::cerr << GOK "SWRTM:   Boundary id=" << ibc << ", Species " << spec_id << ", " << param_type << "=" << ebc.rtbc[ibc].spec_bcfile[spec_id] << std::endl;

                                    // Cross-check: constant BCs (type 2) should use bcval, not bcfile
                                    if (ebc.rtbc[ibc].spec_bctype[spec_id] == RTSW_BC_CON_CONST)
                                    {
                                        if (par.masterproc)
                                        {
                                            std::cerr << RERROR "SWRTM:   Boundary id=" << ibc << ", Species " << spec_id << ": " << param_type << " specified for CONST boundary type (type="
                                                      << ebc.rtbc[ibc].spec_bctype[spec_id] << "). Use bcval instead." << std::endl;
                                        }
                                        return 0;
                                    }
                                }
                            }
                        }
                    }
                    else if (ibc >= 0)
                    {
                        if (par.masterproc)
                        {
                            std::cerr << RERROR << "SWRTM: In rtswbc.input: Key " << pline.key << " not understood." << std::endl;
                            return 0;
                        }
                    }
                }
            } // end while
            fInStream.close();
        }
        else
        {
            if (par.masterproc)
                std::cerr << YEXC << "SWRTM: rtswbc.input not found. Default boundaries used." << std::endl;
        }

        // --- Normalize BC normal vectors to unit length ---
        for (int k = 0; k < ebc.rtbc.size(); k++)
        {
            real mod = sqrt(ebc.rtbc[k].normalx * ebc.rtbc[k].normalx + ebc.rtbc[k].normaly * ebc.rtbc[k].normaly);
            ebc.rtbc[k].normalx /= mod;
            ebc.rtbc[k].normaly /= mod;
        }

        // --- Load per-species time-series files for time-varying BCs ---
        // Only species with bctype == RTSW_BC_CON_T (time-varying concentration)
        // require and trigger timeseries file loading.
        for (int k = 0; k < ebc.rtbc.size(); k++)
        {
            for (int iSpec = 0; iSpec < rtsw.n_mass; iSpec++)
            {
                // Determine whether this species needs a timeseries file
                int readts = 0;
                switch (ebc.rtbc[k].spec_bctype[iSpec])
                {
                case RTSW_BC_CON_T:
                    readts = 1;
                    break;
                }

                // Load the timeseries data if a file is specified
                if (readts && ebc.rtbc[k].has_bcfile[iSpec] && ebc.rtbc[k].spec_bcfile[iSpec].length() > 0)
                {
                    std::string fname = dir + ebc.rtbc[k].spec_bcfile[iSpec];
                    int ndatat = 0;

                    // Verify that the timeseries file exists before attempting to read
                    std::ifstream fts_check(fname);
                    if (!fts_check.good())
                    {
                        if (par.masterproc)
                        {
                            std::cerr << RERROR "SWRTM: RTM timeseries file for boundary id=" << k << ", species " << iSpec << " not found: " << fname << std::endl;
                        }
                        return 0;
                    }
                    fts_check.close();

                    // Parse the timeseries file: first line has "ndata <count>",
                    // followed by <count> lines of "time value" pairs
                    std::ifstream fts(fname);
                    if (fts.is_open())
                    {
                        if (par.masterproc)
                            std::cout << GOK << "SWRTM: Reading RTM timeseries file for boundary id=" << k << ", species " << iSpec << ": " << fname << std::endl;
                        fts.ignore(256, ' ');  // Skip the "ndata" label
                        fts >> ndatat;
                        if (ndatat > 0)
                        {
                            ebc.rtbc[k].spec_ts[iSpec].initialise(ndatat);
                        }
                        for (int i = 0; i < ndatat; i++)
                        {
                            if (!fts.fail() && !fts.eof())
                            {
                                fts >> ebc.rtbc[k].spec_ts[iSpec].time(i) >> ebc.rtbc[k].spec_ts[iSpec].value(i);
                            }
                            else
                            {
                                if (par.masterproc)
                                {
                                    std::cerr << RERROR "SWRTM: Error reading timeseries file for boundary id=" << k << ", species " << iSpec << ": " << ebc.rtbc[k].spec_bcfile[iSpec] << std::endl;
                                    return 0;
                                }
                            }
                        }
                        fts.close();
                    }
                    else
                    {
                        if (par.masterproc)
                        {
                            std::cerr << RERROR "SWRTM: Error opening RTM timeseries file for species " << iSpec << ": " << fname << std::endl;
                            return 0;
                        }
                    }
                }
            }
        }
        if (par.masterproc)
            std::cout << GOK << "SWRTM: rtsw SW boundary file parsed and boundaries set" << std::endl;
        return 1;
    }

    /* ====================================================================== */
    /*       READ SOURCE/SINK TERMS: readRTSSFileSW()                          */
    /* ====================================================================== */

    /**
     * @brief Read surface water reactive transport source/sink terms from rtswss.input.
     *
     * Parses a colon-separated key:value file defining spatial source/sink zones
     * for solute injection or extraction. Each source/sink is defined by:
     *   - A polygon file (determines which computational cells are affected)
     *   - Per-species source/sink type and value (constant or time-varying)
     *   - Optional timeseries files for time-varying source/sink rates
     *
     * **File structure:**
     *   sscount: <number_of_sources>
     *   id: <source_id>
     *   rtsstype: <global_sstype>
     *   polygon: <polygon_filename>
     *   spec_X_sstype: <type_for_species_X>
     *   spec_X_ssval: <constant_value_for_species_X> [mg/L/s or kg/m^3/s]
     *   spec_X_ssfile: <timeseries_filename_for_species_X>
     *
     * @param fNameIn  Full path to the rtswss.input file
     * @param dom      Computational domain [read-only]
     * @param ss       Source/sink data container [modified]
     * @param par      Parallel decomposition info [read-only]
     * @param state    Simulation state [read-only]
     * @param rtsw     Surface water RT state [read-only]
     *
     * @return 1 on success, 0 on failure
     */
    int readRTSSFileSW(std::string fNameIn, Domain &dom, SourceSinkData &ss, Parallel &par, State &state, RTStateSW &rtsw)
    {
        std::ifstream fInStream(fNameIn);
        if (!fInStream.is_open())
        {
            if (par.masterproc)
                std::cout << YEXC << "SWRTM: rtswss.input not found. No SW RT external sources/sinks applied." << std::endl;
            return 1;
        }

        // Extract directory from the input file path
        std::string dir = fNameIn.substr(0, fNameIn.find_last_of('/') + 1);
        std::string line;
        PsLn pline;

        int sscount = 0;
        int sscountFound = 0;
        int iss = -1;
        std::vector<std::string> polygonFile;

        // Parse the source/sink file line by line
        while (std::getline(fInStream, line))
        {
            pline.line = line;
            pline.parse();
            if (!pline.key.empty())
            {
                if (!strcmp("sscount", pline.key.c_str()))
                {
                    // Total number of source/sink zones
                    pline.value >> sscount;
                    sscountFound = 1;
                    if (sscount < 1)
                    {
                        if (par.masterproc)
                            std::cout << YEXC << "SWRTM: rtswss.input indicates zero external sources/sinks." << std::endl;
                        return 1;
                    }
                    ss.rtswss.resize(sscount);
                    polygonFile.resize(sscount);
                }
                else if (!strcmp("id", pline.key.c_str()))
                {
                    // Start a new source/sink block
                    iss++;
                    if (sscount > 0 && iss >= 0 && iss < sscount)
                    {
                        pline.value >> ss.rtswss[iss].id;

                        // Initialize per-species data structures for this source/sink
                        ss.rtswss[iss].spec_sstype.resize(rtsw.n_mass, -999);
                        ss.rtswss[iss].spec_ssval_const.resize(rtsw.n_mass, 0.0);
                        ss.rtswss[iss].spec_ssfile.resize(rtsw.n_mass, "");
                        ss.rtswss[iss].has_ssfile.resize(rtsw.n_mass, 0);
                        ss.rtswss[iss].spec_ts.resize(rtsw.n_mass);
                    }
                }
                else if (!strcmp("rtsstype", pline.key.c_str()) && iss >= 0)
                {
                    // Global/default source/sink type for this zone
                    if (sscount > 0 && iss < sscount)
                        pline.value >> ss.rtswss[iss].rtsstype;
                }
                else if (!strcmp("polygon", pline.key.c_str()) && iss >= 0)
                {
                    // Polygon filename defining the spatial extent of this source/sink
                    if (sscount > 0 && iss < sscount)
                    {
                        std::string polyFile = pline.value.str();
                        polyFile.erase(0, polyFile.find_first_not_of(" \t\r\n"));
                        polyFile.erase(polyFile.find_last_not_of(" \t\r\n") + 1);
                        polygonFile[iss] = polyFile;
                    }
                }
                // --- Per-species source/sink keys (format: spec_X_<param>) ---
                else if (pline.key.length() > 5 && pline.key.substr(0, 5) == "spec_" && iss >= 0)
                {
                    if (sscount > 0 && iss < sscount)
                    {
                        std::string spec_str = pline.key.substr(5);
                        size_t underscore_pos = spec_str.find('_');
                        if (underscore_pos != std::string::npos)
                        {
                            int spec_id = std::stoi(spec_str.substr(0, underscore_pos));
                            std::string param_type = spec_str.substr(underscore_pos + 1);

                            if (spec_id >= 0 && spec_id < rtsw.n_mass)
                            {
                                if (param_type == "sstype")
                                {
                                    // Source/sink type for this species
                                    pline.value >> ss.rtswss[iss].spec_sstype[spec_id];
                                }
                                else if (param_type == "ssval")
                                {
                                    // Constant source/sink rate for this species [mg/L/s or kg/m^3/s]
                                    pline.value >> ss.rtswss[iss].spec_ssval_const[spec_id];
                                }
                                else if (param_type == "ssfile" || param_type == "file")
                                {
                                    // Timeseries filename for time-varying source/sink rates
                                    std::string fname = pline.value.str();
                                    fname.erase(0, fname.find_first_not_of(" \t\r\n"));
                                    fname.erase(fname.find_last_not_of(" \t\r\n") + 1);
                                    ss.rtswss[iss].spec_ssfile[spec_id] = fname;
                                    ss.rtswss[iss].has_ssfile[spec_id] = 1;
                                }
                            }
                        }
                    }
                }
                else if (iss >= 0)
                {
                    if (par.masterproc)
                    {
                        std::cerr << RERROR << "SWRTM: In rtswss.input: Key " << pline.key << " not understood." << std::endl;
                        return 0;
                    }
                }
            }
        }
        fInStream.close();

        // --- Validate that sscount was defined ---
        if (!sscountFound)
        {
            if (par.masterproc)
                std::cerr << RERROR << "SWRTM: Number of source/sinks not defined in rtswss.input. Please define 'sscount'" << std::endl;
            return 0;
        }

        // --- Validate that all declared source/sink blocks were present ---
        if (iss + 1 < sscount)
        {
            if (par.masterproc)
            {
                std::cerr << RERROR << "SWRTM: Expected " << sscount << " source/sink input blocks, but only found " << iss + 1 << " in rtswss.input." << std::endl;
                return 0;
            }
        }

        // --- Read polygon files and determine affected computational cells ---
        // Each polygon file defines the spatial footprint of a source/sink zone.
        // The polygon is used to select which computational cells fall within the zone.
        for (int k = 0; k < sscount; k++)
        {
            if (polygonFile[k].empty())
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: Polygon file missing for source/sink id '" << ss.rtswss[k].id << "'" << std::endl;
                return 0;
            }

            std::string polyPath = dir + polygonFile[k];
            std::ifstream fPoly(polyPath);
            if (!fPoly.is_open())
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: Error opening polygon file: " << polyPath << std::endl;
                return 0;
            }

            // Read number of polygon vertices
            int npts = 0;
            std::string dummy;
            fPoly >> dummy >> npts;
            if (npts <= 0)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: Invalid number of points in polygon file: " << polyPath << std::endl;
                return 0;
            }

            // Read polygon vertex coordinates [m]
            realArr xPoly("xPoly", npts);
            realArr yPoly("yPoly", npts);
            for (int i = 0; i < npts; i++)
            {
                if (!fPoly.fail() && !fPoly.eof())
                {
                    fPoly >> xPoly(i) >> yPoly(i);
                }
                else
                {
                    if (par.masterproc)
                    {
                        std::cerr << RERROR "SWRTM: Error reading polygon file for source/sink id '" << ss.rtswss[k].id << "': " << polyPath << std::endl;
                        return 0;
                    }
                }
            }
            fPoly.close();

            // Use the polygon to identify which computational cells are inside the zone
            if (!ss.rtswss[k].find_icells_sw(state, ss.rtswss[k].id, dom, par, npts, xPoly, yPoly))
            {
                return 0;
            }
        }

        // --- Load timeseries files for time-varying source/sink rates ---
        // For each species that has an associated timeseries file, read the
        // time-value pairs into the spec_ts structure.
        for (int k = 0; k < sscount; k++)
        {
            for (int iSpec = 0; iSpec < rtsw.n_mass; iSpec++)
            {
                if (ss.rtswss[k].has_ssfile[iSpec] && ss.rtswss[k].spec_ssfile[iSpec].length() > 0)
                {
                    std::string fname = dir + ss.rtswss[k].spec_ssfile[iSpec];
                    std::ifstream fts(fname);
                    if (fts.is_open())
                    {
                        int ndatat = 0;
                        std::string dummy;
                        fts >> dummy >> ndatat; // Parse "ndata <count>" header
                        if (ndatat > 0)
                        {
                            ss.rtswss[k].spec_ts[iSpec].initialise(ndatat);
                            for (int i = 0; i < ndatat; i++)
                            {
                                if (!fts.fail() && !fts.eof())
                                {
                                    fts >> ss.rtswss[k].spec_ts[iSpec].time(i) >> ss.rtswss[k].spec_ts[iSpec].value(i);
                                }
                                else
                                {
                                    if (par.masterproc)
                                    {
                                        std::cerr << RERROR "SWRTM: Error reading RT SourceSink timeseries file for species " << iSpec << ": " << fname << std::endl;
                                        return 0;
                                    }
                                }
                            }
                        }
                        fts.close();
                    }
                    else
                    {
                        if (par.masterproc)
                            std::cerr << RERROR "SWRTM: Error opening RT SourceSink timeseries file: " << fname << std::endl;
                        return 0;
                    }
                }
            }
        }

        if (par.masterproc)
            std::cout << GOK << "SWRTM: rtswss.input parsed and SW RT source/sinks set" << std::endl;
        return 1;
    }

    /* ====================================================================== */
    /*       READ REACTION PARAMETERS: readRTReactionFile()                    */
    /* ====================================================================== */

    /**
     * @brief Read surface water reaction parameters from reactionsw.input.
     *
     * Parses a colon-separated key:value file containing all parameters needed
     * by the surface water reaction module. The reaction module is activated
     * when ReactionModule == 1 in the transport parameters file.
     *
     * **Reaction parameter categories:**
     *
     * 1. Global / environmental parameters:
     *    - ReactModel_Type:           Reaction model selector (1 = nitrogen cycle)
     *    - Nitrogen_Cycle_Simulation: Flag for nitrogen cycle simulation
     *    - Temp_Coeff_Theta:          Arrhenius temperature coefficient [dimensionless]
     *    - Opt_Temp:                  Optimal temperature for reactions [deg C]
     *    - Opt_pH:                    Optimal pH for reactions [dimensionless]
     *
     * 2. Decay reaction (per-species, multi-value):
     *    - Lambda_1: First-order decay rate constants [1/s] (one per species)
     *    - Lambda_2: Second-order decay rate constants [1/s] (one per species)
     *
     * 3. Nitrification (NH4+ -> NO3-):
     *    - Rate_Max_Nit:  Maximum nitrification rate [mg/L/s]
     *    - K_Monod_NH4:   Monod half-saturation constant for NH4+ [mg/L]
     *    - K_Monod_DO:    Monod half-saturation constant for dissolved oxygen [mg/L]
     *
     * 4. Denitrification (NO3- -> N2):
     *    - Rate_Max_Denit: Maximum denitrification rate [mg/L/s]
     *    - K_Monod_NO3:    Monod half-saturation constant for NO3- [mg/L]
     *    - K_Monod_DOC:    Monod half-saturation constant for dissolved organic carbon [mg/L]
     *    - Ki_Inhib_DO:    Inhibition constant for dissolved oxygen [mg/L]
     *
     * 5. Aerobic respiration (heterotrophic):
     *    - Rate_Max_Hetero:      Maximum heterotrophic respiration rate [mg/L/s]
     *    - K_Monod_DOC_Aerobic:  Monod half-saturation for DOC [mg/L]
     *    - K_Monod_DO_Aerobic:   Monod half-saturation for DO [mg/L]
     *
     * 6. Mineralization (organic N -> NH4+):
     *    - Rate_Max_Min:   Maximum mineralization rate [mg/L/s]
     *    - K_Monod_DON:    Monod half-saturation for dissolved organic nitrogen [mg/L]
     *    - K_Monod_DO_Min: Monod half-saturation for DO [mg/L]
     *
     * 7. Surface-water-specific source/sink:
     *    - Ka_sw:      Atmospheric reaeration coefficient [1/s]
     *    - DO_sat:     Dissolved oxygen saturation concentration [mg/L]
     *    - DOC_eq_sw:  Equilibrium dissolved organic carbon concentration [mg/L]
     *    - DON_eq_sw:  Equilibrium dissolved organic nitrogen concentration [mg/L]
     *    - K_rel_sw:   Release rate constant for sediment-water exchange [1/s]
     *
     * @note Modeled after RTInitGW::readRTReactionFile, adapted for surface water.
     *
     * @param fNameIn  Full path to the reactionsw.input file
     * @param rtsw     Surface water RT state [modified]
     * @param par      Parallel decomposition info [read-only]
     *
     * @return 1 on success, 0 on failure; exits on missing required parameters
     */
    int readRTReactionFile(std::string fNameIn, RTStateSW &rtsw, Parallel &par)
    {
        // Initialize all scalar parameters to -999 (sentinel for "not yet set")
        rtsw.ReactModel_Type = -999;
        rtsw.Nitrogen_Cycle_Simulation = -999;

        // Decay reaction parameters (per-species arrays)
        rtsw.Lambda_1 = realArr("Lambda_1", rtsw.n_mass);  ///< First-order decay rate [1/s]
        rtsw.Lambda_2 = realArr("Lambda_2", rtsw.n_mass);  ///< Second-order decay rate [1/s]

        // Nitrification parameters (single values)
        rtsw.Rate_Max_Nit = -999;   ///< Maximum nitrification rate [mg/L/s]
        rtsw.K_Monod_NH4 = -999;    ///< Monod half-saturation for NH4+ [mg/L]
        rtsw.K_Monod_DO = -999;     ///< Monod half-saturation for DO [mg/L]

        // Denitrification parameters (single values)
        rtsw.Rate_Max_Denit = -999;  ///< Maximum denitrification rate [mg/L/s]
        rtsw.K_Monod_NO3 = -999;     ///< Monod half-saturation for NO3- [mg/L]
        rtsw.K_Monod_DOC = -999;     ///< Monod half-saturation for DOC [mg/L]
        rtsw.Ki_Inhib_DO = -999;     ///< DO inhibition constant [mg/L]

        // Aerobic respiration parameters (single values)
        rtsw.Rate_Max_Hetero = -999;       ///< Maximum heterotrophic respiration rate [mg/L/s]
        rtsw.K_Monod_DOC_Aerobic = -999;   ///< Monod half-saturation for DOC [mg/L]
        rtsw.K_Monod_DO_Aerobic = -999;    ///< Monod half-saturation for DO [mg/L]

        // Mineralization parameters (single values)
        rtsw.Rate_Max_Min = -999;    ///< Maximum mineralization rate [mg/L/s]
        rtsw.K_Monod_DON = -999;     ///< Monod half-saturation for DON [mg/L]
        rtsw.K_Monod_DO_Min = -999;  ///< Monod half-saturation for DO [mg/L]

        // Surface-water-specific source/sink parameters
        rtsw.Ka_sw = -999;       ///< Atmospheric reaeration coefficient [1/s]
        rtsw.DO_sat = -999;      ///< Dissolved oxygen saturation concentration [mg/L]
        rtsw.DOC_eq_sw = -999;   ///< Equilibrium DOC concentration [mg/L]
        rtsw.DON_eq_sw = -999;   ///< Equilibrium DON concentration [mg/L]
        rtsw.K_rel_sw = -999;    ///< Sediment-water exchange release rate [1/s]

        // Environmental correction factors
        rtsw.Temp_Coeff_Theta = -999;  ///< Arrhenius temperature coefficient [dimensionless]
        rtsw.Opt_Temp = -999;          ///< Optimal temperature [deg C]
        rtsw.Opt_pH = -999;            ///< Optimal pH [dimensionless]

        std::string strAux;

        // Read the colon-separated key:value file line by line
        std::ifstream fInStream(fNameIn);
        std::string line;
        PsLn pline;
        if (fInStream.is_open())
        {
            while (std::getline(fInStream, line))
            {
                pline.line = line;
                pline.parse();
                if (!pline.key.empty())
                {
                    // ---- Global / environmental parameters ----
                    if (!strcmp("ReactModel_Type", pline.key.c_str()))
                    {
                        pline.value >> rtsw.ReactModel_Type;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: ReactModel_Type: " << rtsw.ReactModel_Type << std::endl;
                    }
                    else if (!strcmp("Nitrogen_Cycle_Simulation", pline.key.c_str()))
                    {
                        pline.value >> rtsw.Nitrogen_Cycle_Simulation;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: Nitrogen_Cycle_Simulation: " << rtsw.Nitrogen_Cycle_Simulation << std::endl;
                    }
                    else if (!strcmp("SPECIES_ID", pline.key.c_str()))
                    {
                        // Species identifiers (informational, used for debugging)
                        std::vector<int> species_ids;
                        int id;
                        while (pline.value >> id)
                        {
                            species_ids.push_back(id);
                        }
                        if (par.masterproc)
                        {
                            std::cerr << GOK "SWRTM: SPECIES_ID: ";
                            for (size_t i = 0; i < species_ids.size(); i++)
                            {
                                std::cerr << species_ids[i];
                                if (i < species_ids.size() - 1)
                                    std::cerr << "; ";
                            }
                            std::cerr << std::endl;
                        }
                    }
                    else if (!strcmp("SPECIES_Name", pline.key.c_str()))
                    {
                        // Species names (informational, used for debugging)
                        std::vector<std::string> species_names;
                        std::string name;
                        while (pline.value >> name)
                        {
                            species_names.push_back(name);
                        }
                        if (par.masterproc)
                        {
                            std::cerr << GOK "SWRTM: SPECIES_Name: ";
                            for (size_t i = 0; i < species_names.size(); i++)
                            {
                                std::cerr << species_names[i];
                                if (i < species_names.size() - 1)
                                    std::cerr << "; ";
                            }
                            std::cerr << std::endl;
                        }
                    }

                    // ---- Decay reaction parameters (per-species) ----
                    else if (!strcmp("Lambda_1", pline.key.c_str()))
                    {
                        // First-order decay rate constant for each species [1/s]
                        for (int i = 0; i < rtsw.n_mass; i++)
                        {
                            pline.value >> rtsw.Lambda_1(i);
                        }
                        if (par.masterproc)
                        {
                            std::cerr << GOK "SWRTM: Lambda_1: ";
                            for (int i = 0; i < rtsw.n_mass; i++)
                            {
                                std::cerr << rtsw.Lambda_1(i);
                                if (i < rtsw.n_mass - 1)
                                    std::cerr << "; ";
                            }
                            std::cerr << std::endl;
                        }
                    }
                    else if (!strcmp("Lambda_2", pline.key.c_str()))
                    {
                        // Second-order decay rate constant for each species [1/s]
                        for (int i = 0; i < rtsw.n_mass; i++)
                        {
                            pline.value >> rtsw.Lambda_2(i);
                        }
                        if (par.masterproc)
                        {
                            std::cerr << GOK "SWRTM: Lambda_2: ";
                            for (int i = 0; i < rtsw.n_mass; i++)
                            {
                                std::cerr << rtsw.Lambda_2(i);
                                if (i < rtsw.n_mass - 1)
                                    std::cerr << "; ";
                            }
                            std::cerr << std::endl;
                        }
                    }

                    // ---- Nitrification parameters (single values) ----
                    else if (!strcmp("Rate_Max_Nit", pline.key.c_str()))
                    {
                        pline.value >> rtsw.Rate_Max_Nit;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: Rate_Max_Nit: " << rtsw.Rate_Max_Nit << std::endl;
                    }
                    else if (!strcmp("K_Monod_NH4", pline.key.c_str()))
                    {
                        pline.value >> rtsw.K_Monod_NH4;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: K_Monod_NH4: " << rtsw.K_Monod_NH4 << std::endl;
                    }
                    else if (!strcmp("K_Monod_DO", pline.key.c_str()))
                    {
                        pline.value >> rtsw.K_Monod_DO;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: K_Monod_DO: " << rtsw.K_Monod_DO << std::endl;
                    }

                    // ---- Denitrification parameters (single values) ----
                    else if (!strcmp("Rate_Max_Denit", pline.key.c_str()))
                    {
                        pline.value >> rtsw.Rate_Max_Denit;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: Rate_Max_Denit: " << rtsw.Rate_Max_Denit << std::endl;
                    }
                    else if (!strcmp("K_Monod_NO3", pline.key.c_str()))
                    {
                        pline.value >> rtsw.K_Monod_NO3;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: K_Monod_NO3: " << rtsw.K_Monod_NO3 << std::endl;
                    }
                    else if (!strcmp("K_Monod_DOC", pline.key.c_str()))
                    {
                        pline.value >> rtsw.K_Monod_DOC;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: K_Monod_DOC: " << rtsw.K_Monod_DOC << std::endl;
                    }
                    else if (!strcmp("Ki_Inhib_DO", pline.key.c_str()))
                    {
                        pline.value >> rtsw.Ki_Inhib_DO;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: Ki_Inhib_DO: " << rtsw.Ki_Inhib_DO << std::endl;
                    }

                    // ---- Aerobic respiration parameters (single values) ----
                    else if (!strcmp("Rate_Max_Hetero", pline.key.c_str()))
                    {
                        pline.value >> rtsw.Rate_Max_Hetero;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: Rate_Max_Hetero: " << rtsw.Rate_Max_Hetero << std::endl;
                    }
                    else if (!strcmp("K_Monod_DOC_Aerobic", pline.key.c_str()))
                    {
                        pline.value >> rtsw.K_Monod_DOC_Aerobic;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: K_Monod_DOC_Aerobic: " << rtsw.K_Monod_DOC_Aerobic << std::endl;
                    }
                    else if (!strcmp("K_Monod_DO_Aerobic", pline.key.c_str()))
                    {
                        pline.value >> rtsw.K_Monod_DO_Aerobic;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: K_Monod_DO_Aerobic: " << rtsw.K_Monod_DO_Aerobic << std::endl;
                    }

                    // ---- Mineralization parameters (single values) ----
                    else if (!strcmp("Rate_Max_Min", pline.key.c_str()))
                    {
                        pline.value >> rtsw.Rate_Max_Min;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: Rate_Max_Min: " << rtsw.Rate_Max_Min << std::endl;
                    }
                    else if (!strcmp("K_Monod_DON", pline.key.c_str()))
                    {
                        pline.value >> rtsw.K_Monod_DON;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: K_Monod_DON: " << rtsw.K_Monod_DON << std::endl;
                    }
                    else if (!strcmp("K_Monod_DO_Min", pline.key.c_str()))
                    {
                        pline.value >> rtsw.K_Monod_DO_Min;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: K_Monod_DO_Min: " << rtsw.K_Monod_DO_Min << std::endl;
                    }

                    // ---- Surface-water-specific source/sink parameters ----
                    else if (!strcmp("Ka_sw", pline.key.c_str()))
                    {
                        pline.value >> rtsw.Ka_sw;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: Ka_sw: " << rtsw.Ka_sw << std::endl;
                    }
                    else if (!strcmp("DO_sat", pline.key.c_str()))
                    {
                        pline.value >> rtsw.DO_sat;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: DO_sat: " << rtsw.DO_sat << std::endl;
                    }
                    else if (!strcmp("DOC_eq_sw", pline.key.c_str()))
                    {
                        pline.value >> rtsw.DOC_eq_sw;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: DOC_eq_sw: " << rtsw.DOC_eq_sw << std::endl;
                    }
                    else if (!strcmp("DON_eq_sw", pline.key.c_str()))
                    {
                        pline.value >> rtsw.DON_eq_sw;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: DON_eq_sw: " << rtsw.DON_eq_sw << std::endl;
                    }
                    else if (!strcmp("K_rel_sw", pline.key.c_str()))
                    {
                        pline.value >> rtsw.K_rel_sw;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: K_rel_sw: " << rtsw.K_rel_sw << std::endl;
                    }

                    // ---- Environmental correction factors ----
                    else if (!strcmp("Temp_Coeff_Theta", pline.key.c_str()))
                    {
                        pline.value >> rtsw.Temp_Coeff_Theta;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: Temp_Coeff_Theta: " << rtsw.Temp_Coeff_Theta << std::endl;
                    }
                    else if (!strcmp("Opt_Temp", pline.key.c_str()))
                    {
                        pline.value >> rtsw.Opt_Temp;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: Opt_Temp: " << rtsw.Opt_Temp << std::endl;
                    }
                    else if (!strcmp("Opt_pH", pline.key.c_str()))
                    {
                        pline.value >> rtsw.Opt_pH;
                        if (par.masterproc)
                            std::cerr << GOK "SWRTM: Opt_pH: " << rtsw.Opt_pH << std::endl;
                    }
                }
            }
        }
        else
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "SWRTM: Unable to open " << fNameIn << "\n";
                return 0;
            }
        }

        // ===== Mandatory parameter validation =====
        if (rtsw.ReactModel_Type == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "SWRTM: key ReactModel_Type not set." << std::endl;
            exit(-1);
        }
        if (rtsw.Nitrogen_Cycle_Simulation == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "SWRTM: key Nitrogen_Cycle_Simulation not set." << std::endl;
            exit(-1);
        }

        // ===== Additional parameter validation when ReactModel_Type == 1 =====
        if (rtsw.ReactModel_Type == 1)
        {
            // Validate nitrification parameters
            if (rtsw.Rate_Max_Nit == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key Rate_Max_Nit not set." << std::endl;
                exit(-1);
            }
            if (rtsw.K_Monod_NH4 == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key K_Monod_NH4 not set." << std::endl;
                exit(-1);
            }
            if (rtsw.K_Monod_DO == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key K_Monod_DO not set." << std::endl;
                exit(-1);
            }

            // Validate denitrification parameters
            if (rtsw.Rate_Max_Denit == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key Rate_Max_Denit not set." << std::endl;
                exit(-1);
            }
            if (rtsw.K_Monod_NO3 == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key K_Monod_NO3 not set." << std::endl;
                exit(-1);
            }
            if (rtsw.K_Monod_DOC == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key K_Monod_DOC not set." << std::endl;
                exit(-1);
            }
            if (rtsw.Ki_Inhib_DO == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key Ki_Inhib_DO not set." << std::endl;
                exit(-1);
            }

            // Validate aerobic respiration parameters
            if (rtsw.Rate_Max_Hetero == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key Rate_Max_Hetero not set." << std::endl;
                exit(-1);
            }
            if (rtsw.K_Monod_DOC_Aerobic == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key K_Monod_DOC_Aerobic not set." << std::endl;
                exit(-1);
            }
            if (rtsw.K_Monod_DO_Aerobic == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key K_Monod_DO_Aerobic not set." << std::endl;
                exit(-1);
            }

            // Validate mineralization parameters
            if (rtsw.Rate_Max_Min == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key Rate_Max_Min not set." << std::endl;
                exit(-1);
            }
            if (rtsw.K_Monod_DON == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key K_Monod_DON not set." << std::endl;
                exit(-1);
            }
            if (rtsw.K_Monod_DO_Min == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key K_Monod_DO_Min not set." << std::endl;
                exit(-1);
            }

            // Validate surface-water-specific source/sink parameters
            if (rtsw.Ka_sw == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key Ka_sw not set." << std::endl;
                exit(-1);
            }
            if (rtsw.DO_sat == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key DO_sat not set." << std::endl;
                exit(-1);
            }
            if (rtsw.DOC_eq_sw == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key DOC_eq_sw not set." << std::endl;
                exit(-1);
            }
            if (rtsw.DON_eq_sw == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key DON_eq_sw not set." << std::endl;
                exit(-1);
            }
            if (rtsw.K_rel_sw == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key K_rel_sw not set." << std::endl;
                exit(-1);
            }

            // Validate environmental correction factors
            if (rtsw.Temp_Coeff_Theta == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key Temp_Coeff_Theta not set." << std::endl;
                exit(-1);
            }
            if (rtsw.Opt_Temp == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key Opt_Temp not set." << std::endl;
                exit(-1);
            }
            if (rtsw.Opt_pH == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "SWRTM: key Opt_pH not set." << std::endl;
                exit(-1);
            }
        }

        // --- Print summary of all reaction parameters ---
        if (par.masterproc)
        {
            std::cerr << GOK "SWRTM: RTM Reaction parameters read\n";
            std::cout << "SWRTM: DEBUG: === Reaction Parameters Check ===" << std::endl;
            std::cout << "SWRTM: DEBUG: ReactModel_Type = " << rtsw.ReactModel_Type << std::endl;
            std::cout << "SWRTM: DEBUG: Nitrogen_Cycle_Simulation = " << rtsw.Nitrogen_Cycle_Simulation << std::endl;
            std::cout << "SWRTM: DEBUG: n_mass = " << rtsw.n_mass << std::endl;
            std::cout << "SWRTM: DEBUG: Rate_Max_Nit = " << rtsw.Rate_Max_Nit << std::endl;
            std::cout << "SWRTM: DEBUG: K_Monod_NH4 = " << rtsw.K_Monod_NH4 << std::endl;
            std::cout << "SWRTM: DEBUG: Rate_Max_Min = " << rtsw.Rate_Max_Min << std::endl;
            std::cout << "SWRTM: DEBUG: Ka_sw = " << rtsw.Ka_sw << std::endl;
            std::cout << "SWRTM: DEBUG: DO_sat = " << rtsw.DO_sat << std::endl;
            std::cout << "SWRTM: DEBUG: ===============================" << std::endl;
        }

        return 1;
    }
};
#endif
