/* -*- mode: c++; c-default-style: "linux" -*- */
/*!
 * @file RTInitGW.h
 * @brief Initialization routines for groundwater reactive transport modeling (GWRTM).
 *
 * This header defines the RTInitgw class, which is responsible for reading
 * all input files and setting up the initial state of the subsurface reactive
 * transport model within the SERGHEI framework.  The class extends Initializer
 * and is conditionally compiled when SERGHEI_SUBSURFACE_TRANSPORT is enabled.
 *
 * The initialization pipeline proceeds in the following order:
 *   1. Read transport parameters       (rttransportgw.input)
 *   2. Allocate state arrays and MPI buffers
 *   3. Read boundary conditions         (rtgwbc.input)
 *   4. Read source/sink terms           (rtgwss.input)
 *   5. Read reaction parameters         (reactiongw.input)  -- if ReactionModule == 1
 *   6. Set initial solute concentrations (aquifer + solid phase)
 *   7. Perform MPI halo exchange and apply boundary conditions
 *   8. Write initial output
 *
 * Supported features:
 *   - Multi-species transport (arbitrary number of dissolved and adsorbed species)
 *   - Per-species initial conditions via constant value or spatial file
 *   - Linear / Freundlich sorption isotherms
 *   - Nitrogen-cycle reaction module (nitrification, denitrification,
 *     aerobic respiration, mineralization)
 *   - Environmental correction factors (temperature, pH)
 *
 * @note All code logic, variable names, and function signatures remain unchanged.
 *       Only comments and documentation annotations have been added or translated.
 */

#ifndef _RT_INITGW_H_
#define _RT_INITGW_H_

#include "define.h"
#include "Domain.h"
#include "FileIO.h"
#include "GwDomain.h"
#include "GwState.h"
#include "Parallel.h"
#include "Parser.h"
#include "SourceSink.h"
#include "RTBCGW.h"
#include "RTStateGW.h"

#if SERGHEI_SUBSURFACE_TRANSPORT

// ===========================================================================
// Class: RTInitgw
// ===========================================================================
/*!
 * @class RTInitgw
 * @brief Initializes the groundwater reactive transport system.
 *
 * Extends the base Initializer class.  Contains a lightweight line parser
 * (PsLn) and the top-level entry point initialize_rt(), as well as helper
 * methods for reading each input file and setting initial concentrations.
 */
class RTInitgw : public Initializer
{

    // =======================================================================
    // Inner helper class: PsLn (Parse Line)
    // =======================================================================
    /*!
     * @class PsLn
     * @brief A lightweight colon-separated key-value line parser.
     *
     * Parses input files of the form:
     *   key : value  // optional trailing comment
     *
     * Semicolons in the value field are replaced with spaces to facilitate
     * multi-value parsing (e.g., lists of numbers separated by ';').
     * This class mirrors the functionality of Parser::ParserLine.
     */
    class PsLn
    {
    public:
        std::string line;      //!< Raw input line
        std::string key;       //!< Extracted key (left of the colon)
        std::stringstream value; //!< Extracted value (right of the colon)

        //! Convert the entire line to lowercase (useful for case-insensitive matching)
        void lowercase()
        {
            std::for_each(line.begin(), line.end(), [](char &c)
                          { c = ::tolower(c); });
        }

        //! Parse the stored line into key and value fields.
        /*!
         * Lines beginning with "//" are treated as comments and skipped.
         * Leading/trailing whitespace and tabs in the key are stripped.
         * Any trailing comment after "//" in the value portion is removed.
         * Semicolons in the value are replaced with spaces for stream extraction.
         */
        void parse()
        {
            // Clear previous state in case the object is reused
            key.clear();
            value.clear();
            // If the line isn't empty and doesn't begin with a comment specifier, split it based on the colon
            if (!line.empty() && line.find("//", 0) != 0)
            {
                uint splitloc = line.find(':', 0);
                key = line.substr(0, splitloc);
                key.erase(std::remove(key.begin(), key.end(), ' '), key.end());
                key.erase(std::remove(key.begin(), key.end(), '\t'), key.end());
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
                // Replace semicolons with spaces to allow proper parsing
                std::replace(strloc.begin(), strloc.end(), ';', ' ');
                value.str(strloc);
            }
        }

        //! Print the parsed line, key, and value to stdout (for debugging)
        void print() { std::cout << "line: " << line << "\tkey: " << key << "\tvalue: " << value.str() << std::endl; }
    };

public:
    // =======================================================================
    // initialize_rt()
    // =======================================================================
    /*!
     * @brief Main entry point for groundwater reactive transport initialization.
     *
     * Orchestrates the full initialization sequence:
     *   1. Read transport parameters from rttransportgw.input
     *   2. Allocate RT state arrays and MPI communication buffers
     *   3. Read boundary conditions from rtgwbc.input
     *   4. Read source/sink terms from rtgwss.input
     *   5. Read reaction parameters from reactiongw.input (if ReactionModule == 1)
     *   6. Set initial solute concentrations (aquifer liquid + solid phase)
     *   7. Copy concentrations from timeIdx=1 to timeIdx=0
     *   8. Perform MPI halo exchange for both time levels
     *   9. Apply boundary conditions for all species
     *  10. Write initial output files
     *
     * @param[in]  rt         Reactive transport state (concentrations, parameters)
     * @param[in]  gw         Groundwater state (hydraulic head, etc.)
     * @param[in]  gdom       Groundwater domain (grid, geometry)
     * @param[in]  rtgbc      Reactive transport boundary condition container
     * @param[in]  gmpi       Groundwater MPI communication handler
     * @param[in]  par        Parallel decomposition information
     * @param[in]  io         File I/O handler
     * @param[in]  ss         Source/sink term container
     * @param[in]  inFolder   Input directory path (trailing slash expected)
     * @param[in]  outFolder  Output directory path (trailing slash expected)
     * @return 1 on success, 0 on failure
     */
    int initialize_rt(RTStateGW &rt, GwState &gw, GwDomain &gdom, RTSubsurfaceBoundaries &rtgbc, GwMPI &gmpi, Parallel &par, FileIO &io, SourceSink &ss, std::string inFolder, std::string outFolder)
    {
        int flag = -1;
        int ii, jj, kk, idx, iGlob, iGlobSW;

        // --- Step 1: Read solute transport input file ---
        std::string fNameIn = inFolder + "rttransportgw.input";
        if (!readRTFile(fNameIn, rt, par))
        {
            std::cerr << GOK << " GWRTM: Reading in reactive transport input file failed." << std::endl;
            return 0;
        }

        // --- Step 2: Allocate subsurface state variable arrays and MPI buffers ---
        rt.allocate(gdom);
        gmpi.allocate_rtgw(gdom);

        // --- Step 3: Read boundary conditions for reactive transport ---
        fNameIn = inFolder + "rtgwbc.input";
        if (!readRTBCFile(fNameIn, gdom, rtgbc, par, rt, gw))
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "GWRTM: Unable to read transport BC file" << std::endl;
                return 0;
            }
        }

        // --- Step 4: Read groundwater RT source/sink terms ---
        fNameIn = inFolder + "rtgwss.input";
        if (!readRTSSFile(fNameIn, gdom, par, rt, gw, ss))
        {
            if (par.masterproc)
            {
                std::cerr << YEXC << " GWRTM: transport source/sink file (rtgwss.input) not found or failed to read. Assuming no custom RT source sinks." << std::endl;
            }
        }

        // --- Step 5: Read solute physico-chemical reaction parameters ---
        if (rt.ReactionModule == 1)
        {
            fNameIn = inFolder + "reactiongw.input";
            if (!readRTReactionFile(fNameIn, rt, par))
            {
                if (par.masterproc)
                {
                    std::cerr << RERROR "GWRTM: Unable to read transport Reaction file" << std::endl;
                    return 0;
                }
            }
        }

        // --- Step 6: Set initial conditions for reactive transport ---
        fNameIn = inFolder;
        if (!setRtState(fNameIn, rt, gw, gdom, par, io))
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "GWRTM: Unable to read RTM transport IC file" << std::endl;
                return 0;
            }
        }

        // --- Step 7: Copy concentrations from the new time level (timeIdx=1)
        //              to the old time level (timeIdx=0) for all species ---
        Kokkos::parallel_for(
            gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
                for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
                {
                    rt.c(iSpec, iGlob, 0) = rt.c(iSpec, iGlob, 1);
                    rt.c_solid(iSpec, iGlob, 0) = rt.c_solid(iSpec, iGlob, 1);
                }
            });

        // --- Step 8a: MPI halo exchange for the current time step (timeIdx = 1) ---
        gmpi.mpi_sendrecv_rt(rt.c, gdom, par, rt.n_mass, 1);

        // --- Step 8b: Also synchronize timeIdx = 0 to prevent uninitialized halo cells
        //              at the initial time step ---
        gmpi.mpi_sendrecv_rt(rt.c, gdom, par, rt.n_mass, 0);

        // --- Step 9: Apply boundary conditions (multi-species support) ---
        for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
        {
            for (int k = 0; k < rtgbc.rtgwbc.size(); k++)
            {
                rtgbc.rtgwbc[k].applyConcentrationBC(rt, gw, gdom, par, iSpec);
            }
        }

        // --- Step 10: Write initial output for reactive transport ---
        io.outputIniRT(rt, gdom, par, outFolder);

        flag = 1;
        return flag;
    }

    // =======================================================================
    // readRTFile()
    // =======================================================================
    /*!
     * @brief Read solute transport parameters from rttransportgw.input.
     *
     * Parses a colon-separated key-value file containing:
     *   - n_mass:              Number of dissolved species
     *   - rt_scheme:           Numerical scheme selector for transport
     *   - Up_Weighting_vplus:  Upwind weighting factor (positive velocity direction)
     *   - Up_Weighting_vminus: Upwind weighting factor (negative velocity direction)
     *   - ReactionModule:      Flag indicating whether reaction module is active (1) or not
     *   - Per-species blocks (indexed by "id") containing:
     *       - aq_mode / aq_val / aq_file:         Liquid-phase initial condition mode and value
     *       - solid_mode / solid_val / solid_file:  Solid-phase initial condition mode and value
     *       - diffusion_molecular:  Molecular diffusion coefficient [m^2/s]
     *       - alpha_L:              Longitudinal dispersivity [m]
     *       - alpha_T:              Transverse dispersivity [m]
     *
     * After reading, validates that all mandatory parameters and per-species
     * initial conditions have been set, and allocates Kokkos device arrays
     * for the transport parameters.
     *
     * @param[in]  fNameIn  Full path to the input file
     * @param[in]  rt       Reactive transport state (populated on output)
     * @param[in]  par      Parallel decomposition (used for masterproc output)
     * @return 1 on success, 0 on failure.  Calls exit(-1) on missing mandatory keys.
     */
    int readRTFile(std::string fNameIn, RTStateGW &rt, Parallel &par)
    {
        // Initialize all scalar parameters to sentinel value -999
        rt.ReactionModule = -999;
        rt.n_mass = -999;
        rt.rt_scheme = -999;
        rt.Up_Weighting_vplus = -999;
        rt.Up_Weighting_vminus = -999;

        // Initialize per-species, per-phase initial condition arrays
        // (these will be re-allocated once n_mass is known)
        rt.RT_Init_Aq_Mode.clear();
        rt.RT_Init_Aq_Const.clear();
        rt.RT_Init_Aq_File.clear();
        rt.RT_Init_Solid_Mode.clear();
        rt.RT_Init_Solid_Const.clear();
        rt.RT_Init_Solid_File.clear();

        // Track the species ID currently being parsed
        int current_species_id = -1;

        // Temporary host-side storage for per-species transport parameters
        std::vector<real> temp_diffusion_molecular;
        std::vector<real> temp_alpha_L;
        std::vector<real> temp_alpha_T;

        // Read in colon-separated key: value file line by line
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
                        pline.value >> rt.n_mass;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: RTM Subsurface n_mass:" << rt.n_mass << std::endl;

                        // Allocate per-species initial condition arrays now that n_mass is known
                        rt.RT_Init_Aq_Mode.resize(rt.n_mass, -999);
                        rt.RT_Init_Aq_Const.resize(rt.n_mass, 0.0);
                        rt.RT_Init_Aq_File.resize(rt.n_mass, "");
                        rt.RT_Init_Solid_Mode.resize(rt.n_mass, -999);
                        rt.RT_Init_Solid_Const.resize(rt.n_mass, 0.0);
                        rt.RT_Init_Solid_File.resize(rt.n_mass, "");

                        // Pre-allocate transport parameter arrays
                        temp_diffusion_molecular.resize(rt.n_mass, 0.0);
                        temp_alpha_L.resize(rt.n_mass, 0.0);
                        temp_alpha_T.resize(rt.n_mass, 0.0);
                    }
                    else if (!strcmp("id", pline.key.c_str()))
                    {
                        // Read species ID -- subsequent per-species keys apply to this ID
                        pline.value >> current_species_id;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: Reading species " << current_species_id << " parameters" << std::endl;
                    }
                    else if (!strcmp("aq_mode", pline.key.c_str()))
                    {
                        // Liquid-phase initial condition mode: 0 = constant value, 1 = read from file
                        if (current_species_id >= 0 && current_species_id < rt.n_mass)
                        {
                            pline.value >> rt.RT_Init_Aq_Mode[current_species_id];
                        }
                    }
                    else if (!strcmp("aq_val", pline.key.c_str()))
                    {
                        // Liquid-phase constant concentration value [mg/L or user-defined]
                        if (current_species_id >= 0 && current_species_id < rt.n_mass)
                        {
                            pline.value >> rt.RT_Init_Aq_Const[current_species_id];
                            rt.RT_Init_Aq_File[current_species_id] = "";
                        }
                    }
                    else if (!strcmp("aq_file", pline.key.c_str()))
                    {
                        // Liquid-phase spatial initial condition file name
                        if (current_species_id >= 0 && current_species_id < rt.n_mass)
                        {
                            pline.value >> rt.RT_Init_Aq_File[current_species_id];
                            rt.RT_Init_Aq_Const[current_species_id] = 0.0;
                        }
                    }
                    else if (!strcmp("solid_mode", pline.key.c_str()))
                    {
                        // Solid-phase initial condition mode: 0 = constant value, 1 = read from file
                        if (current_species_id >= 0 && current_species_id < rt.n_mass)
                        {
                            pline.value >> rt.RT_Init_Solid_Mode[current_species_id];
                            if (par.masterproc)
                                std::cerr << GOK "GWRTM:   solid_mode=" << rt.RT_Init_Solid_Mode[current_species_id] << " for species " << current_species_id << std::endl;
                        }
                    }
                    else if (!strcmp("solid_val", pline.key.c_str()))
                    {
                        // Solid-phase constant concentration value [mg/kg or user-defined]
                        if (current_species_id >= 0 && current_species_id < rt.n_mass)
                        {
                            pline.value >> rt.RT_Init_Solid_Const[current_species_id];
                            rt.RT_Init_Solid_File[current_species_id] = "";
                            if (par.masterproc)
                                std::cerr << GOK "GWRTM:   solid_val=" << rt.RT_Init_Solid_Const[current_species_id] << " for species " << current_species_id << std::endl;
                        }
                    }
                    else if (!strcmp("solid_file", pline.key.c_str()))
                    {
                        // Solid-phase spatial initial condition file name
                        if (current_species_id >= 0 && current_species_id < rt.n_mass)
                        {
                            pline.value >> rt.RT_Init_Solid_File[current_species_id];
                            rt.RT_Init_Solid_Const[current_species_id] = 0.0;
                            if (par.masterproc)
                                std::cerr << GOK "GWRTM:   solid_file=" << rt.RT_Init_Solid_File[current_species_id] << " for species " << current_species_id << std::endl;
                        }
                    }
                    else if (!strcmp("diffusion_molecular", pline.key.c_str()))
                    {
                        // Molecular diffusion coefficient for the current species [m^2/s]
                        if (current_species_id >= 0 && current_species_id < rt.n_mass)
                        {
                            pline.value >> temp_diffusion_molecular[current_species_id];
                        }
                    }
                    else if (!strcmp("alpha_L", pline.key.c_str()))
                    {
                        // Longitudinal dispersivity for the current species [m]
                        if (current_species_id >= 0 && current_species_id < rt.n_mass)
                        {
                            pline.value >> temp_alpha_L[current_species_id];
                        }
                    }
                    else if (!strcmp("alpha_T", pline.key.c_str()))
                    {
                        // Transverse dispersivity for the current species [m]
                        if (current_species_id >= 0 && current_species_id < rt.n_mass)
                        {
                            pline.value >> temp_alpha_T[current_species_id];
                        }
                    }
                    else if (!strcmp("rt_scheme", pline.key.c_str()))
                    {
                        pline.value >> rt.rt_scheme;
                    }
                    else if (!strcmp("Up_Weighting_vplus", pline.key.c_str()))
                    {
                        pline.value >> rt.Up_Weighting_vplus;
                    }
                    else if (!strcmp("Up_Weighting_vminus", pline.key.c_str()))
                    {
                        pline.value >> rt.Up_Weighting_vminus;
                    }
                    else if (!strcmp("ReactionModule", pline.key.c_str()))
                    {
                        pline.value >> rt.ReactionModule;
                    }
                }
            }
        }
        else
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "GWRTM: Unable to open " << fNameIn << "\n";
                return 0;
            }
        }

        // --- Validation: ensure all mandatory scalar parameters were set ---
        if (rt.n_mass == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "GWRTM: key n_mass not set." << std::endl;
            exit(-1);
        }

        if (rt.rt_scheme == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "GWRTM: key rt_scheme not set." << std::endl;
            exit(-1);
        }
        if (rt.Up_Weighting_vplus == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "GWRTM: key Up_Weighting_vplus not set." << std::endl;
            exit(-1);
        }
        if (rt.Up_Weighting_vminus == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "GWRTM: key Up_Weighting_vminus not set." << std::endl;
            exit(-1);
        }
        if (rt.ReactionModule == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "GWRTM: key ReactionModule not set." << std::endl;
            exit(-1);
        }

        // --- Validation: ensure liquid-phase initial condition mode is set for every species ---
        for (int i = 0; i < rt.n_mass; i++)
        {
            if (rt.RT_Init_Aq_Mode[i] == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: aq_mode for species " << i << " not set." << std::endl;
                exit(-1);
            }
        }

        // --- Validation: ensure solid-phase initial condition mode is set for every species
        //     when the reaction module is enabled ---
        if (rt.ReactionModule == 1)
        {
            for (int i = 0; i < rt.n_mass; i++)
            {
                if (rt.RT_Init_Solid_Mode[i] == -999)
                {
                    if (par.masterproc)
                        std::cerr << RERROR "GWRTM: solid_mode for species " << i << " not set." << std::endl;
                    exit(-1);
                }
            }
        }

        // --- Print summary of parsed initial condition settings ---
        if (par.masterproc)
        {
            std::cerr << GOK "GWRTM: Initial conditions:" << std::endl;
            for (int i = 0; i < rt.n_mass; i++)
            {
                std::cerr << "  Species " << i << ": ";
                // Liquid-phase info
                if (rt.RT_Init_Aq_Mode[i] == 0)
                    std::cerr << "aq_mode=0(const), aq_val=" << rt.RT_Init_Aq_Const[i];
                else if (rt.RT_Init_Aq_Mode[i] == 1)
                    std::cerr << "aq_mode=1(file), aq_file=" << rt.RT_Init_Aq_File[i];
                // Solid-phase info (if set)
                if (rt.RT_Init_Solid_Mode[i] != -999)
                {
                    std::cerr << ", ";
                    if (rt.RT_Init_Solid_Mode[i] == 0)
                        std::cerr << "solid_mode=0(const), solid_val=" << rt.RT_Init_Solid_Const[i];
                    else if (rt.RT_Init_Solid_Mode[i] == 1)
                        std::cerr << "solid_mode=1(file), solid_file=" << rt.RT_Init_Solid_File[i];
                }
                std::cerr << std::endl;
            }
        }

        // --- Allocate Kokkos device arrays for transport parameters and copy values ---
        rt.diffusion_molecular = realArr("diffusion_molecular", rt.n_mass);
        rt.alpha_L = realArr("alpha_L", rt.n_mass);
        rt.alpha_T = realArr("alpha_T", rt.n_mass);

        for (int i = 0; i < rt.n_mass; i++)
        {
            rt.diffusion_molecular(i) = temp_diffusion_molecular[i];
            rt.alpha_L(i) = temp_alpha_L[i];
            rt.alpha_T(i) = temp_alpha_T[i];
        }

        // --- Print transport parameters summary ---
        if (par.masterproc)
        {
            std::cerr << GOK "GWRTM: Transport properties:" << std::endl;
            for (int i = 0; i < rt.n_mass; i++)
            {
                std::cerr << "  Species " << i << ": diffusion_molecular=" << rt.diffusion_molecular(i)
                          << ", alpha_L=" << rt.alpha_L(i)
                          << ", alpha_T=" << rt.alpha_T(i) << std::endl;
            }
            std::cerr << GOK "GWRTM: Transport parameters read\n";
        }
        return 1;
    }

    // =======================================================================
    // readRTSSFile()
    // =======================================================================
    /*!
     * @brief Read groundwater reactive transport source/sink terms from rtgwss.input.
     *
     * Each source/sink entry is defined by an ID and a polygon file that
     * determines the spatial extent.  Per-species source/sink types and
     * values (or time-series files) can be specified using the syntax:
     *   spec_X_sstype : <type>
     *   spec_X_ssval  : <constant_value>
     *   spec_X_ssfile : <time_series_filename>
     * where X is the zero-based species index.
     *
     * After parsing the input file, this function:
     *   1. Reads polygon geometry and maps internal grid cells via find_icells()
     *   2. Reads time-series files for species with time-varying source/sink
     *
     * @param[in]     fNameIn  Full path to the source/sink input file
     * @param[in]     gdom     Groundwater domain
     * @param[in]     par      Parallel decomposition
     * @param[in,out] rt       Reactive transport state (used for n_mass)
     * @param[in]     gw       Groundwater state
     * @param[in,out] ss       Source/sink container (populated on output)
     * @return 1 on success, 0 on failure
     */
    int readRTSSFile(std::string fNameIn, GwDomain &gdom, Parallel &par, RTStateGW &rt, GwState &gw, SourceSink &ss)
    {
        std::ifstream fInStream(fNameIn);
        std::string dir = fNameIn.substr(0, fNameIn.length() - 12); // "rtgwss.input" length = 12
        std::vector<std::string> polygonFile;
        std::vector<std::string> fullPathPoly;
        std::string line;
        PsLn pline;
        int nPoly;
        int sscount = 0, sscountFound = 0, iss = -2;

        if (fInStream.is_open())
        {
            while (std::getline(fInStream, line))
            {
                pline.line = line;
                pline.parse();
                if (!pline.key.empty())
                {
                    if (!strcmp("sscount", pline.key.c_str()))
                    {
                        pline.value >> sscount;
                        sscountFound = 1;
                        if (sscount < 1)
                            return 1;
                        ss.rtgwss.resize(sscount);
                        polygonFile.resize(sscount);
                        fullPathPoly.resize(sscount);
                        iss++;
                    }
                    else if (!strcmp("id", pline.key.c_str()))
                    {
                        iss++;
                        if (sscount > 0 && iss >= 0)
                        {
                            pline.value >> ss.rtgwss[iss].id;
                            // Pre-allocate per-species storage for this source/sink entry
                            ss.rtgwss[iss].spec_sstype.resize(rt.n_mass, -999);
                            ss.rtgwss[iss].spec_ssval_const.resize(rt.n_mass, 0.0);
                            ss.rtgwss[iss].spec_ssfile.resize(rt.n_mass, "");
                            ss.rtgwss[iss].has_ssfile.resize(rt.n_mass, 0);
                            ss.rtgwss[iss].spec_ts.resize(rt.n_mass);
                        }
                    }
                    else if (!strcmp("rtsstype", pline.key.c_str()) && iss >= 0)
                    {
                        if (sscount > 0)
                            pline.value >> ss.rtgwss[iss].rtsstype;
                    }
                    else if (!strcmp("polygon", pline.key.c_str()) && iss >= 0)
                    {
                        if (sscount > 0)
                            pline.value >> polygonFile[iss];
                    }
                    // --- Parse per-species source/sink parameters (format: spec_X_<param>) ---
                    else if (pline.key.length() > 5 && pline.key.substr(0, 5) == "spec_" && iss >= 0)
                    {
                        std::string spec_str = pline.key.substr(5);
                        size_t underscore_pos = spec_str.find('_');
                        if (underscore_pos != std::string::npos)
                        {
                            int spec_id = std::stoi(spec_str.substr(0, underscore_pos));
                            std::string param_type = spec_str.substr(underscore_pos + 1);

                            if (spec_id >= 0 && spec_id < rt.n_mass)
                            {
                                if (param_type == "sstype")
                                {
                                    pline.value >> ss.rtgwss[iss].spec_sstype[spec_id];
                                }
                                else if (param_type == "ssval")
                                {
                                    pline.value >> ss.rtgwss[iss].spec_ssval_const[spec_id];
                                }
                                else if (param_type == "ssfile" || param_type == "file")
                                {
                                    pline.value >> ss.rtgwss[iss].spec_ssfile[spec_id];
                                    ss.rtgwss[iss].has_ssfile[spec_id] = 1;
                                }
                            }
                        }
                    }
                }
            }
            fInStream.close();

            if (!sscountFound)
                return 0;

            // --- Parse polygon files to identify internal grid cells for each source/sink ---
            for (int k = 0; k < sscount; k++)
            {
                fullPathPoly[k] = dir + polygonFile[k];
                std::ifstream fPoly(fullPathPoly[k]);
                if (fPoly.is_open())
                {
                    fPoly.ignore(256, ' ');
                    fPoly >> nPoly;
                    realArr xPoly = realArr("xPoly", nPoly);
                    realArr yPoly = realArr("yPoly", nPoly);
                    realArr zPoly = realArr("zPoly", nPoly);
                    for (int i = 0; i < nPoly; i++)
                    {
                        fPoly >> xPoly(i) >> yPoly(i) >> zPoly(i);
                    }
                    int find_result = ss.rtgwss[k].find_icells(gw, ss.rtgwss[k].id, gdom, par, nPoly, xPoly, yPoly, zPoly);
                    if (find_result != 1)
                    {
                        if (par.masterproc)
                        {
                            std::cerr << RERROR << "GWRTM: Failed to find internal cells for RT source/sink " << ss.rtgwss[k].id << std::endl;
                        }
                        return 0;
                    }
                    fPoly.close();
                }
                else
                {
                    if (par.masterproc)
                    {
                        std::cerr << YEXC << "GWRTM: Warning - could not open polygon file " << fullPathPoly[k] << " for RT source/sink " << ss.rtgwss[k].id << std::endl;
                    }
                }
            }

            // --- Read time-series files for species with time-varying injection ---
            for (int k = 0; k < sscount; k++)
            {
                for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
                {
                    if (ss.rtgwss[k].spec_sstype[iSpec] == 1 && ss.rtgwss[k].has_ssfile[iSpec])
                    {
                        std::string fname = dir + ss.rtgwss[k].spec_ssfile[iSpec];
                        std::ifstream fts(fname);
                        if (fts.is_open())
                        {
                            int ndatat = 0;
                            fts.ignore(256, ' ');
                            fts >> ndatat;
                            if (ndatat > 0)
                                ss.rtgwss[k].spec_ts[iSpec].initialise(ndatat);
                            for (int i = 0; i < ndatat; i++)
                            {

                                fts >> ss.rtgwss[k].spec_ts[iSpec].time(i) >> ss.rtgwss[k].spec_ts[iSpec].value(i);
                            }
                            fts.close();
                        }
                    }
                }
            }
            if (par.masterproc)
                std::cout << GOK << "GWRTM: RT Source/Sink file parsed successfully." << std::endl;
        }
        return 1;
    }

    // =======================================================================
    // readRTReactionFile()
    // =======================================================================
    /*!
     * @brief Read reaction parameters from reactiongw.input.
     *
     * Parses the reaction parameter file which contains:
     *
     * Global parameters:
     *   - ReactModel_Type:          Reaction model selector (0 = sorption only, 1 = full nitrogen cycle)
     *   - Bulk_Density (rho_b):     Bulk density of the porous medium [kg/L or g/cm^3]
     *   - Nitrogen_Cycle_Simulation: Flag to enable nitrogen cycle reactions
     *
     * Per-species sorption parameters (one value per species):
     *   - Sorption_Type:  Sorption model (0 = none, 1 = linear/Kd, 2 = Freundlich, 3 = Langmuir)
     *   - Kd:             Linear distribution coefficient [L/kg]
     *   - Kf:             Freundlich affinity constant
     *   - Nf:             Freundlich exponent
     *   - Alpha_D:        Dual-domain mass transfer coefficient [1/s]
     *   - Beta_D:         Dual-domain immobile fraction [-]
     *   - Beta:           Langmuir sorption capacity parameter
     *   - Lambda_1:       First-order decay rate in liquid phase [1/s]
     *   - Lambda_2:       First-order decay rate in solid phase [1/s]
     *
     * Nitrification parameters (single value):
     *   - Rate_Max_Nit:   Maximum nitrification rate [mg/L/s]
     *   - K_Monod_NH4:    Half-saturation constant for NH4+ [mg/L]
     *   - K_Monod_DO:     Half-saturation constant for dissolved oxygen [mg/L]
     *
     * Denitrification parameters (single value):
     *   - Rate_Max_Denit: Maximum denitrification rate [mg/L/s]
     *   - K_Monod_NO3:    Half-saturation constant for NO3- [mg/L]
     *   - K_Monod_DOC:    Half-saturation constant for DOC [mg/L]
     *   - Ki_Inhib_DO:    DO inhibition constant [mg/L]
     *
     * Aerobic respiration parameters (single value):
     *   - Rate_Max_Hetero:       Maximum heterotrophic respiration rate [mg/L/s]
     *   - K_Monod_DOC_Aerobic:   Half-saturation constant for DOC [mg/L]
     *   - K_Monod_DO_Aerobic:    Half-saturation constant for DO [mg/L]
     *
     * Mineralization parameters (single value):
     *   - Rate_Max_Min:    Maximum mineralization rate [mg/L/s]
     *   - K_Monod_DON:     Half-saturation constant for DON [mg/L]
     *   - K_Monod_DO_Min:  Half-saturation constant for DO [mg/L]
     *
     * Reaeration and release parameters:
     *   - Ka_DO:     Reaeration rate coefficient [1/s]
     *   - DO_sat:    Saturated dissolved oxygen concentration [mg/L]
     *   - DOC_eq:    Equilibrium DOC concentration [mg/L]
     *   - DON_eq:    Equilibrium DON concentration [mg/L]
     *   - K_rel_DOC: DOC release rate coefficient [1/s]
     *   - K_rel_DON: DON release rate coefficient [1/s]
     *
     * Environmental correction factors:
     *   - Temp_Coeff_Theta: Arrhenius temperature correction coefficient [-]
     *   - Opt_Temp:         Optimal temperature [deg C]
     *   - Opt_pH:           Optimal pH [-]
     *
     * @param[in]     fNameIn  Full path to the reaction input file
     * @param[in,out] rt       Reactive transport state (populated on output)
     * @param[in]     par      Parallel decomposition (used for masterproc output)
     * @return 1 on success, 0 on failure
     */
    int readRTReactionFile(std::string fNameIn, RTStateGW &rt, Parallel &par)
    {
        // Initialize all scalar reaction parameters to sentinel value -999
        rt.ReactModel_Type = -999;
        rt.rho_b = -999;
        rt.Nitrogen_Cycle_Simulation = -999;

        // Nitrification reaction parameters
        rt.Rate_Max_Nit = -999;
        rt.K_Monod_NH4 = -999;
        rt.K_Monod_DO = -999;

        // Denitrification reaction parameters
        rt.Rate_Max_Denit = -999;
        rt.K_Monod_NO3 = -999;
        rt.K_Monod_DOC = -999;
        rt.Ki_Inhib_DO = -999;

        // Aerobic respiration reaction parameters
        rt.Rate_Max_Hetero = -999;
        rt.K_Monod_DOC_Aerobic = -999;
        rt.K_Monod_DO_Aerobic = -999;

        // Mineralization reaction parameters
        rt.Rate_Max_Min = -999;
        rt.K_Monod_DON = -999;
        rt.K_Monod_DO_Min = -999;

        // Reaeration, DOC release, and DON release parameters
        rt.Ka_DO = -999;
        rt.DO_sat = -999;
        rt.DOC_eq = -999;
        rt.DON_eq = -999;
        rt.K_rel_DOC = -999;
        rt.K_rel_DON = -999;

        // Environmental correction factors
        rt.Temp_Coeff_Theta = -999;
        rt.Opt_Temp = -999;
        rt.Opt_pH = -999;

        // Initialize all per-species array values to sentinel value -999
        for (int i = 0; i < rt.n_mass; i++)
        {
            rt.Sorption_Type(i) = -999;
            rt.Kd(i) = -999;
            rt.Kf(i) = -999;
            rt.Nf(i) = -999;
            rt.Alpha_D(i) = -999;
            rt.Beta_D(i) = -999;
            rt.Beta(i) = -999;
            rt.Lambda_1(i) = -999;
            rt.Lambda_2(i) = -999;
        }

        std::string strAux;
        // Read in colon-separated key: value file line by line
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
                    // ---- Global physical / environmental parameters ----
                    if (!strcmp("ReactModel_Type", pline.key.c_str()))
                    {
                        pline.value >> rt.ReactModel_Type;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: ReactModel_Type: " << rt.ReactModel_Type << std::endl;
                    }
                    else if (!strcmp("Bulk_Density", pline.key.c_str()))
                    {
                        pline.value >> rt.rho_b;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: Bulk_Density: " << rt.rho_b << std::endl;
                    }
                    else if (!strcmp("Nitrogen_Cycle_Simulation", pline.key.c_str()))
                    {
                        pline.value >> rt.Nitrogen_Cycle_Simulation;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: Nitrogen_Cycle_Simulation: " << rt.Nitrogen_Cycle_Simulation << std::endl;
                    }
                    else if (!strcmp("SPECIES_ID", pline.key.c_str()))
                    {
                        // Read species IDs (informational, used for debug output)
                        std::vector<int> species_ids;
                        int id;
                        while (pline.value >> id)
                        {
                            species_ids.push_back(id);
                        }
                        if (par.masterproc)
                        {
                            std::cerr << GOK "GWRTM: SPECIES_ID: ";
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
                        // Read species names (informational, used for debug output)
                        std::vector<std::string> species_names;
                        std::string name;
                        while (pline.value >> name)
                        {
                            // Remove semicolons from the parsed name
                            name.erase(std::remove(name.begin(), name.end(), ';'), name.end());
                            species_names.push_back(name);
                        }
                        if (par.masterproc)
                        {
                            std::cerr << GOK "GWRTM: SPECIES_Name: ";
                            for (size_t i = 0; i < species_names.size(); i++)
                            {
                                std::cerr << species_names[i];
                                if (i < species_names.size() - 1)
                                    std::cerr << "; ";
                            }
                            std::cerr << std::endl;
                        }
                    }
                    // ---- Sorption parameters (per-species, one value per species) ----
                    else if (!strcmp("Sorption_Type", pline.key.c_str()))
                    {
                        for (int i = 0; i < rt.n_mass; i++)
                        {
                            pline.value >> rt.Sorption_Type(i);
                        }
                        if (par.masterproc)
                        {
                            std::cerr << GOK "GWRTM: Sorption_Type: ";
                            for (int i = 0; i < rt.n_mass; i++)
                            {
                                std::cerr << rt.Sorption_Type(i);
                                if (i < rt.n_mass - 1)
                                    std::cerr << "; ";
                            }
                            std::cerr << std::endl;
                        }
                    }
                    else if (!strcmp("Kd", pline.key.c_str()))
                    {
                        for (int i = 0; i < rt.n_mass; i++)
                        {
                            pline.value >> rt.Kd(i);
                        }
                        if (par.masterproc)
                        {
                            std::cerr << GOK "GWRTM: Kd: ";
                            for (int i = 0; i < rt.n_mass; i++)
                            {
                                std::cerr << rt.Kd(i);
                                if (i < rt.n_mass - 1)
                                    std::cerr << "; ";
                            }
                            std::cerr << std::endl;
                        }
                    }
                    else if (!strcmp("Kf", pline.key.c_str()))
                    {
                        for (int i = 0; i < rt.n_mass; i++)
                        {
                            pline.value >> rt.Kf(i);
                        }
                        if (par.masterproc)
                        {
                            std::cerr << GOK "GWRTM: Kf: ";
                            for (int i = 0; i < rt.n_mass; i++)
                            {
                                std::cerr << rt.Kf(i);
                                if (i < rt.n_mass - 1)
                                    std::cerr << "; ";
                            }
                            std::cerr << std::endl;
                        }
                    }
                    else if (!strcmp("Nf", pline.key.c_str()))
                    {
                        for (int i = 0; i < rt.n_mass; i++)
                        {
                            pline.value >> rt.Nf(i);
                        }
                        if (par.masterproc)
                        {
                            std::cerr << GOK "GWRTM: Nf: ";
                            for (int i = 0; i < rt.n_mass; i++)
                            {
                                std::cerr << rt.Nf(i);
                                if (i < rt.n_mass - 1)
                                    std::cerr << "; ";
                            }
                            std::cerr << std::endl;
                        }
                    }
                    else if (!strcmp("Alpha_D", pline.key.c_str()))
                    {
                        for (int i = 0; i < rt.n_mass; i++)
                        {
                            pline.value >> rt.Alpha_D(i);
                        }
                        if (par.masterproc)
                        {
                            std::cerr << GOK "GWRTM: Alpha_D: ";
                            for (int i = 0; i < rt.n_mass; i++)
                            {
                                std::cerr << rt.Alpha_D(i);
                                if (i < rt.n_mass - 1)
                                    std::cerr << "; ";
                            }
                            std::cerr << std::endl;
                        }
                    }
                    else if (!strcmp("Beta_D", pline.key.c_str()))
                    {
                        for (int i = 0; i < rt.n_mass; i++)
                        {
                            pline.value >> rt.Beta_D(i);
                        }
                        if (par.masterproc)
                        {
                            std::cerr << GOK "GWRTM: Beta_D: ";
                            for (int i = 0; i < rt.n_mass; i++)
                            {
                                std::cerr << rt.Beta_D(i);
                                if (i < rt.n_mass - 1)
                                    std::cerr << "; ";
                            }
                            std::cerr << std::endl;
                        }
                    }
                    else if (!strcmp("Beta", pline.key.c_str()))
                    {
                        for (int i = 0; i < rt.n_mass; i++)
                        {
                            pline.value >> rt.Beta(i);
                        }
                        if (par.masterproc)
                        {
                            std::cerr << GOK "GWRTM: Beta: ";
                            for (int i = 0; i < rt.n_mass; i++)
                            {
                                std::cerr << rt.Beta(i);
                                if (i < rt.n_mass - 1)
                                    std::cerr << "; ";
                            }
                            std::cerr << std::endl;
                        }
                    }
                    else if (!strcmp("Lambda_1", pline.key.c_str()))
                    {
                        for (int i = 0; i < rt.n_mass; i++)
                        {
                            pline.value >> rt.Lambda_1(i);
                        }
                        if (par.masterproc)
                        {
                            std::cerr << GOK "GWRTM: Lambda_1: ";
                            for (int i = 0; i < rt.n_mass; i++)
                            {
                                std::cerr << rt.Lambda_1(i);
                                if (i < rt.n_mass - 1)
                                    std::cerr << "; ";
                            }
                            std::cerr << std::endl;
                        }
                    }
                    else if (!strcmp("Lambda_2", pline.key.c_str()))
                    {
                        for (int i = 0; i < rt.n_mass; i++)
                        {
                            pline.value >> rt.Lambda_2(i);
                        }
                        if (par.masterproc)
                        {
                            std::cerr << GOK "GWRTM: Lambda_2: ";
                            for (int i = 0; i < rt.n_mass; i++)
                            {
                                std::cerr << rt.Lambda_2(i);
                                if (i < rt.n_mass - 1)
                                    std::cerr << "; ";
                            }
                            std::cerr << std::endl;
                        }
                    }
                    // ---- Nitrification parameters (single values) ----
                    else if (!strcmp("Rate_Max_Nit", pline.key.c_str()))
                    {
                        pline.value >> rt.Rate_Max_Nit;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: Rate_Max_Nit: " << rt.Rate_Max_Nit << std::endl;
                    }
                    else if (!strcmp("K_Monod_NH4", pline.key.c_str()))
                    {
                        pline.value >> rt.K_Monod_NH4;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: K_Monod_NH4: " << rt.K_Monod_NH4 << std::endl;
                    }
                    else if (!strcmp("K_Monod_DO", pline.key.c_str()))
                    {
                        pline.value >> rt.K_Monod_DO;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: K_Monod_DO: " << rt.K_Monod_DO << std::endl;
                    }

                    // ---- Denitrification parameters (single values) ----
                    else if (!strcmp("Rate_Max_Denit", pline.key.c_str()))
                    {
                        pline.value >> rt.Rate_Max_Denit;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: Rate_Max_Denit: " << rt.Rate_Max_Denit << std::endl;
                    }
                    else if (!strcmp("K_Monod_NO3", pline.key.c_str()))
                    {
                        pline.value >> rt.K_Monod_NO3;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: K_Monod_NO3: " << rt.K_Monod_NO3 << std::endl;
                    }
                    else if (!strcmp("K_Monod_DOC", pline.key.c_str()))
                    {
                        pline.value >> rt.K_Monod_DOC;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: K_Monod_DOC: " << rt.K_Monod_DOC << std::endl;
                    }
                    else if (!strcmp("Ki_Inhib_DO", pline.key.c_str()))
                    {
                        pline.value >> rt.Ki_Inhib_DO;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: Ki_Inhib_DO: " << rt.Ki_Inhib_DO << std::endl;
                    }

                    // ---- Aerobic respiration parameters ----
                    else if (!strcmp("Rate_Max_Hetero", pline.key.c_str()))
                    {
                        pline.value >> rt.Rate_Max_Hetero;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: Rate_Max_Hetero: " << rt.Rate_Max_Hetero << std::endl;
                    }
                    else if (!strcmp("K_Monod_DOC_Aerobic", pline.key.c_str()))
                    {
                        pline.value >> rt.K_Monod_DOC_Aerobic;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: K_Monod_DOC_Aerobic: " << rt.K_Monod_DOC_Aerobic << std::endl;
                    }
                    else if (!strcmp("K_Monod_DO_Aerobic", pline.key.c_str()))
                    {
                        pline.value >> rt.K_Monod_DO_Aerobic;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: K_Monod_DO_Aerobic: " << rt.K_Monod_DO_Aerobic << std::endl;
                    }

                    // ---- Mineralization parameters ----
                    else if (!strcmp("Rate_Max_Min", pline.key.c_str()))
                    {
                        pline.value >> rt.Rate_Max_Min;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: Rate_Max_Min: " << rt.Rate_Max_Min << std::endl;
                    }
                    else if (!strcmp("K_Monod_DON", pline.key.c_str()))
                    {
                        pline.value >> rt.K_Monod_DON;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: K_Monod_DON: " << rt.K_Monod_DON << std::endl;
                    }
                    else if (!strcmp("K_Monod_DO_Min", pline.key.c_str()))
                    {
                        pline.value >> rt.K_Monod_DO_Min;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: K_Monod_DO_Min: " << rt.K_Monod_DO_Min << std::endl;
                    }

                    // ---- Reaeration, DOC release, and DON release parameters ----
                    else if (!strcmp("Ka_DO", pline.key.c_str()))
                    {
                        pline.value >> rt.Ka_DO;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: Ka_DO: " << rt.Ka_DO << std::endl;
                    }
                    else if (!strcmp("DO_sat", pline.key.c_str()))
                    {
                        pline.value >> rt.DO_sat;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: DO_sat: " << rt.DO_sat << std::endl;
                    }
                    else if (!strcmp("DOC_eq", pline.key.c_str()))
                    {
                        pline.value >> rt.DOC_eq;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: DOC_eq: " << rt.DOC_eq << std::endl;
                    }
                    else if (!strcmp("DON_eq", pline.key.c_str()))
                    {
                        pline.value >> rt.DON_eq;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: DON_eq: " << rt.DON_eq << std::endl;
                    }
                    else if (!strcmp("K_rel_DOC", pline.key.c_str()))
                    {
                        pline.value >> rt.K_rel_DOC;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: K_rel_DOC: " << rt.K_rel_DOC << std::endl;
                    }
                    else if (!strcmp("K_rel_DON", pline.key.c_str()))
                    {
                        pline.value >> rt.K_rel_DON;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: K_rel_DON: " << rt.K_rel_DON << std::endl;
                    }

                    // ---- Environmental correction factors (single values) ----
                    else if (!strcmp("Temp_Coeff_Theta", pline.key.c_str()))
                    {
                        pline.value >> rt.Temp_Coeff_Theta;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: Temp_Coeff_Theta: " << rt.Temp_Coeff_Theta << std::endl;
                    }
                    else if (!strcmp("Opt_Temp", pline.key.c_str()))
                    {
                        pline.value >> rt.Opt_Temp;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: Opt_Temp: " << rt.Opt_Temp << std::endl;
                    }
                    else if (!strcmp("Opt_pH", pline.key.c_str()))
                    {
                        pline.value >> rt.Opt_pH;
                        if (par.masterproc)
                            std::cerr << GOK "GWRTM: Opt_pH: " << rt.Opt_pH << std::endl;
                    }
                }
            }
        }
        else
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "GWRTM: Unable to open " << fNameIn << "\n";
                return 0;
            }
        }

        // ====== Mandatory parameter checks ======

        // Check essential global parameters
        if (rt.ReactModel_Type == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "GWRTM: key ReactModel_Type not set." << std::endl;
            return 0;
        }
        if (rt.rho_b == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "GWRTM: key Bulk_Density not set." << std::endl;
            return 0;
        }
        if (rt.Nitrogen_Cycle_Simulation == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "GWRTM: key Nitrogen_Cycle_Simulation not set." << std::endl;
            return 0;
        }

        // Check per-species Sorption_Type
        for (int i = 0; i < rt.n_mass; i++)
        {
            if (rt.Sorption_Type(i) == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: Sorption_Type[" << i << "] not set." << std::endl;
                return 0;
            }
        }

        // ====== Reaction parameter checks (when ReactModel_Type == 1, i.e., full nitrogen cycle) ======
        if (rt.ReactModel_Type == 1)
        {
            // Nitrification parameter validation
            if (rt.Rate_Max_Nit == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key Rate_Max_Nit not set." << std::endl;
                return 0;
            }
            if (rt.K_Monod_NH4 == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key K_Monod_NH4 not set." << std::endl;
                return 0;
            }
            if (rt.K_Monod_DO == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key K_Monod_DO not set." << std::endl;
                return 0;
            }

            // Denitrification parameter validation
            if (rt.Rate_Max_Denit == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key Rate_Max_Denit not set." << std::endl;
                return 0;
            }
            if (rt.K_Monod_NO3 == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key K_Monod_NO3 not set." << std::endl;
                return 0;
            }
            if (rt.K_Monod_DOC == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key K_Monod_DOC not set." << std::endl;
                return 0;
            }
            if (rt.Ki_Inhib_DO == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key Ki_Inhib_DO not set." << std::endl;
                return 0;
            }

            // Aerobic respiration parameter validation
            if (rt.Rate_Max_Hetero == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key Rate_Max_Hetero not set." << std::endl;
                return 0;
            }
            if (rt.K_Monod_DOC_Aerobic == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key K_Monod_DOC_Aerobic not set." << std::endl;
                return 0;
            }
            if (rt.K_Monod_DO_Aerobic == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key K_Monod_DO_Aerobic not set." << std::endl;
                return 0;
            }

            // Mineralization parameter validation
            if (rt.Rate_Max_Min == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key Rate_Max_Min not set." << std::endl;
                return 0;
            }
            if (rt.K_Monod_DON == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key K_Monod_DON not set." << std::endl;
                return 0;
            }
            if (rt.K_Monod_DO_Min == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key K_Monod_DO_Min not set." << std::endl;
                return 0;
            }

            // Reaeration, DOC/DON release parameter validation
            if (rt.Ka_DO == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key Ka_DO not set." << std::endl;
                return 0;
            }
            if (rt.DO_sat == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key DO_sat not set." << std::endl;
                return 0;
            }
            if (rt.DOC_eq == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key DOC_eq not set." << std::endl;
                return 0;
            }
            if (rt.DON_eq == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key DON_eq not set." << std::endl;
                return 0;
            }
            if (rt.K_rel_DOC == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key K_rel_DOC not set." << std::endl;
                return 0;
            }
            if (rt.K_rel_DON == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key K_rel_DON not set." << std::endl;
                return 0;
            }

            // Environmental correction factor validation
            if (rt.Temp_Coeff_Theta == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key Temp_Coeff_Theta not set." << std::endl;
                return 0;
            }
            if (rt.Opt_Temp == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key Opt_Temp not set." << std::endl;
                return 0;
            }
            if (rt.Opt_pH == -999)
            {
                if (par.masterproc)
                    std::cerr << RERROR "GWRTM: key Opt_pH not set." << std::endl;
                return 0;
            }
        }

        // --- Print debug summary of reaction parameters ---
        if (par.masterproc)
        {
            std::cerr << GOK "GWRTM: RTM Reaction parameters read\n";
            std::cout << "GWRTM: DEBUG: === Reaction Parameters Check ===" << std::endl;
            std::cout << "GWRTM: DEBUG: ReactModel_Type = " << rt.ReactModel_Type << std::endl;
            std::cout << "GWRTM: DEBUG: Nitrogen_Cycle_Simulation = " << rt.Nitrogen_Cycle_Simulation << std::endl;
            std::cout << "GWRTM: DEBUG: n_mass = " << rt.n_mass << std::endl;
            std::cout << "GWRTM: DEBUG: Rate_Max_Nit(0) = " << rt.Rate_Max_Nit << std::endl;
            std::cout << "GWRTM: DEBUG: K_Monod_NH4(0) = " << rt.K_Monod_NH4 << std::endl;
            std::cout << "GWRTM: DEBUG: ===============================" << std::endl;
        }

        return 1;
    }

    // =======================================================================
    // readRTBCFile()
    // =======================================================================
    /*!
     * @brief Read reactive transport boundary conditions from rtgwbc.input.
     *
     * Parses a colon-separated key-value file that defines one or more
     * boundary conditions.  Each boundary is specified by an ID block
     * containing:
     *   - rtbctype:  Global boundary type (e.g., Dirichlet, Neumann, Cauchy)
     *   - polygon:   Polygon file defining the spatial extent of the boundary
     *   - direction: Face direction (1=+x, 2=-x, 3=+y, 4=-y, 5=+z, 6=-z)
     *   - bcvals:    Constant boundary value (applied uniformly)
     *
     * Per-species boundary conditions use the syntax:
     *   spec_X_bctype : <boundary_type>
     *   spec_X_bcval  : <constant_value>     (for CONST-type boundaries)
     *   spec_X_bcfile : <filename>            (for time-varying or spatial BCs)
     * where X is the zero-based species index.
     *
     * After parsing, this function:
     *   1. Reads polygon geometry and maps boundary cells via find_bcells()
     *   2. Reads global time-series files for time-varying boundaries
     *   3. Reads spatially distributed boundary data files
     *   4. Reads per-species spatial boundary files
     *   5. Reads per-species time-series boundary files
     *
     * @param[in]     fNameIn  Full path to the boundary condition input file
     * @param[in]     gdom     Groundwater domain (grid, geometry)
     * @param[in,out] rtgbc    Reactive transport boundary condition container
     * @param[in]     par      Parallel decomposition
     * @param[in,out] rt       Reactive transport state (used for n_mass)
     * @param[in]     gw       Groundwater state (used for cell mapping)
     * @return 1 on success, 0 on failure
     */
    int readRTBCFile(std::string fNameIn, GwDomain &gdom, RTSubsurfaceBoundaries &rtgbc, Parallel &par, RTStateGW &rt, GwState &gw)
    {
        std::ifstream fInStream(fNameIn);
        std::string dir;
        std::vector<std::string> polygonFile;
        std::vector<std::string> fullPathPoly;
        std::vector<std::string> tsFile, bcFile;
        std::vector<int> hasbcfile;
        std::string line;
        PsLn pline;
        int nPoly, nPoly3D;
        dir = fNameIn.substr(0, fNameIn.length() - 12);         // Extract directory: "rtgwbc.input" is 12 chars
        int bccount = 0, bccountFound = 0, ibc = -2, ndata = 0; // ibc starts at -2; incremented to -1 on bccount, then 0,1,... on each id
        real val;

        // --- Parse the rtgwbc.input file ---
        if (fInStream.is_open())
        {
            while (std::getline(fInStream, line))
            {
                pline.line = line;
                pline.parse();
                if (!pline.key.empty())
                {
                    // Read the number of boundaries
                    if (!strcmp("bccount", pline.key.c_str()))
                    {
                        pline.value >> bccount;
                        bccountFound = 1;
                        if (bccount < 1)
                        {
                            std::cout << YEXC << "GWRTM: RTM rtgwbc.input indicates zero external boundaries." << std::endl;
                            return 1;
                        }
                        rtgbc.rtgwbc.resize(bccount);
                        rtgbc.id.resize(bccount);
                        polygonFile.resize(bccount);
                        fullPathPoly.resize(bccount);
                        tsFile.resize(bccount);
                        bcFile.resize(bccount);
                        hasbcfile.resize(bccount);
                        ibc++; // ibc should be set to -1
                    }
                    else if (!strcmp("id", pline.key.c_str()))
                    {
                        ibc++; // Advances to 0 for the first boundary, 1 for the second, etc.
                        hasbcfile[ibc] = 0;
                        if (bccount > 0 && ibc >= 0)
                        {
                            pline.value >> rtgbc.id[ibc];
                            // Initialize multi-species boundary condition data structures
                            rtgbc.rtgwbc[ibc].spec_bctype.resize(rt.n_mass, -999);
                            rtgbc.rtgwbc[ibc].spec_bcvals.resize(rt.n_mass);
                            rtgbc.rtgwbc[ibc].spec_bcfile.resize(rt.n_mass, "");
                            rtgbc.rtgwbc[ibc].spec_ts.resize(rt.n_mass);
                            rtgbc.rtgwbc[ibc].spec_bcval_const.resize(rt.n_mass, 0.0);
                            rtgbc.rtgwbc[ibc].has_bcfile.resize(rt.n_mass, 0);
                        }
                    }
                    else if (!strcmp("rtbctype", pline.key.c_str()) && ibc >= 0)
                    {
                        if (bccount > 0)
                            pline.value >> rtgbc.rtgwbc[ibc].rtbctype;
                    }
                    else if (!strcmp("polygon", pline.key.c_str()) && ibc >= 0)
                    {
                        if (bccount > 0)
                            pline.value >> polygonFile[ibc];
                    }
                    else if (!strcmp("direction", pline.key.c_str()) && ibc >= 0)
                    {
                        if (bccount > 0)
                        {
                            pline.value >> rtgbc.rtgwbc[ibc].direction;
                            // Determine the number of boundary data points based on face direction
                            // direction: 1=+x, 2=-x, 3=+y, 4=-y, 5=+z, 6=-z
                            if (rtgbc.rtgwbc[ibc].direction == 1 || rtgbc.rtgwbc[ibc].direction == 2)
                            {
                                ndata = gdom.ny_glob * gdom.nz;
                            }
                            else if (rtgbc.rtgwbc[ibc].direction == 3 || rtgbc.rtgwbc[ibc].direction == 4)
                            {
                                ndata = gdom.nx_glob * gdom.nz;
                            }
                            else if (rtgbc.rtgwbc[ibc].direction == 5 || rtgbc.rtgwbc[ibc].direction == 6)
                            {
                                ndata = gdom.nx_glob * gdom.ny_glob;
                            }
                            else
                            {
                                std::cerr << RERROR << "GWRTM: In rtgwbc.input: direction must be 1, 2, 3, 4, 5 or 6 " << std::endl;
                            }
                        }
                    }
                    else if (!strcmp("bcvals", pline.key.c_str()) && ibc >= 0)
                    {
                        if (bccount > 0)
                        {
                            if (ndata <= 0)
                            {
                                std::cerr << RERROR << "GWRTM: In rtgwbc.input: direction should on top of bcvals " << std::endl;
                            }
                            rtgbc.rtgwbc[ibc].bcvals = realArr("bcvals", ndata);
                            pline.value >> val;
                            for (int idx = 0; idx < ndata; idx++)
                            {
                                rtgbc.rtgwbc[ibc].bcvals(idx) = val;
                            }
                        }
                    }
                    // --- Parse per-species boundary condition parameters (format: spec_X_<param>) ---
                    else if (pline.key.length() > 5 && pline.key.substr(0, 5) == "spec_" && ibc >= 0)
                    {
                        // Parse species ID and parameter type from the key
                        // Format: spec_X_bctype, spec_X_bcval, spec_X_bcfile
                        std::string spec_str = pline.key.substr(5); // Strip "spec_" prefix
                        size_t underscore_pos = spec_str.find('_');
                        if (underscore_pos != std::string::npos)
                        {
                            int spec_id = std::stoi(spec_str.substr(0, underscore_pos));
                            std::string param_type = spec_str.substr(underscore_pos + 1);

                            if (spec_id >= 0 && spec_id < rt.n_mass)
                            {
                                if (param_type == "bctype")
                                {
                                    pline.value >> rtgbc.rtgwbc[ibc].spec_bctype[spec_id];
                                    if (par.masterproc)
                                        std::cerr << GOK "GWRTM:   Boundary id=" << rtgbc.id[ibc] << ", Species " << spec_id << ", bctype=" << rtgbc.rtgwbc[ibc].spec_bctype[spec_id] << std::endl;
                                }
                                else if (param_type == "bcval")
                                {
                                    // Validate that bcval contains a numeric value, not a filename
                                    std::string val_str = pline.value.str();
                                    if (val_str.find('.') != std::string::npos &&
                                        (val_str.find(".input") != std::string::npos

                                         ))
                                    {
                                        if (par.masterproc)
                                        {
                                            std::cerr << RERROR "GWRTM:   Boundary id=" << rtgbc.id[ibc] << ", Species " << spec_id
                                                      << ": bcval should be a numeric value, not a filename (" << val_str << ")." << std::endl;
                                        }
                                        return 0;
                                    }

                                    pline.value >> rtgbc.rtgwbc[ibc].spec_bcval_const[spec_id];
                                    if (par.masterproc)
                                        std::cerr << GOK "GWRTM:   Boundary id=" << rtgbc.id[ibc] << ", Species " << spec_id << ", bcval=" << rtgbc.rtgwbc[ibc].spec_bcval_const[spec_id] << std::endl;

                                    // Validate that time-varying boundary types do not use bcval (should use bcfile)
                                    if (rtgbc.rtgwbc[ibc].spec_bctype[spec_id] == SUB_RT_BC_Dirichlet_T ||
                                        rtgbc.rtgwbc[ibc].spec_bctype[spec_id] == SUB_RT_BC_Neumann_T ||
                                        rtgbc.rtgwbc[ibc].spec_bctype[spec_id] == SUB_RT_BC_Cauchy_T)
                                    {
                                        if (par.masterproc)
                                        {
                                            std::cerr << RERROR "GWRTM:   Boundary id=" << rtgbc.id[ibc] << ", Species " << spec_id << ": bcval specified for T boundary type (type="
                                                      << rtgbc.rtgwbc[ibc].spec_bctype[spec_id] << "). Use bcfile instead." << std::endl;
                                        }
                                        return 0;
                                    }
                                }
                                else if (param_type == "bcfile" || param_type == "file")
                                {
                                    pline.value >> rtgbc.rtgwbc[ibc].spec_bcfile[spec_id];
                                    rtgbc.rtgwbc[ibc].has_bcfile[spec_id] = 1;
                                    if (par.masterproc)
                                        std::cerr << GOK "GWRTM:   Boundary id=" << rtgbc.id[ibc] << ", Species " << spec_id << ", " << param_type << "=" << rtgbc.rtgwbc[ibc].spec_bcfile[spec_id] << std::endl;

                                    // Validate that constant boundary types do not use bcfile (should use bcval)
                                    if (rtgbc.rtgwbc[ibc].spec_bctype[spec_id] == SUB_RT_BC_Dirichlet_CONST ||
                                        rtgbc.rtgwbc[ibc].spec_bctype[spec_id] == SUB_RT_BC_Neumann_CONST ||
                                        rtgbc.rtgwbc[ibc].spec_bctype[spec_id] == SUB_RT_BC_Cauchy_CONST)
                                    {
                                        if (par.masterproc)
                                        {
                                            std::cerr << RERROR "GWRTM:   Boundary id=" << rtgbc.id[ibc] << ", Species " << spec_id << ": " << param_type << " specified for CONST boundary type (type="
                                                      << rtgbc.rtgwbc[ibc].spec_bctype[spec_id] << "). Use bcval instead." << std::endl;
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
                            std::cerr << RERROR << "GWRTM: In rtgwbc.input: Key " << pline.key << " not understood." << std::endl;
                            return 0;
                        }
                    }
                    if (ibc < -1)
                    {
                        if (par.masterproc)
                        {
                            std::cerr << RERROR << "GWRTM: No boundaries defined in rtgwbc.input, number of boundaries not defined, or 'id' key not found." << std::endl;
                            return 0;
                        }
                    }
                }
            } // end while
            fInStream.close();
            if (!bccountFound)
            {
                std::cerr << RERROR << "GWRTM: Number of boundaries not defined in rtgwbc.input. Please define 'bccount'" << std::endl;
                return 0;
            }
            else
            {
                if (par.masterproc)
                {
                    std::cerr << GOK "GWRTM: RTM Subsurface BC set\n";
                }
            }

            // --- Allocate and initialize per-species boundary value arrays ---
            for (int k = 0; k < bccount; k++)
            {
                for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
                {
                    // If a boundary file is specified, allocate the array (to be populated later)
                    if (rtgbc.rtgwbc[k].has_bcfile[iSpec] && ndata > 0)
                    {
                        rtgbc.rtgwbc[k].spec_bcvals[iSpec] = realArr("spec_bcvals", ndata);
                    }
                    // Otherwise, if a constant boundary type is set, fill the array with the constant value
                    else if (rtgbc.rtgwbc[k].spec_bctype[iSpec] != -999 && ndata > 0)
                    {
                        rtgbc.rtgwbc[k].spec_bcvals[iSpec] = realArr("spec_bcvals", ndata);
                        real const_val = rtgbc.rtgwbc[k].spec_bcval_const[iSpec];
                        for (int idx = 0; idx < ndata; idx++)
                        {
                            rtgbc.rtgwbc[k].spec_bcvals[iSpec](idx) = const_val;
                        }
                    }
                }
            }
        }
        else
        {
            if (par.masterproc)
            {
                std::cerr << YEXC << "GWRTM: rtgwbc.input not found. Default boundaries used." << std::endl;
            }
        }

        // --- Read polygon files to identify boundary cells ---
        for (int k = 0; k < polygonFile.size(); k++)
        {
            fullPathPoly[k] = dir + polygonFile[k];

            std::ifstream fPoly(fullPathPoly[k]);
            // Read in the k-th polygon
            if (fPoly.is_open())
            {
                fPoly.ignore(256, ' ');
                fPoly >> nPoly;
                realArr xPoly = realArr("xPoly", nPoly);
                realArr yPoly = realArr("yPoly", nPoly);
                for (int i = 0; i < nPoly; i++)
                {
                    if (!fPoly.fail() && !fPoly.eof())
                    {
                        fPoly >> xPoly(i) >> yPoly(i);
#if SERGHEI_DEBUG_BOUNDARY
                        std::cout << GGD << "GWRTM: RTM subbc polygon " << k << ". Point " << i << "/" << nPoly << "\t" << xPoly(i) << "\t" << yPoly(i) << std::endl;
#endif
                    }
                    else
                    {
                        if (par.masterproc)
                        {
                            std::cerr << RERROR "GWRTM: Error reading RTM subsurface boundary polygon file " << rtgbc.id[k] << ": " << fullPathPoly[k] << std::endl;
                            return 0;
                        }
                    }
                }
                if (!rtgbc.rtgwbc[k].find_bcells(gw, rtgbc.id[k], gdom, par, nPoly, xPoly, yPoly))
                    return 0;
            }
            else
            {
                if (par.masterproc)
                {
                    std::cerr << RERROR "GWRTM: RTM Polygon file " << k << ": " << fullPathPoly[k] << " not found." << std::endl;
                    return 0;
                }
            }
            fPoly.close();

        } // end for read in of the k-th polygon

        // --- Read global time-series boundary condition files ---
        for (int k = 0; k < tsFile.size(); k++)
        {
            std::string fname = dir + tsFile[k];
            std::ifstream fts(fname);
            int ndatat = 0, readts = 0;
            // Determine if a time-series file is needed based on boundary type
            switch (rtgbc.rtgwbc[k].rtbctype)
            {
            case SUB_RT_BC_Dirichlet_T:
            case SUB_RT_BC_Neumann_T:
            case SUB_RT_BC_Cauchy_T:
                readts = 1;
                break;
            }
            if (readts)
            {
                if (fts.is_open())
                {
                    fts.ignore(256, ' ');
                    fts >> ndatat;
                    if (ndatat > 0)
                    {
                        rtgbc.rtgwbc[k].ts.initialise(ndatat);
                    }
                    for (int i = 0; i < ndatat; i++)
                    {
                        if (!fts.fail() && !fts.eof())
                        {
                            fts >> rtgbc.rtgwbc[k].ts.time(i) >> rtgbc.rtgwbc[k].ts.value(i);
                        }
                        else
                        {
                            if (par.masterproc)
                            {
                                std::cerr << RERROR "GWRTM: Error reading timeseries file for RTM boundary " << rtgbc.id[k] << ": " << tsFile[k] << std::endl;
                                return 0;
                            }
                        }
                    } // end for ndata
                    fts.close();
                }
                else
                {
                    if (par.masterproc)
                    {
                        std::cerr << RERROR "GWRTM: Error opening RTM timeseries file " << fname << std::endl;
                        return 0;
                    }
                }
            }
        } // end for timeseries files

        // --- Read spatially distributed boundary condition files ---
        for (int k = 0; k < bcFile.size(); k++)
        {
            std::string fname = dir + bcFile[k];
            std::ifstream fbc(fname);
            int nx, ny, nz, readbc = 0;
            // Determine if a spatial BC file is needed based on boundary type
            if (bcFile[k].length() > 0)
            {
                switch (rtgbc.rtgwbc[k].rtbctype)
                {
                case SUB_RT_BC_Neumann_CONST:
                    if (gdom.isRain)
                    {
                        break;
                    }
                    else
                    {
                        readbc = 1;
                    }
                case SUB_RT_BC_Dirichlet_CONST:
                case SUB_RT_BC_Cauchy_CONST:
                    readbc = 1;
                    break;
                }
            }

            if (hasbcfile[k] && readbc && fbc.good())
            {
                // Read the spatial boundary data file
                if (fbc.is_open())
                {
                    std::cout << GOK << "GWRTM: Reading RTM subsurface boundary file : " << fname << std::endl;
                    fbc.ignore(256, ' ');
                    fbc >> nx >> ny >> nz;
                    if (nx * ny * nz != ndata)
                    {
                        std::cerr << RERROR << nx << ny << nz << ndata << std::endl;
                        std::cerr << RERROR "GWRTM: Error reading RTM bc data file " << bcFile[k] << " nx*ny*nz != ndata! " << std::endl;
                        return 0;
                    }
                    if (ndata > 0)
                    {
                        int idx = 0;
                        for (int kk = 0; kk < nz; kk++)
                        {
                            for (int jj = 0; jj < ny; jj++)
                            {
                                for (int ii = 0; ii < nx; ii++)
                                {
                                    if (!fbc.fail() && !fbc.eof())
                                    {
                                        fbc >> rtgbc.rtgwbc[k].bcvals(idx);
                                        idx += 1;
                                    }
                                    else
                                    {
                                        if (par.masterproc)
                                        {
                                            std::cerr << RERROR "GWRTM: Error reading RTM bc file for boundary " << rtgbc.id[k] << ": " << bcFile[k] << std::endl;
                                            return 0;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    fbc.close();
                }
                else
                {
                    if (par.masterproc)
                    {
                        std::cerr << RERROR "GWRTM: Error opening RTM bc file " << fname << std::endl;
                        return 0;
                    }
                }
            }
        }

        // --- Read per-species spatially distributed boundary files ---
        for (int k = 0; k < bccount; k++)
        {
            for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
            {
                // Read spatial boundary data file for this species
                if (rtgbc.rtgwbc[k].has_bcfile[iSpec] && rtgbc.rtgwbc[k].spec_bcfile[iSpec].length() > 0)
                {
                    std::string fname = dir + rtgbc.rtgwbc[k].spec_bcfile[iSpec];
                    int nx, ny, nz, readbc = 0;

                    // Verify the file exists
                    std::ifstream fbc_check(fname);
                    if (!fbc_check.good())
                    {
                        if (par.masterproc)
                        {
                            std::cerr << RERROR "GWRTM: RTM bc file for boundary id=" << rtgbc.id[k] << ", species " << iSpec << " not found: " << fname << std::endl;
                        }
                        return 0;
                    }
                    fbc_check.close();

                    std::ifstream fbc(fname);
                    // Determine whether to read based on boundary type
                    switch (rtgbc.rtgwbc[k].spec_bctype[iSpec])
                    {
                    case SUB_RT_BC_Neumann_CONST:
                        if (gdom.isRain)
                        {
                            break;
                        }
                        else
                        {
                            readbc = 1;
                        }
                    case SUB_RT_BC_Dirichlet_CONST:
                    case SUB_RT_BC_Cauchy_CONST:
                        readbc = 1;
                        break;
                    }

                    if (readbc && fbc.good())
                    {
                        if (fbc.is_open())
                        {
                            if (par.masterproc)
                                std::cout << GOK << "GWRTM: Reading RTM subsurface boundary file for boundary id=" << rtgbc.id[k] << ", species " << iSpec << ": " << fname << std::endl;
                            fbc.ignore(256, ' ');
                            fbc >> nx >> ny >> nz;
                            if (nx * ny * nz != ndata)
                            {
                                std::cerr << RERROR << nx << " " << ny << " " << nz << " " << ndata << std::endl;
                                std::cerr << RERROR "GWRTM: Error reading RTM bc data file for species " << iSpec << ": " << rtgbc.rtgwbc[k].spec_bcfile[iSpec] << " nx*ny*nz != ndata! " << std::endl;
                                return 0;
                            }
                            if (ndata > 0)
                            {
                                int idx = 0;
                                for (int kk = 0; kk < nz; kk++)
                                {
                                    for (int jj = 0; jj < ny; jj++)
                                    {
                                        for (int ii = 0; ii < nx; ii++)
                                        {
                                            if (!fbc.fail() && !fbc.eof())
                                            {
                                                fbc >> rtgbc.rtgwbc[k].spec_bcvals[iSpec](idx);
                                                idx += 1;
                                            }
                                            else
                                            {
                                                if (par.masterproc)
                                                {
                                                    std::cerr << RERROR "GWRTM: Error reading RTM bc file for boundary " << rtgbc.id[k] << ", species " << iSpec << ": " << rtgbc.rtgwbc[k].spec_bcfile[iSpec] << std::endl;
                                                    return 0;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            fbc.close();
                        }
                        else
                        {
                            if (par.masterproc)
                            {
                                std::cerr << RERROR "GWRTM: Error opening RTM bc file for species " << iSpec << ": " << fname << std::endl;
                                return 0;
                            }
                        }
                    }
                }
            }
        }

        // --- Read per-species time-series boundary condition files ---
        for (int k = 0; k < bccount; k++)
        {
            for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
            {
                // Determine if a time-series file is needed for this species
                int readts = 0;
                switch (rtgbc.rtgwbc[k].spec_bctype[iSpec])
                {
                case SUB_RT_BC_Dirichlet_T:
                case SUB_RT_BC_Neumann_T:
                case SUB_RT_BC_Cauchy_T:
                    readts = 1;
                    break;
                }

                // Read the time-series file if a boundary file is specified
                if (readts && rtgbc.rtgwbc[k].spec_bcfile[iSpec].length() > 0)
                {
                    std::string fname = dir + rtgbc.rtgwbc[k].spec_bcfile[iSpec];
                    int ndatat = 0;

                    // Verify the file exists
                    std::ifstream fts_check(fname);
                    if (!fts_check.good())
                    {
                        if (par.masterproc)
                        {
                            std::cerr << RERROR "GWRTM: RTM timeseries file for boundary id=" << rtgbc.id[k] << ", species " << iSpec << " not found: " << fname << std::endl;
                        }
                        return 0;
                    }
                    fts_check.close();

                    std::ifstream fts(fname);
                    if (fts.is_open())
                    {
                        if (par.masterproc)
                            std::cout << GOK << "GWRTM: Reading RTM timeseries file for species " << iSpec << ": " << fname << std::endl;
                        fts.ignore(256, ' ');
                        fts >> ndatat;
                        if (ndatat > 0)
                        {
                            rtgbc.rtgwbc[k].spec_ts[iSpec].initialise(ndatat);
                        }
                        for (int i = 0; i < ndatat; i++)
                        {
                            if (!fts.fail() && !fts.eof())
                            {
                                fts >> rtgbc.rtgwbc[k].spec_ts[iSpec].time(i) >> rtgbc.rtgwbc[k].spec_ts[iSpec].value(i);
                            }
                            else
                            {
                                if (par.masterproc)
                                {
                                    std::cerr << RERROR "GWRTM: Error reading timeseries file for RTM boundary " << rtgbc.id[k] << ", species " << iSpec << ": " << rtgbc.rtgwbc[k].spec_bcfile[iSpec] << std::endl;
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
                            std::cerr << RERROR "GWRTM: Error opening RTM timeseries file for species " << iSpec << ": " << fname << std::endl;
                            return 0;
                        }
                    }
                }
            }
        }

        if (par.masterproc)
            std::cout << GOK << "GWRTM: RTM Subsurface boundary file parsed and boundaries set" << std::endl;
        return 1;
    }

    // =======================================================================
    // setRtState()
    // =======================================================================
    /*!
     * @brief Set initial solute concentrations for all species in both liquid and solid phases.
     *
     * This function initializes the concentration arrays rt.c (liquid/aqueous phase)
     * and rt.c_solid (solid/adsorbed phase) for all grid cells and all species.
     * Two initialization modes are supported per species, controlled by RT_Init_Aq_Mode
     * and RT_Init_Solid_Mode:
     *   - Mode 0: Constant value -- all cells set to the specified constant
     *   - Mode 1: File-based     -- spatial data read from a separate file
     *
     * Solid-phase initialization is only performed when ReactionModule == 1,
     * since adsorbed concentrations are not needed in pure transport mode.
     *
     * @param[in]     inFolder  Input directory path
     * @param[in,out] rt        Reactive transport state (concentrations populated on output)
     * @param[in]     gw        Groundwater state (hydraulic head, used for optional physical corrections)
     * @param[in]     gdom      Groundwater domain (grid dimensions)
     * @param[in]     par       Parallel decomposition
     * @param[in]     io        File I/O handler
     * @return 1 on success; calls exit(-1) on invalid mode
     */
    int setRtState(std::string inFolder, RTStateGW &rt, GwState &gw, GwDomain &gdom, Parallel &par, FileIO &io)
    {
        int ii, jj, kk, idx, ivg, iGlob, iGlobSW;

        std::ifstream fInStream(inFolder + "rttransportgw.input");
        std::string line;
        PsLn pline;
        std::string tempStr;

        // --- Zero-initialize all concentration arrays ---
        for (iGlob = 0; iGlob < gdom.nCellMem; iGlob++)
        {
            for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
            {
                rt.c(iSpec, iGlob, 0) = 0.0;
                rt.c(iSpec, iGlob, 1) = 0.0;
                rt.c_solid(iSpec, iGlob, 0) = 0.0;
                rt.c_solid(iSpec, iGlob, 1) = 0.0;
            }
        }

        // --- Set liquid-phase initial conditions for each species independently ---
        for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
        {
            if (rt.RT_Init_Aq_Mode[iSpec] == 0)
            {
                // Constant mode: set all cells to the specified constant value
                real const_val = rt.RT_Init_Aq_Const[iSpec];
                for (iGlob = 0; iGlob < gdom.nCellMem; iGlob++)
                {
                    rt.c(iSpec, iGlob, 0) = const_val;
                    rt.c(iSpec, iGlob, 1) = const_val;

                    // The following commented block demonstrates an optional physical
                    // correction for a hillslope test case: only setting concentrations
                    // in the saturated zone (where pressure head > 0).
                    // Retained for reference but not active.
                    real pressure_head = gw.h(iGlob, 0);
                }
                if (par.masterproc)
                    std::cerr << GOK "GWRTM: Species " << iSpec << ": aq_mode=0(const), set to aq_val=" << const_val << std::endl;
            }
            else if (rt.RT_Init_Aq_Mode[iSpec] == 1)
            {
                // File mode: read spatial initial conditions from file
                tempStr = rt.RT_Init_Aq_File[iSpec];
                if (par.masterproc)
                    std::cerr << GOK "GWRTM: Species " << iSpec << ": aq_mode=1(file), reading from aq_file=" << tempStr << std::endl;
                readRtICFile(tempStr, inFolder, rt, gdom, par, iSpec);
            }
            else
            {
                std::cerr << RERROR "GWRTM: RT_Init_Aq_" << iSpec << ": mode must be 0(const) or 1(file)" << std::endl;
                exit(-1);
            }
        }

        // --- Set solid-phase initial conditions for each species (only if reaction module is active) ---
        if (rt.ReactionModule == 1)
        {
            for (int iSpec = 0; iSpec < rt.n_mass; iSpec++)
            {
                if (rt.RT_Init_Solid_Mode[iSpec] == 0)
                {
                    // Constant mode: set all cells to the specified constant value
                    real const_val = rt.RT_Init_Solid_Const[iSpec];
                    for (iGlob = 0; iGlob < gdom.nCellMem; iGlob++)
                    {
                        rt.c_solid(iSpec, iGlob, 0) = const_val;
                        rt.c_solid(iSpec, iGlob, 1) = const_val;
                    }
                    if (par.masterproc)
                        std::cerr << GOK "GWRTM: Species " << iSpec << ": solid_mode=0(const), set to solid_val=" << const_val << std::endl;
                }
                else if (rt.RT_Init_Solid_Mode[iSpec] == 1)
                {
                    // File mode: read spatial initial conditions from file
                    tempStr = rt.RT_Init_Solid_File[iSpec];
                    if (par.masterproc)
                        std::cerr << GOK "GWRTM: Species " << iSpec << ": solid_mode=1(file), reading from solid_file=" << tempStr << std::endl;
                    readRt_Solid_ICFile(tempStr, inFolder, rt, gdom, par, iSpec);
                }
                else
                {
                    std::cerr << RERROR "GWRTM: RT_Init_Solid_" << iSpec << ": mode must be 0(const) or 1(file)" << std::endl;
                    exit(-1);
                }
            }
        }

        if (par.masterproc)
        {
            std::cerr << GOK "GWRTM: Subsurface Reactive transport initial condition set" << std::endl;
        }
        return 1;
    }

    // =======================================================================
    // readRtICFile()
    // =======================================================================
    /*!
     * @brief Read aquifer (liquid-phase) initial solute concentration from a spatial file.
     *
     * The file format is:
     *   <header_token> <nx>
     *   <header_token> <ny>
     *   <header_token> <nodata_value>
     *   <data values...>   (nx * ny * nz values, row-major order)
     *
     * The grid dimensions (nx, ny) in the file must match the global domain.
     * Data is mapped from the global grid to the local MPI partition and
     * stored in rt.c(iSpec, :, timeIdx) for both time levels (0 and 1).
     *
     * @param[in]     fNameIn   File name (relative to input directory)
     * @param[in]     fDirIn    Input directory path
     * @param[in,out] rt        Reactive transport state (concentrations populated)
     * @param[in]     gdom      Groundwater domain (grid dimensions)
     * @param[in]     par       Parallel decomposition (for partition offset)
     * @param[in]     iSpec     Species index to read data for (-1 for legacy mode)
     * @return 1 on success, 0 on failure; calls exit(-1) if file cannot be opened
     */
    int readRtICFile(std::string fNameIn, std::string fDirIn, RTStateGW &rt, GwDomain &gdom, Parallel &par, int iSpec = -1)
    {
        std::string fname = fDirIn + fNameIn;
        std::ifstream fInStream(fname);
        std::string line;
        int tnx, tny, iGlob, iGlobSW, idx, ivg, ii, jj, kk, ii2;
        real tmp, nodata_value;
        int ndata = gdom.ny_glob * gdom.nx_glob * gdom.nz_glob;
        intArr tmpVar = intArr("var", ndata);
        std::string str;
        if (fInStream.is_open())
        {
            std::getline(fInStream, str, ' ');
            std::getline(fInStream, str);
            std::stringstream(str) >> tnx;
            std::getline(fInStream, str, ' ');
            std::getline(fInStream, str);
            std::stringstream(str) >> tny;
            std::getline(fInStream, str, ' ');
            std::getline(fInStream, str);
            std::stringstream(str) >> nodata_value;
            // Compare the file dimensions with the DEM to ensure consistency
            if (gdom.ny_glob != tny || gdom.nx_glob != tnx)
            {
                if (par.masterproc)
                {
                    std::cerr << RERROR "GWRTM: RTM IC file parameters don't match DEM parameters. Unable to continue\n";
                    return 0;
                    if (par.masterproc)
                    {
                        std::cerr << BDASH "GWRTM: nx_glob: " << gdom.nx_glob << tnx << "\n";
                        std::cerr << BDASH "GWRTM: ny_glob: " << gdom.ny_glob << tny << "\n";
                    }
                    return 0;
                }
            }
            // Read all data into a temporary view
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
                        std::cerr << RERROR "GWRTM: Error reading RTM IC file. Not enough data\n";
                        return 0;
                    }
                }
            }
            fInStream.close();
        }
        else if (par.masterproc)
        {
            std::cerr << RERROR "GWRTM: Unable to open " << fNameIn << "\n";
            exit(-1);
        }

        // --- Copy data from global temporary array to local partition for the specified species ---
        // Supports arbitrary file names (no longer restricted to hardcoded names)
        for (idx = 0; idx < gdom.nCell; idx++)
        {
            gdom.unpackIndices(idx, kk, jj, ii);
            // Compute the local memory index including halo cells
            iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;

            // Read initial concentration for the specified species
            // Map local indices to the global data array using MPI partition offsets
            ii2 = kk * gdom.nx_glob * gdom.ny_glob + (par.j_beg + jj) * (gdom.nx_glob) + par.i_beg + ii;
            if (ii2 < ndata)
            {
                rt.c(iSpec, iGlob, 1) = tmpVar(ii2);
            }
            else
            {
                rt.c(iSpec, iGlob, 1) = 0.0;
            }
            // Copy from new time level to old time level
            rt.c(iSpec, iGlob, 0) = rt.c(iSpec, iGlob, 1);
        }
        if (par.masterproc)
        {
            std::cerr << GOK "GWRTM: RTM Subsurface Initial Concentration set\n";
        }
        return 1;
    }

    // =======================================================================
    // readRt_Solid_ICFile()
    // =======================================================================
    /*!
     * @brief Read solid-phase initial solute concentration from a spatial file.
     *
     * Supports two modes:
     *   1. Per-species mode (iSpec >= 0): Reads solid-phase concentration for a single
     *      species from the file.  The file contains nx * ny * nz values.
     *   2. Legacy multi-species mode (iSpec == -1, file named "concen_solid.input"):
     *      Reads all species from a single file.  For multi-species (n_mass > 1),
     *      the file layout is expected to be nx * ny * nz * n_mass values with
     *      species interleaved.
     *
     * The file format is the same as readRtICFile:
     *   <header_token> <nx>
     *   <header_token> <ny>
     *   <header_token> <nodata_value>
     *   <data values...>
     *
     * @param[in]     fNameIn   File name (relative to input directory)
     * @param[in]     fDirIn    Input directory path
     * @param[in,out] rt        Reactive transport state (solid-phase concentrations populated)
     * @param[in]     gdom      Groundwater domain (grid dimensions)
     * @param[in]     par       Parallel decomposition (for partition offset)
     * @param[in]     iSpec     Species index (-1 triggers legacy multi-species mode)
     * @return 1 on success, 0 on failure; calls exit(-1) if file cannot be opened
     */
    int readRt_Solid_ICFile(std::string fNameIn, std::string fDirIn, RTStateGW &rt, GwDomain &gdom, Parallel &par, int iSpec = -1)
    {
        std::string fname = fDirIn + fNameIn;
        std::ifstream fInStream(fname);
        std::string line;
        int tnx, tny, iGlob, iGlobSW, idx, ivg, ii, jj, kk, ii2;
        real tmp, nodata_value;
        int ndata = gdom.ny_glob * gdom.nx_glob * gdom.nz_glob;
        intArr tmpVar = intArr("var", ndata);
        std::string str;
        if (fInStream.is_open())
        {
            std::getline(fInStream, str, ' ');
            std::getline(fInStream, str);
            std::stringstream(str) >> tnx;
            std::getline(fInStream, str, ' ');
            std::getline(fInStream, str);
            std::stringstream(str) >> tny;
            std::getline(fInStream, str, ' ');
            std::getline(fInStream, str);
            std::stringstream(str) >> nodata_value;
            // Compare the file dimensions with the DEM to ensure consistency
            if (gdom.ny_glob != tny || gdom.nx_glob != tnx)
            {
                if (par.masterproc)
                {
                    std::cerr << RERROR "GWRTM: RTM rt_Solid_IC file parameters don't match DEM parameters. Unable to continue\n";
                    if (par.masterproc)
                    {
                        std::cerr << BDASH "GWRTM: nx_glob: " << gdom.nx_glob << tnx << "\n";
                        std::cerr << BDASH "GWRTM: ny_glob: " << gdom.ny_glob << tny << "\n";
                    }
                    return 0;
                }
            }
            // Read all data into a temporary view
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
                        std::cerr << RERROR "GWRTM: Error reading rt_Solid_IC file. Not enough data\n";
                        return 0;
                    }
                }
            }
            fInStream.close();
        }
        else
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "GWRTM: Unable to open " << fNameIn << "\n";
                exit(-1);
            }
        }

        // --- Copy data from global temporary array to local partition ---
        // Supports arbitrary file names (no longer restricted to hardcoded names)
        for (idx = 0; idx < gdom.nCell; idx++)
        {
            gdom.unpackIndices(idx, kk, jj, ii);
            // Compute the local memory index including halo cells
            iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;

            if (iSpec >= 0)
            {
                // Per-species mode: read solid-phase initial condition for a single species
                // Map local indices to the global data array using MPI partition offsets
                ii2 = kk * gdom.nx_glob * gdom.ny_glob + (par.j_beg + jj) * (gdom.nx_glob) + par.i_beg + ii;
                if (ii2 < ndata)
                {
                    rt.c_solid(iSpec, iGlob, 1) = tmpVar(ii2);
                }
                else
                {
                    rt.c_solid(iSpec, iGlob, 1) = 0.0;
                }
                // Copy from new time level to old time level
                rt.c_solid(iSpec, iGlob, 0) = rt.c_solid(iSpec, iGlob, 1);
            }
            else
            {
                // Legacy backward-compatible mode: read for all species at once
                // Only triggered when the file name is exactly "concen_solid.input"
                if (!strcmp(fNameIn.c_str(), "concen_solid.input"))
                {
                    // Read solid-phase concentration for each species
                    for (int iSpec_local = 0; iSpec_local < rt.n_mass; iSpec_local++)
                    {
                        // Map local indices to the global data array
                        ii2 = kk * gdom.nx_glob * gdom.ny_glob + (par.j_beg + jj) * (gdom.nx_glob) + par.i_beg + ii;

                        if (rt.n_mass > 1)
                        {
                            // Multi-species: species values are interleaved in the file
                            int file_idx = ii2 * rt.n_mass + iSpec_local;
                            if (file_idx < ndata)
                            {
                                rt.c_solid(iSpec_local, iGlob, 1) = tmpVar(file_idx);
                            }
                            else
                            {
                                rt.c_solid(iSpec_local, iGlob, 1) = 0.0;
                            }
                        }
                        else
                        {
                            // Single species: direct mapping
                            rt.c_solid(iSpec_local, iGlob, 1) = tmpVar(ii2);
                        }

                        // Copy from new time level to old time level
                        rt.c_solid(iSpec_local, iGlob, 0) = rt.c_solid(iSpec_local, iGlob, 1);
                    }
                }
                else
                {
                    if (par.masterproc)
                    {
                        std::cerr << RERROR "GWRTM: Error reading RTM Solid concentration IC file. File name might be wrong.\n";
                        return 0;
                    }
                }
            }
        }
        if (par.masterproc)
        {
            std::cerr << GOK "GWRTM: RTM Soils Subsurface concentration content set\n";
        }
        return 1;
    }

}; // end class RTInitgw

#endif
#endif
