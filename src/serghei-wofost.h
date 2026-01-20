/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _SERGHEI_H_
#define _SERGHEI_H_

#include "stdlib.h"
#include <iostream>
#include <string>
#include "define.h"
#include "Domain.h"
#include "Parallel.h"
#include "Parser.h"
#include "State.h"
#include "BC.h"
#include "Initializer.h"
#include "TimeIntegrator.h"
#include "FileIO.h"
#include "Exchange.h"
#include "SourceSink.h"
#include "DomainIntegrator.h"
#include "Vegetation.h"
#include "ParticleTracking.h"
#include "tools.h"
#include <thread> 

#include <cstdio>
#include <cstdlib>

#if SERGHEI_SUBSURFACE_MODEL
#include "GwDomain.h"
#include "GwFunction.h"
#include "GwInit.h"
#include "GwMPI.h"
#include "GwMatrix.h"
#include "GwState.h"
#include "GwSolver.h"
#include "GwIntegrator.h"
#endif

// --- WOFOST Integration ---
#if CROP_GROWTH_MODEL
#include "Wofost72.h"
#include "CropState.h"
#include "CropInit.h"
#include "MeteoState.h"
#include "MeteoInit.h"
// Helper function to compute root zone moisture from SERGHEI state
// This function needs to be implemented (likely in a new utility file or GwFunction)
// void compute_root_zone_moisture(const GwState& gw, const GwDomain& gdom, 
//                                 const realArr& RD, realArr& root_SM);
// Assuming it exists for now.
#endif
// --------------------------

// ... (Other includes for Transport modules) ...

class SERGHEI
{
public:
	Parallel par;
	Initializer init;
	FileIO io;

	State state;
	Domain dom;
#if SERGHEI_SUBSURFACE_MODEL
	GwState gw;
	GwDomain gdom;
	GwInit ginit;
	GwMPI gmpi;
	GwFunction gwf;
	SubsurfaceBoundaries gbc;
	GwMatrix A;
	GwIntegrator gint;
#ifdef __NVCC__
	GwSolver<Kokkos::Cuda> gsolver;
#else
	GwSolver<Kokkos::OpenMP> gsolver;
#endif
#endif

// ... (Transport modules) ...

	SourceSink ss;
	Exchange exch;
	ExternalBoundaries ebc;

    // --- WOFOST Objects ---
#if CROP_GROWTH_MODEL
    Wofost72 wofost;
    CropState cropParam;
    MeteoState meteo;
    CropInit cropInit;
    MeteoInit meteoInit;
    
    // Intermediate variable for coupling
    realArr current_root_zone_SM; 
    
    // Timer for crop step
    double time_since_last_crop_step = 0.0;
    const double CROP_STEP_SIZE = 86400.0; // 1 Day
#endif
    // ----------------------

private:
	Parser parser;
	TimeIntegrator tint;
	surfaceIntegrator sint;
	boundaryIntegrator bint;
// ... (Observations, etc) ...

	double oldVolume, newVolume, diffVolume;
	double accumDt = 0.0;

	// Kokkos objects
	Kokkos::Timer timer;
	Kokkos::InitializationSettings kokkosSettings;

public:
	std::string inFolder, outFolder;

	////////////// METHODS ///////////////
public:
	int start(int argc, char **argv)
	{
        // ... (Standard SERGHEI initialization code: MPI, Kokkos, Domain) ...
        // (Assuming standard initialization code remains here)
        
        init.initializeMPI(&argc, &argv, par);
        // ... set settings ...
        Kokkos::initialize(kokkosSettings);
        // ... init model ...
        if (!init.initialize(state, ss.swss, ebc, dom, par, tint, sint, bint, parser, exch, io, inFolder, outFolder)) return 0;

#if SERGHEI_SUBSURFACE_MODEL
		if (!ginit.initialize_gw(gw, gdom, state, dom, gbc, gmpi, gint, par, io, ss, inFolder, outFolder)) return 0;
		A.init(gdom);
		gsolver.init(A, gdom);
#endif

// ... (Transport Init) ...

        // --- WOFOST Initialization ---
#if CROP_GROWTH_MODEL
        if (par.masterproc) std::cout << GOK "Initializing WOFOST Crop Model..." << std::endl;

        // 1. Read Crop Parameters
        if (!cropInit.readCropParameters(inFolder + "crop.input", cropParam, par)) {
            std::cerr << RERROR "Failed to read crop.input" << std::endl;
            return 0;
        }

        // 2. Read Meteo Data
        if (!meteoInit.readMeteoFile(inFolder + "meteo.input", meteo, par)) {
            std::cerr << RERROR "Failed to read meteo.input" << std::endl;
            return 0;
        }

        // 3. Allocate Wofost Memory
        // Wofost runs on 2D surface cells (assuming crop grows on surface)
        // Or 3D cells if we track root depth per column? 
        // Typically crop models are 1D (vertical) per surface grid cell.
        // SERGHEI Surface domain size: dom.nCellMem (includes halo) or dom.nCell (physical)
        // Let's use gdom.nCellMem (Subsurface size) or dom.nCellMem (Surface size).
        // Since root growth interacts with subsurface columns, let's assume one crop per SURFACE cell.
        // We need to map surface index to subsurface column index.
        // For simplicity, let's assume we use the surface grid size: dom.nCellMem
        
        int n_crop_cells = dom.nCellMem; 
        
        // However, Wofost72 constructor expects an int to allocate Views.
        // We need to construct it here. 
        // Note: Wofost72 wofost; declared in class is default constructed.
        // We need an init/allocate method in Wofost72 or use placement new / pointer.
        // Let's assume Wofost72 has an `allocate(int n)` method or we assign it.
        // BETTER: Make wofost a pointer `Wofost72* wofost;` in class definition
        // wofost = new Wofost72(n_crop_cells);
        // OR: Add `wofost.allocate(n_crop_cells)` to Wofost72 class. 
        // Let's assume the latter for stack allocation safety.
        wofost.allocate(n_crop_cells); // Need to add this to Wofost72.h!

        current_root_zone_SM = realArr("RootZoneSM", n_crop_cells);

        // 4. Compute Initial Soil Moisture for WOFOST
        // compute_root_zone_moisture(gw, gdom, wofost.s.RD, current_root_zone_SM); 
        // Placeholder: Initialize with Field Capacity or from gw state
        Kokkos::deep_copy(current_root_zone_SM, 0.3); // Dummy init

        // 5. Initialize WOFOST State
        int start_doy = (int)(dom.startTime / 86400.0) % 365 + 1;
        wofost.initialize(cropParam, meteo, current_root_zone_SM, start_doy);

        if (par.masterproc) std::cout << GOK "WOFOST Initialized." << std::endl;
#endif
        // -----------------------------

        // ... (Rest of start function: Observations, Initial Output) ...
        
        return 1;
    }

	int compute()
	{
        // ... (Setup dt) ...
#if SERGHEI_SUBSURFACE_MODEL
		gdom.dt = gdom.dt_init;
#endif

		// Main Time Loop
		while (dom.etime < dom.endTime)
		{
            // ... (Surface Integrator, etc.) ...

			// run subsurface model
#if SERGHEI_SUBSURFACE_MODEL
            // ... (Rainfall copy) ...

            // --- WOFOST Coupling Step ---
#if CROP_GROWTH_MODEL
            // Update the crop timer
            time_since_last_crop_step += dom.dt; // Use current step dt

            // Check if it's time to run WOFOST (Daily)
            if (time_since_last_crop_step >= CROP_STEP_SIZE) {
                
                if (par.masterproc) std::cout << GOK "Running WOFOST Step at t=" << dom.etime << std::endl;

                // 1. Calculate average root zone moisture from SERGHEI state
                // This function needs to map 3D gw.wc to 2D surface crop cells based on RD
                // compute_root_zone_moisture(gw, gdom, wofost.s.RD, current_root_zone_SM);
                // Placeholder:
                Kokkos::deep_copy(current_root_zone_SM, 0.3); // TODO: Implement Kernel

                // 2. Prepare time variables
                int current_doy = (int)(dom.etime / 86400.0) % 365 + 1;

                // 3. Run WOFOST
                wofost.calc_rates_and_integrate(cropParam, meteo, dom.etime, 
                                              current_root_zone_SM, current_doy, 1.0);

                // 4. Feedback to SERGHEI
                // Update parameters in SourceSink (ss) or Evaporation BC based on new LAI/RD
                // Example: Update transpiration coefficient based on LAI
                // update_serghei_parameters(ss, wofost.s.LAI, wofost.s.RD);
                
                // Reset timer
                time_since_last_crop_step = 0.0;
            }
#endif
            // ----------------------------

            // ... (Standard SERGHEI GW Solver Loop) ...
            // Asynchronous coupling logic...
            // Solver calls...
#endif

            // ... (Volume check, Output writing) ...
            
            // --- Output WOFOST Variables ---
#if CROP_GROWTH_MODEL
            // If output interval reached
			if (fabs(dom.etime - dom.startTime - io.numOut * io.outFreq) <= TOL12) {
                // You might want to write a specific crop output file
                // io.writeCropOutput(wofost, dom, par, outFolder);
            }
#endif
            
            // ... (Particle tracking, etc) ...

		} // end of time loop

		return 1;
	}

	int finalise()
	{
        // ... (Timers, IO close) ...
#if CROP_GROWTH_MODEL
        wofost.finalize(); // Calculate HI
        // Write final crop summary
#endif
		return 1;
	}
};

#endif