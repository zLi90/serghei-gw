/* -*- mode: c++; c-default-style: "linux" -*- */

/**
 * @file RTStateSW.h
 * @brief State variables for surface water reactive (solute) transport.
 *
 * Defines the RTStateSW class which holds all concentration fields,
 * transport parameters, reaction parameters, and initial condition
 * configuration for multi-species solute transport in surface water.
 * Supports the nitrogen cycle (NH4, NO3, DO, DOC, DON) via
 * multi-component Monod kinetics.
 */

#ifndef _RT_STATE_SW_H_
#define _RT_STATE_SW_H_

#if SERGHEI_SURFACE_TRANSPORT

#include "const.h"
#include "define.h"
#include "Domain.h"
#include "SArray.h"
#include "Indexing.h"
#include "RTStateGW.h"
#include <set>
#include <math.h>

class RTStateSW
{

public:
	// =====================================================================
	// Time stepping
	// =====================================================================
	real dt;  ///< Transport sub-timestep size [s]

	// =====================================================================
	// Primary state variables (multi-species, Kokkos arrays)
	// =====================================================================
	realArr3 c;           ///< Liquid-phase concentration [mg/L or kg/m³]. Indexed as c(iSpec, iGlob, iTime) where iTime: 0=previous, 1=current
	realArr3 solute_mass; ///< Solute mass per cell [mg or kg]. Indexed as solute_mass(iSpec, iGlob, iTime), used for mass-conservative updates
	realArr2 csw4gw;      ///< Surface water concentration copied to groundwater as a boundary condition. Indexed as csw4gw(iSpec, iGlob) [mg/L]
	realArr2 ConQss;      ///< Surface-subsurface solute exchange flux q*c. Indexed as ConQss(iSpec, iGlob) [mg/(m²·s) or kg/(m²·s)]

	// =====================================================================
	// Hydraulic / grid variables
	// =====================================================================
	realArr Pe;  ///< Cell Peclet number [-], ratio of advective to diffusive transport rate
	realArr h;   ///< Water depth [m], copied from State for transport calculations

	// =====================================================================
	// Configuration flags
	// =====================================================================
	int RtInitialModeSW;  ///< Initial condition mode for reactive transport (deprecated, kept for backward compatibility)
	int InitialCon;       ///< Initial concentration flag (deprecated)
	int n_mass;           ///< Number of solute species (chemical components) [-]
	int Advection_Scheme; ///< Advection scheme selector: 1=1st-order upwind, 2=TVD Van Leer limiter, 3=TVD Superbee limiter
	int east, west, north, south; ///< Index offsets for the four neighboring cells
	real c0;  ///< Initial concentration (deprecated, kept for backward compatibility) [mg/L]

	// =====================================================================
	// Transport parameters (per-species, allocated after n_mass is known)
	// =====================================================================
	realArr diffusion_molecular; ///< Molecular diffusion coefficient per species [m²/s]. Typical range: 1e-10 to 1e-8
	realArr alpha_L;             ///< Longitudinal dispersivity per species [m]. Typical range: 0.01 to 100
	realArr alpha_T;             ///< Transverse dispersivity per species [m]. Typical range: 0.001 to 10

	/// Temporary storage for transport parameters before Kokkos allocation
	std::vector<real> temp_diffusion_molecular;  ///< Pre-allocation staging [m²/s]
	std::vector<real> temp_alpha_L;              ///< Pre-allocation staging [m]
	std::vector<real> temp_alpha_T;              ///< Pre-allocation staging [m]

	// =====================================================================
	// Reaction module configuration
	// =====================================================================
	real ReactionModule;            ///< Reaction module switch: 0=no reactions, 1=reactions enabled [-]
	real ReactModel_Type;           ///< Reaction model type: 0=no reaction, 1=multi-component Monod kinetics [-]
	real Nitrogen_Cycle_Simulation; ///< Nitrogen cycle simulation switch: 0=off, 1=on [-]

	// =====================================================================
	// Decay reaction parameters (per-species)
	// =====================================================================
	realArr Lambda_1; ///< First-order decay rate in liquid phase per species [1/s]. Typical range: 0 to 1e-5
	realArr Lambda_2; ///< First-order decay rate in solid phase per species [1/s]. Typical range: 0 to 1e-5

	// =====================================================================
	// Nitrification parameters (single values, Monod kinetics)
	// NH4+ + 2O2 -> NO3- + 2H+ + H2O
	// =====================================================================
	real Rate_Max_Nit;  ///< Maximum nitrification rate [mg/L/s]. Typical: 1e-6 to 1e-4
	real K_Monod_NH4;   ///< Half-saturation constant for NH4 in nitrification [mg/L]. Typical: 0.5 to 5.0
	real K_Monod_DO;    ///< Half-saturation constant for dissolved oxygen in nitrification [mg/L]. Typical: 0.1 to 2.0

	// =====================================================================
	// Denitrification parameters (single values, Monod kinetics with DO inhibition)
	// NO3- -> NO2- -> N2
	// =====================================================================
	real Rate_Max_Denit; ///< Maximum denitrification rate [mg/L/s]. Typical: 1e-6 to 1e-4
	real K_Monod_NO3;    ///< Half-saturation constant for NO3 in denitrification [mg/L]. Typical: 0.1 to 5.0
	real K_Monod_DOC;    ///< Half-saturation constant for DOC (carbon source) in denitrification [mg C/L]. Typical: 0.5 to 10.0
	real Ki_Inhib_DO;    ///< DO inhibition constant for denitrification [mg/L]. Typical: 0.1 to 2.0. Higher DO inhibits denitrification

	// =====================================================================
	// Aerobic respiration parameters (DOC consumption with O2)
	// =====================================================================
	real Rate_Max_Hetero;      ///< Maximum aerobic respiration rate [mg/L/s]. Typical: 1e-6 to 1e-4
	real K_Monod_DOC_Aerobic;  ///< Half-saturation constant for DOC in aerobic respiration [mg C/L]. Typical: 0.5 to 10.0
	real K_Monod_DO_Aerobic;   ///< Half-saturation constant for O2 in aerobic respiration [mg/L]. Typical: 0.1 to 2.0

	// =====================================================================
	// Mineralization parameters (DON -> NH4)
	// =====================================================================
	real Rate_Max_Min;    ///< Maximum mineralization rate [mg/L/s]. Typical: 1e-7 to 1e-5
	real K_Monod_DON;     ///< Half-saturation constant for DON in mineralization [mg N/L]. Typical: 0.1 to 5.0
	real K_Monod_DO_Min;  ///< Half-saturation constant for DO in mineralization [mg/L]. Typical: 0.1 to 2.0

	// =====================================================================
	// Surface-water specific source/sink parameters
	// =====================================================================
	real Ka_sw;      ///< Surface water reaeration rate constant [1/s]. Typical: 1e-6 to 1e-4
	real DO_sat;     ///< Saturated dissolved oxygen concentration [mg/L]. Typical: 8-12 at 15-25°C
	real DOC_eq_sw;  ///< Equilibrium DOC concentration in surface water [mg C/L]. Typical: 2-20
	real DON_eq_sw;  ///< Equilibrium DON concentration in surface water [mg N/L]. Typical: 0.1-5.0
	real K_rel_sw;   ///< Solute relative release rate constant in surface water [1/s]. Typical: 1e-7 to 1e-5

	// =====================================================================
	// Environmental correction factors (single values)
	// =====================================================================
	real Temp_Coeff_Theta; ///< Temperature correction coefficient (theta) [-]. Typical: 1.02-1.10
	real Opt_Temp;         ///< Optimal (reference) temperature [°C]. Typical: 20-25
	real Opt_pH;           ///< Optimal pH for biological reactions [-]. Typical: 6.5-8.5

	// =====================================================================
	// Initial condition configuration (per-species)
	// =====================================================================
	/// Liquid-phase initial conditions: each species can have an independent mode and value
	std::vector<int> RT_Init_Aq_Mode;          ///< Initial mode per species: 0=constant value, 1=read from file [-]
	std::vector<real> RT_Init_Aq_Const;        ///< Constant initial concentration per species [mg/L]
	std::vector<std::string> RT_Init_Aq_File;  ///< Initial condition filename per species

	// =====================================================================
	// Rainfall concentration parameters (per-species)
	// =====================================================================
	std::vector<real> RainCon;   ///< Rainfall solute concentration per species (host vector) [mg/L]
	realArr RainCon_arr;         ///< Rainfall solute concentration per species (Kokkos device array) [mg/L]

	// =====================================================================
	// Allocate state variables for surface water transport
	// Must be called after n_mass is set
	// =====================================================================
	inline void
	allocate(Domain &dom)
	{

		Pe = realArr("Pe", dom.nCellMem);
		/// solute_mass(iSpec, iGlob, iTime): iTime=0 previous, iTime=1 current
		solute_mass = realArr3("solute_mass", n_mass, dom.nCellMem, 2);

		c = realArr3("c", n_mass, dom.nCellMem, 2);

		h = realArr("h", dom.nCellMem);

		csw4gw = realArr2("csw4gw", n_mass, dom.nCellMem);
		ConQss = realArr2("ConQss", n_mass, dom.nCell);

		// Allocate per-species transport parameter arrays
		diffusion_molecular = realArr("diffusion_molecular", n_mass);
		alpha_L = realArr("alpha_L", n_mass);
		alpha_T = realArr("alpha_T", n_mass);

		// Allocate rainfall concentration array
		RainCon_arr = realArr("RainCon_arr", n_mass);

		// Allocate per-species reaction parameter arrays
		Lambda_1 = realArr("Lambda_1", n_mass);
		Lambda_2 = realArr("Lambda_2", n_mass);
	}
};

#endif

#endif
