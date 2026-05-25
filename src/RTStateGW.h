/* -*- mode: c++; c-default-style: "linux" -*- */

/**
 * @file RTStateGW.h
 * @brief State variables for groundwater reactive (solute) transport.
 *
 * Defines the RTStateGW class which holds all concentration fields,
 * dispersion tensors, sorption parameters, reaction parameters (Monod
 * nitrogen cycle kinetics), and initial condition configuration for
 * multi-species solute transport in groundwater (subsurface).
 *
 * Supports: multiple sorption models (Linear, Freundlich, Langmuir),
 * nitrogen cycle simulation (nitrification, denitrification, aerobic
 * respiration, mineralization), and surface-subsurface solute exchange.
 */

#ifndef _RT_STATEGW_H_
#define _RT_STATEGW_H_

#if SERGHEI_SUBSURFACE_TRANSPORT

#include "const.h"
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "GwDomain.h"
#include <set>
#include <math.h>

class RTStateGW
{

public:
	// =====================================================================
	// Time stepping
	// =====================================================================
	real dt_init, dt_max, dt, dtOld;  ///< Transport time step: initial, maximum, current, previous [s]

	// =====================================================================
	// Upstream weighting coefficients for advection discretization
	// =====================================================================
	real Up_Weighting_vplus;   ///< Upstream weighting factor for positive velocity direction [-]
	real Up_Weighting_vminus;  ///< Upstream weighting factor for negative velocity direction [-]
	realArr Up_Weighting_xp, Up_Weighting_xm,  ///< Per-cell upstream weights in +x and -x directions [-]
	        Up_Weighting_yp, Up_Weighting_ym,  ///< Per-cell upstream weights in +y and -y directions [-]
	        Up_Weighting_zp, Up_Weighting_zm;  ///< Per-cell upstream weights in +z and -z directions [-]

	// =====================================================================
	// Species count
	// =====================================================================
	int n_mass;  ///< Number of solute species (chemical components) [-]

	// =====================================================================
	// Transport parameters
	// =====================================================================
	realArr diffusion_molecular;  ///< Molecular diffusion coefficient per species [m²/s]. Typical: 1e-10 to 1e-8
	realArr alpha_T;              ///< Transverse dispersivity per species [m]. Typical: 0.001 to 10
	realArr alpha_L;              ///< Longitudinal dispersivity per species [m]. Typical: 0.01 to 100
	realArr2 q_dispersion;        ///< Interface flow velocity for dispersion tensor calculation [m/s]. Indexed as q_dispersion(iGlob, 12)
	real Pe;                      ///< Global Peclet number [-] (deprecated scalar, kept for compatibility)

	// =====================================================================
	// Water content for reaction term computation
	// =====================================================================
	realArr wc_mid,          ///< Mid-step water content used for reaction terms [-]. Range: 0 to 1
	       wc_original_old,   ///< Water content at previous time step [-]
	       wc_original_new;   ///< Water content at current time step [-]

	// =====================================================================
	// Reaction module switches
	// =====================================================================
	real ReactionModule;            ///< Reaction module switch: 0=no reactions, 1=reactions enabled [-]
	real rho_b;                     ///< Bulk density of porous medium [kg/m³]. Typical: 1300-1800

	real ReactModel_Type;           ///< Reaction model type: 0=no reaction, 1=multi-component Monod kinetics [-]
	realArr Sorption_Type;          ///< Sorption type per species: 0=None, 1=Linear, 2=Freundlich, 3=Langmuir [-]
	real Nitrogen_Cycle_Simulation; ///< Nitrogen cycle simulation switch: 0=off, 1=on [-]

	// =====================================================================
	// Single-species parameters (backward compatibility)
	// =====================================================================
	real Kl;     ///< Langmuir equilibrium constant (single-species legacy) [L/mg]
	real eta;    ///< Langmuir sorption capacity (single-species legacy) [mg/kg]
	real lambda; ///< Decay coefficient (single-species legacy) [1/s]

	// =====================================================================
	// Multi-species sorption parameters
	// =====================================================================
	realArr Kd;       ///< Distribution coefficient (Linear isotherm) per species [L/kg or m³/kg]. Typical: 0 to 100
	realArr Kf;       ///< Freundlich constant per species [(mg/kg)/(mg/L)^Nf]. Typical: 0.01 to 10
	realArr Nf;       ///< Freundlich nonlinearity exponent per species [-]. Typical: 0.5 to 1.2
	realArr Alpha_D;  ///< Langmuir equilibrium constant per species [L/mg]. Typical: 0.01 to 10
	realArr Beta_D;   ///< Langmuir sorption capacity per species [mg/kg]. Typical: 100 to 10000
	realArr Beta;     ///< Non-equilibrium sorption parameter per species [1/s]. Typical: 0 to 1e-4
	realArr Lambda_1; ///< First-order liquid-phase decay rate per species [1/s]. Typical: 0 to 1e-5
	realArr Lambda_2; ///< First-order solid-phase decay rate per species [1/s]. Typical: 0 to 1e-5

	// =====================================================================
	// Decay reaction parameters (per-species array)
	// =====================================================================
	realArr Decay_Rate;  ///< Decay rate per species [1/s]. General-purpose array for decay

	// =====================================================================
	// Nitrification parameters (single values, Monod kinetics)
	// NH4+ + 2O2 -> NO3- + 2H+ + H2O
	// =====================================================================
	real Rate_Max_Nit;  ///< Maximum nitrification rate [mg/L/s]. Typical: 1e-6 to 1e-4
	real K_Monod_NH4;   ///< Half-saturation constant for NH4 [mg/L]. Typical: 0.5 to 5.0
	real K_Monod_DO;    ///< Half-saturation constant for dissolved oxygen [mg/L]. Typical: 0.1 to 2.0
	real Y_O2_NH4;      ///< Stoichiometric yield coefficient: O2 consumed per NH4 [mg/mg]

	// =====================================================================
	// Denitrification parameters (single values, Monod with DO inhibition)
	// NO3- -> NO2- -> N2
	// =====================================================================
	real Rate_Max_Denit; ///< Maximum denitrification rate [mg/L/s]. Typical: 1e-6 to 1e-4
	real K_Monod_NO3;    ///< Half-saturation constant for NO3 [mg/L]. Typical: 0.1 to 5.0
	real K_Monod_DOC;    ///< Half-saturation constant for DOC [mg C/L]. Typical: 0.5 to 10.0
	real Ki_Inhib_DO;    ///< DO inhibition constant [mg/L]. Typical: 0.1 to 2.0

	// =====================================================================
	// Aerobic respiration parameters (DOC + O2 -> CO2)
	// =====================================================================
	real Rate_Max_Hetero;      ///< Maximum aerobic respiration rate [mg/L/s]. Typical: 1e-6 to 1e-4
	real K_Monod_DOC_Aerobic;  ///< Half-saturation constant for DOC [mg C/L]. Typical: 0.5 to 10.0
	real K_Monod_DO_Aerobic;   ///< Half-saturation constant for O2 [mg/L]. Typical: 0.1 to 2.0

	// =====================================================================
	// Mineralization parameters (DON -> NH4)
	// =====================================================================
	real Rate_Max_Min;    ///< Maximum mineralization rate [mg/L/s]. Typical: 1e-7 to 1e-5
	real K_Monod_DON;     ///< Half-saturation constant for DON [mg N/L]. Typical: 0.1 to 5.0
	real K_Monod_DO_Min;  ///< Half-saturation constant for DO [mg/L]. Typical: 0.1 to 2.0

	// =====================================================================
	// Reaeration, DOC/DON release parameters (groundwater-specific)
	// =====================================================================
	real Ka_DO;      ///< Reaeration rate constant [1/s]. Typical: 1e-7 to 1e-5
	real DO_sat;     ///< Saturated dissolved oxygen concentration [mg/L]. Typical: 8-12
	real DOC_eq;     ///< DOC release equilibrium concentration [mg C/L]. Typical: 2-20
	real DON_eq;     ///< DON release equilibrium concentration [mg N/L]. Typical: 0.1-5.0
	real K_rel_DOC;  ///< DOC release rate constant [1/s]. Typical: 1e-7 to 1e-5
	real K_rel_DON;  ///< DON release rate constant [1/s]. Typical: 1e-7 to 1e-5

	// =====================================================================
	// Environmental correction factors
	// =====================================================================
	real Temp_Coeff_Theta; ///< Temperature correction coefficient (theta) [-]. Typical: 1.02-1.10
	real Opt_Temp;         ///< Optimal (reference) temperature [°C]. Typical: 20-25
	real Opt_pH;           ///< Optimal pH for biological reactions [-]. Typical: 6.5-8.5

	realArr Rf;  ///< Retardation factor per cell [-]. Rf = 1 + rho_b * Kd / theta. Range: >= 1

	// =====================================================================
	// Primary concentration fields (multi-species)
	// =====================================================================
	/// Liquid-phase concentration [mg/L]. Indexed as c(iSpec, iGlob, iTime) where iTime: 0=previous, 1=current
	realArr3 c;
	/// Surface water solute concentration for boundary coupling [mg/L]. Indexed as csw(iSpec, iGlob)
	realArr2 csw;
	/// Surface-subsurface solute exchange flux q*c [mg/(m²·s)]. Indexed as ConQss(iSpec, iGlob)
	realArr2 ConQss;
	/// Rainfall solute concentration (Kokkos device array) [mg/L]. Copied from surface water module
	realArr RainCon_arr;

	/// Solid-phase (adsorbed) concentration [mg/kg]. Indexed as c_solid(iSpec, iGlob, iTime)
	realArr3 c_solid;

	// =====================================================================
	// Index access helper functions (Host-side)
	// =====================================================================
	/// Get liquid-phase concentration value for a given cell, species, and time index
	inline real get_c(int iGlob, int iSpec, int iTime) const
	{
		return c(iGlob, iSpec, iTime);
	}

	/// Set liquid-phase concentration value for a given cell, species, and time index
	inline void set_c(int iGlob, int iSpec, int iTime, real value)
	{
		c(iGlob, iSpec, iTime) = value;
	}

	/// Get solid-phase concentration value for a given cell, species, and time index
	inline real get_c_solid(int iGlob, int iSpec, int iTime) const
	{
		return c_solid(iGlob, iSpec, iTime);
	}

	/// Set solid-phase concentration value for a given cell, species, and time index
	inline void set_c_solid(int iGlob, int iSpec, int iTime, real value)
	{
		c_solid(iGlob, iSpec, iTime) = value;
	}

	// =====================================================================
	// Dispersion tensor (multi-species)
	// =====================================================================
	/// Dispersion tensor components per species. Indexed as dcal(iSpec, iGlob, iComponent)
	/// iComponent maps: 0=d_xx, 1=d_yy, 2=d_zz, 3=d_xy, 4=d_xz, 5=d_yx, 6=d_yz, 7=d_zx, 8=d_zy [m²/s]
	realArr3 dcal;

	// =====================================================================
	// Initial condition mode (legacy, per-group)
	// =====================================================================
	int RT_Aquifer_initialMode,  ///< Aquifer initial concentration mode (deprecated) [-]
	    RT_Solid_initialMode,    ///< Solid-phase initial concentration mode (deprecated) [-]
	    rt_scheme;               ///< Reactive transport numerical scheme selector [-]

	// =====================================================================
	// Per-species initial condition configuration
	// =====================================================================
	/// Liquid-phase initial conditions per species
	std::vector<int> RT_Init_Aq_Mode;          ///< Mode per species: 0=constant, 1=file [-]
	std::vector<real> RT_Init_Aq_Const;        ///< Constant initial value per species [mg/L]
	std::vector<std::string> RT_Init_Aq_File;  ///< IC filename per species

	/// Solid-phase initial conditions per species
	std::vector<int> RT_Init_Solid_Mode;          ///< Mode per species: 0=constant, 1=file [-]
	std::vector<real> RT_Init_Solid_Const;        ///< Constant initial value per species [mg/kg]
	std::vector<std::string> RT_Init_Solid_File;  ///< IC filename per species

	// =====================================================================
	// Diagnostic / debug / miscellaneous variables
	// =====================================================================
	real total_mass;   ///< Total solute mass in domain [mg or kg]
	realArr2 aveVB;    ///< Boundary face velocities for mass balance [m/s]. Indexed as aveVB(iGlob, 4)
	realArr2 RTcoef;   ///< Matrix coefficients for implicit solver [various]. Indexed as RTcoef(iGlob, 30)
	realArr tau;       ///< Tortuosity per cell [-]. Range: 0 to 1

	// =====================================================================
	// Debugging / testing arrays for per-component flux decomposition
	// =====================================================================
	realArr2 c_advxx, c_advyy, c_advzz;  ///< Advection flux components per species [mg/(m²·s)]
	realArr2 c_difxx, c_difyy, c_difzz;  ///< Diffusion/dispersion diagonal components per species [mg/(m²·s)]
	realArr2 c_difxy, c_difyx;           ///< Diffusion/dispersion cross-term xy per species [mg/(m²·s)]
	realArr2 c_difxz, c_difzx;           ///< Diffusion/dispersion cross-term xz per species [mg/(m²·s)]
	realArr2 c_difyz, c_difzy;           ///< Diffusion/dispersion cross-term yz per species [mg/(m²·s)]
	realArr wz;  ///< Z-direction grid coefficient [-]

	// =====================================================================
	// Allocate all state variables for groundwater transport
	// Must be called after n_mass is set
	// =====================================================================
	inline void
	allocate(GwDomain &gdom)
	{
		// Concentration: (n_mass x nCellMem x 2) for current and previous timestep
		c = realArr3("c", n_mass, gdom.nCellMem, 2);
		c_solid = realArr3("c_solid", n_mass, gdom.nCellMem, 2);
		csw = realArr2("csw", n_mass, gdom.nCellSwMem);
		ConQss = realArr2("ConQss", n_mass, gdom.nCellSw);
		RainCon_arr = realArr("RainCon_arr", n_mass);

		// Dispersion tensor: 9 components per cell per species
		dcal = realArr3("dcal", n_mass, gdom.nCellMem, 9);

		// Boundary velocities and dispersion interface velocities
		aveVB = realArr2("aveVB", gdom.nCellMem, 4);
		q_dispersion = realArr2("q_dispersion", gdom.nCellMem, 12);
		tau = realArr("tau", gdom.nCellMem);

		RTcoef = realArr2("RTcoef", gdom.nCell, 30);

		// Debug flux decomposition arrays
		c_advxx = realArr2("c_advxx", n_mass, gdom.nCellMem);
		c_advyy = realArr2("c_advyy", n_mass, gdom.nCellMem);
		c_advzz = realArr2("c_advzz", n_mass, gdom.nCellMem);
		c_difxx = realArr2("c_difxx", n_mass, gdom.nCellMem);
		c_difyy = realArr2("c_difyy", n_mass, gdom.nCellMem);
		c_difzz = realArr2("c_difzz", n_mass, gdom.nCellMem);
		c_difxy = realArr2("c_difxy", n_mass, gdom.nCellMem);
		c_difxz = realArr2("c_difxz", n_mass, gdom.nCellMem);
		c_difyz = realArr2("c_difyz", n_mass, gdom.nCellMem);
		c_difyx = realArr2("c_difyx", n_mass, gdom.nCellMem);
		c_difzx = realArr2("c_difzx", n_mass, gdom.nCellMem);
		c_difzy = realArr2("c_difzy", n_mass, gdom.nCellMem);

		Rf = realArr("Rf", gdom.nCellMem);
		wz = realArr("wz", gdom.nCellMem);

		// Upstream weighting per cell per direction
		Up_Weighting_xp = realArr("Up_Weighting_xp", gdom.nCellMem);
		Up_Weighting_xm = realArr("Up_Weighting_xm", gdom.nCellMem);
		Up_Weighting_yp = realArr("Up_Weighting_yp", gdom.nCellMem);
		Up_Weighting_ym = realArr("Up_Weighting_ym", gdom.nCellMem);
		Up_Weighting_zp = realArr("Up_Weighting_zp", gdom.nCellMem);
		Up_Weighting_zm = realArr("Up_Weighting_zm", gdom.nCellMem);

		// Allocate multi-species reaction parameter arrays
		Kd = realArr("Kd", n_mass);
		Kf = realArr("Kf", n_mass);
		Nf = realArr("Nf", n_mass);
		Alpha_D = realArr("Alpha_D", n_mass);
		Beta_D = realArr("Beta_D", n_mass);
		Beta = realArr("Beta", n_mass);
		Lambda_1 = realArr("Lambda_1", n_mass);
		Lambda_2 = realArr("Lambda_2", n_mass);

		Decay_Rate = realArr("Decay_Rate", n_mass);

		Sorption_Type = realArr("Sorption_Type", n_mass);

		// Water content arrays for reaction term computation
		wc_mid = realArr("wc_mid", gdom.nCellMem);
		wc_original_new = realArr("wc_original_new", gdom.nCellMem);
		wc_original_old = realArr("wc_original_old", gdom.nCellMem);
	}
};

#endif

#endif
