/* ==========================================================================
 * RTBCGW.h - Boundary Conditions for Groundwater Reactive (Solute) Transport
 * ==========================================================================
 * This header defines the boundary condition framework for the subsurface
 * reactive transport (solute transport) module. It supports multiple boundary
 * types (Dirichlet, Neumann, Cauchy, surface-water exchange, free drainage),
 * multi-species transport, and parallel (MPI + Kokkos) computation.
 *
 * Key classes:
 *   RTBCGW              - A single boundary condition object (one per BC region)
 *   RTSubsurfaceBoundaries - Container holding all RTBCGW objects
 * ========================================================================== */

#ifndef _RTBCGW_H_
#define _RTBCGW_H_

#if SERGHEI_SUBSURFACE_TRANSPORT
#include "const.h"
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "GwDomain.h"
#include "GwMPI.h"
#include "GwState.h"
#include "Parallel.h"
#include "RTStateGW.h"

/* ==========================================================================
 * Section 1: Boundary Type Macros
 * ==========================================================================
 * These macros define integer identifiers for each type of reactive
 * transport boundary condition supported by the groundwater module.
 *
 * Physical interpretation:
 *   Dirichlet  - Prescribed concentration at the boundary (first-kind BC).
 *                Units: [mg/L] or [kg/m^3] depending on model convention.
 *   Neumann    - Prescribed concentration gradient / diffusive flux (second-kind BC).
 *                Units: [mg/(L*m)] or equivalent.
 *   Cauchy     - Prescribed total (advective + diffusive) mass flux (third-kind BC).
 *                Units: [mg/(m^2*s)] or equivalent.
 *   _CONST     - The boundary value is a constant (time-invariant).
 *   _T         - The boundary value varies in time, read from a time series.
 *   SWE        - Surface-water / groundwater exchange boundary (dynamic coupling).
 *   FD         - Free drainage boundary (zero-gradient / outflow condition).
 * ========================================================================== */

#define SUB_RT_BC_NOFLOW 1			/* No-flow boundary: neither water nor solute crosses the boundary */
#define SUB_RT_BC_Dirichlet_CONST 2 /* Constant-concentration (Dirichlet) boundary */
#define SUB_RT_BC_Neumann_CONST 3	/* Constant-flux (Neumann) boundary */
#define SUB_RT_BC_Cauchy_CONST 4	/* Constant mixed (Cauchy) boundary */
#define SUB_RT_BC_Dirichlet_T 5		/* Time-varying Dirichlet boundary (time-series concentration) */
#define SUB_RT_BC_Neumann_T 6		/* Time-varying Neumann boundary */
#define SUB_RT_BC_Cauchy_T 7		/* Time-varying Cauchy boundary */
#define SUB_RT_BC_SWE 8				/* Surface-water exchange boundary (coupled with overland flow) */
#define SUB_RT_BC_FD 9				/* Free drainage boundary (zero-gradient outflow) */

/* ==========================================================================
 * Section 2: Boundary Direction Macros
 * ==========================================================================
 * These macros identify which face of the computational domain the boundary
 * condition is applied to.  The coordinate convention is:
 *   X  - horizontal, increasing to the right
 *   Y  - horizontal, increasing upward (or northward)
 *   Z  - vertical,   increasing upward
 *
 *   XPLUS  = right   face (maximum X)
 *   XMINUS = left    face (minimum X)
 *   YPLUS  = back    face (maximum Y)
 *   YMINUS = front   face (minimum Y)
 *   ZPLUS  = direction 5, lower/bottom face in the groundwater layer order
 *   ZMINUS = direction 6, upper/surface face in the groundwater layer order
 * ========================================================================== */

#define XPLUS 1	 /* Right  face (X maximum) */
#define XMINUS 2 /* Left   face (X minimum) */
#define YPLUS 3	 /* Back   face (Y maximum) */
#define YMINUS 4 /* Front  face (Y minimum) */
#define ZPLUS 5	 /* Direction 5: lower/bottom face in current GW convention */
#define ZMINUS 6 /* Direction 6: upper/surface face in current GW convention */

/* ==========================================================================
 * Section 3: RTBCGW Class
 * ==========================================================================
 * Represents a single reactive-transport boundary condition applied to a
 * polygonal region on one face of the subsurface domain.
 *
 * Thread safety: This class is safe to invoke inside a parallel region
 * (e.g., Kokkos parallel_for / parallel_reduce) because it does not
 * modify shared state during such invocations.
 * ========================================================================== */
class RTBCGW
{

	// this class is safe to invoke in a parallel region
public:
	/* ------------------------------------------------------------------
	 * Member Data: Cell Counts
	 * ------------------------------------------------------------------
	 * ncellsBC   - Number of boundary cells found in this subdomain
	 *              for this BC region.  Range: 0 to total domain cells.
	 *
	 * ncellsIT   - Number of internal source/sink cells (e.g., injection
	 *              or extraction wells) associated with this BC.
	 * ------------------------------------------------------------------ */
	int ncellsBC = 0; // number of boundary cells found in this subdomain
	int ncellsIT = 0; // number of internal source/sink cells

	/* ------------------------------------------------------------------
	 * Member Data: Cell Index Arrays
	 * ------------------------------------------------------------------
	 * bcells     - Kokkos integer array of global indices (with halo) of
	 *              boundary cells that belong to this BC region.
	 *              These are interior cells on the domain edge.
	 *
	 * gcells     - Kokkos integer array of global indices of the ghost
	 *              cells immediately adjacent to each boundary cell (used
	 *              for setting ghost-cell concentration values).
	 *
	 * rtgwcells  - Kokkos integer array of global indices used to address
	 *              the reactive-transport state array (maps boundary cells
	 *              to their corresponding RT ghost cells).
	 *
	 * icells     - Kokkos integer array of global indices for internal
	 *              source/sink cells (e.g., wells).
	 * ------------------------------------------------------------------ */
	intArr bcells, gcells, rtgwcells, icells;

	/* ------------------------------------------------------------------
	 * Member Data: Surface-Groundwater Exchange Type
	 * ------------------------------------------------------------------
	 * swgw_type  - For each boundary cell, classifies the surface-
	 *              groundwater exchange regime:
	 *                0 = No ponding, no surface-water / groundwater exchange
	 *                1 = Ponding with large water depth (significant exchange)
	 *                2 = Ponding with small water depth (limited exchange)
	 *              Only allocated when the BC type is SUB_RT_BC_SWE.
	 * ------------------------------------------------------------------ */
	intArr swgw_type;

	/* ------------------------------------------------------------------
	 * Member Data: Boundary Metadata
	 * ------------------------------------------------------------------
	 * location    - Integer identifier for the geometric location or
	 *               region to which this BC applies (typically matches
	 *               a polygon ID from the input file).
	 *
	 * rtbctype    - The default (global) reactive-transport boundary type
	 *               for this BC region.  One of the SUB_RT_BC_* macros.
	 *               Individual species may override this via spec_bctype.
	 *
	 * isInDomain  - Flag indicating whether this BC is active in the
	 *               current subdomain (1 = active, 0 = not found).
	 *
	 * direction   - The face of the domain to which the BC is applied.
	 *               One of the XPLUS/XMINUS/YPLUS/YMINUS/ZPLUS/ZMINUS macros.
	 * ------------------------------------------------------------------ */
	int location, rtbctype, isInDomain, direction;

	/* ------------------------------------------------------------------
	 * Member Data: Boundary Values and Time Series
	 * ------------------------------------------------------------------
	 * bcvals  - Kokkos real array of boundary values for each BC cell.
	 *           For Dirichlet BCs this holds concentration [mg/L or kg/m^3];
	 *           for Neumann BCs this holds flux values [mg/(L*m)].
	 *           Size: ncellsBC (spatially varying) or 1 (uniform).
	 *
	 * bcdata  - General-purpose real array for auxiliary boundary data.
	 *
	 * ts      - Default (global) TimeSeries object used for time-varying
	 *           BC types (_T variants).  Contains (time, value) pairs
	 *           that are linearly interpolated to the current simulation
	 *           time [s].
	 * ------------------------------------------------------------------ */
	realArr bcvals, bcdata;
	TimeSeries ts;

	/* ------------------------------------------------------------------
	 * Member Data: Mass Flux Accounting
	 * ------------------------------------------------------------------
	 * Total, inflow, and outflow mass fluxes across this BC, summed over
	 * all species and all boundary cells.  Computed during
	 * applyConcentrationBC().  Units: [mg/s] or [kg/s].
	 *
	 * QMassTot      - Net total mass flux (positive = out of domain).
	 * QMassInflow   - Total mass flux entering the domain (positive value).
	 * QMassOutflow  - Total mass flux leaving the domain (positive value).
	 * ------------------------------------------------------------------ */
	real QMassTot = 0, QMassInflow = 0, QMassOutflow = 0;

	/* ------------------------------------------------------------------
	 * Member Data: Multi-Species Boundary Condition Support
	 * ------------------------------------------------------------------
	 * The following vectors allow each chemical species to have its own
	 * boundary type, values, and time series, overriding the global
	 * defaults (rtbctype, bcvals, ts).
	 *
	 * spec_bctype      - Boundary type for each species (indexed by
	 *                    species number).  Size: number of species.
	 *                    Values: one of the SUB_RT_BC_* macros.
	 *
	 * spec_bcvals      - Boundary value array for each species.  Each
	 *                    element is a realArr of size ncellsBC (spatially
	 *                    varying) or 1 (uniform).  Units: [mg/L] for
	 *                    concentration BCs, [mg/(L*m)] for flux BCs.
	 *
	 * spec_bcfile      - File name string for spatially distributed
	 *                    boundary data for each species (e.g., a raster
	 *                    or CSV file).
	 *
	 * spec_ts          - TimeSeries object for each species, used when
	 *                    the species has a time-varying BC type (_T).
	 *
	 * spec_bcval_const - Constant boundary value for each species
	 *                    (scalar).  Used when the BC type is _CONST and
	 *                    the value is spatially uniform.  Units: [mg/L].
	 *
	 * has_bcfile       - Flag per species: 1 = a boundary data file was
	 *                    provided for this species, 0 = use constant or
	 *                    time-series value.
	 * ------------------------------------------------------------------ */
	std::vector<int> spec_bctype;		  // boundary type for each species
	std::vector<realArr> spec_bcvals;	  // spatially distributed boundary values per species [mg/L]
	std::vector<std::string> spec_bcfile; // boundary data file name per species
	std::vector<TimeSeries> spec_ts;	  // time series data per species
	std::vector<real> spec_bcval_const;	  // constant boundary value per species [mg/L]
	std::vector<int> has_bcfile;		  // flag: 1 if species has a boundary data file, 0 otherwise
	std::vector<real> spec_mass_in_rate;  // per-species boundary mass inflow [mg/s]
	std::vector<real> spec_mass_out_rate; // per-species boundary mass outflow [mg/s]

	KOKKOS_INLINE_FUNCTION
	void mergeRTCoef(RTStateGW const &rt, int idom, int exterior, int mapped) const
	{
		rt.RTcoef(idom, mapped) += rt.RTcoef(idom, exterior);
		rt.RTcoef(idom, exterior) = 0.0;
	}

	KOKKOS_INLINE_FUNCTION
	void moveRTCoefToRHS(RTStateGW const &rt, int idom, int exterior, real c_external) const
	{
		rt.RTcoef(idom, 7) -= rt.RTcoef(idom, exterior) * c_external;
		rt.RTcoef(idom, exterior) = 0.0;
	}

	KOKKOS_INLINE_FUNCTION
	void clearRTCoef(RTStateGW const &rt, int idom, int exterior) const
	{
		rt.RTcoef(idom, exterior) = 0.0;
	}

	KOKKOS_INLINE_FUNCTION
	void applyConservativeFDElimination(RTStateGW const &rt, int idom, int bc_direction) const
	{
		if (bc_direction == XPLUS)
		{
			mergeRTCoef(rt, idom, 1, 0);
			mergeRTCoef(rt, idom, 9, 3);
			mergeRTCoef(rt, idom, 11, 4);
			mergeRTCoef(rt, idom, 15, 5);
			mergeRTCoef(rt, idom, 17, 6);
		}
		else if (bc_direction == XMINUS)
		{
			mergeRTCoef(rt, idom, 2, 0);
			mergeRTCoef(rt, idom, 12, 3);
			mergeRTCoef(rt, idom, 13, 4);
			mergeRTCoef(rt, idom, 18, 5);
			mergeRTCoef(rt, idom, 19, 6);
		}
		else if (bc_direction == YPLUS)
		{
			mergeRTCoef(rt, idom, 3, 0);
			mergeRTCoef(rt, idom, 9, 1);
			mergeRTCoef(rt, idom, 12, 2);
			mergeRTCoef(rt, idom, 22, 5);
			mergeRTCoef(rt, idom, 23, 6);
		}
		else if (bc_direction == YMINUS)
		{
			mergeRTCoef(rt, idom, 4, 0);
			mergeRTCoef(rt, idom, 11, 1);
			mergeRTCoef(rt, idom, 13, 2);
			mergeRTCoef(rt, idom, 24, 5);
			mergeRTCoef(rt, idom, 25, 6);
		}
		else if (bc_direction == ZPLUS)
		{
			mergeRTCoef(rt, idom, 5, 0);
			mergeRTCoef(rt, idom, 15, 1);
			mergeRTCoef(rt, idom, 18, 2);
			mergeRTCoef(rt, idom, 22, 3);
			mergeRTCoef(rt, idom, 24, 4);
		}
		else if (bc_direction == ZMINUS)
		{
			mergeRTCoef(rt, idom, 6, 0);
			mergeRTCoef(rt, idom, 17, 1);
			mergeRTCoef(rt, idom, 19, 2);
			mergeRTCoef(rt, idom, 23, 3);
			mergeRTCoef(rt, idom, 25, 4);
		}
	}

	KOKKOS_INLINE_FUNCTION
	void eliminateFDInflowGhost(RTStateGW const &rt, int idom, int bc_direction, real c_external) const
	{
		if (bc_direction == XPLUS)
		{
			moveRTCoefToRHS(rt, idom, 1, c_external);
			clearRTCoef(rt, idom, 9);
			clearRTCoef(rt, idom, 11);
			clearRTCoef(rt, idom, 15);
			clearRTCoef(rt, idom, 17);
		}
		else if (bc_direction == XMINUS)
		{
			moveRTCoefToRHS(rt, idom, 2, c_external);
			clearRTCoef(rt, idom, 12);
			clearRTCoef(rt, idom, 13);
			clearRTCoef(rt, idom, 18);
			clearRTCoef(rt, idom, 19);
		}
		else if (bc_direction == YPLUS)
		{
			moveRTCoefToRHS(rt, idom, 3, c_external);
			clearRTCoef(rt, idom, 9);
			clearRTCoef(rt, idom, 12);
			clearRTCoef(rt, idom, 22);
			clearRTCoef(rt, idom, 23);
		}
		else if (bc_direction == YMINUS)
		{
			moveRTCoefToRHS(rt, idom, 4, c_external);
			clearRTCoef(rt, idom, 11);
			clearRTCoef(rt, idom, 13);
			clearRTCoef(rt, idom, 24);
			clearRTCoef(rt, idom, 25);
		}
		else if (bc_direction == ZPLUS)
		{
			moveRTCoefToRHS(rt, idom, 5, c_external);
			clearRTCoef(rt, idom, 15);
			clearRTCoef(rt, idom, 18);
			clearRTCoef(rt, idom, 22);
			clearRTCoef(rt, idom, 24);
		}
		else if (bc_direction == ZMINUS)
		{
			moveRTCoefToRHS(rt, idom, 6, c_external);
			clearRTCoef(rt, idom, 17);
			clearRTCoef(rt, idom, 19);
			clearRTCoef(rt, idom, 23);
			clearRTCoef(rt, idom, 25);
		}
	}

	/* ------------------------------------------------------------------
	 * Member Data: MPI Communicator
	 * ------------------------------------------------------------------
	 * comm - MPI communicator that includes only the ranks (subdomains)
	 *        that contain at least one boundary cell for this BC region.
	 *        Used for collective operations specific to this BC (e.g.,
	 *        summing mass fluxes across all participating ranks).
	 * ------------------------------------------------------------------ */
	MPI_Comm comm;

	/* ==================================================================
	 * Method: find_bcells
	 * ==================================================================
	 * Locates all cells within the current MPI subdomain that lie inside
	 * the specified polygonal boundary region, and identifies the
	 * corresponding ghost cells in the halo layer.  Sets up the MPI
	 * communicator for this BC region.
	 *
	 * Parameters:
	 *   gw     - (in) Groundwater state (used for context; not modified).
	 *   rtid   - (in) String identifier for this boundary condition
	 *            (used in error messages).
	 *   gdom   - (in) Groundwater domain geometry (grid spacing,
	 *            dimensions, halo width, origin coordinates).
	 *   par    - (in) Parallel decomposition info (MPI rank, subdomain
	 *            ranges, number of ranks).
	 *   nPoly  - (in) Number of vertices in the boundary polygon.
	 *   xPoly  - (in) X-coordinates of polygon vertices [m].
	 *   yPoly  - (in) Y-coordinates of polygon vertices [m].
	 *
	 * Returns:
	 *   1 on success (boundary cells found across all subdomains),
	 *   0 on failure (no boundary cells found anywhere).
	 *
	 * Notes:
	 *   - For top/bottom boundaries (ZPLUS/ZMINUS), only the outermost
	 *     layer of cells is selected.
	 *   - For lateral boundaries, all vertical layers are included.
	 *   - Allocates and populates bcells, gcells, and rtgwcells arrays.
	 *   - Creates an MPI subgroup communicator (comm) for participating
	 *     ranks.
	 * ================================================================== */
	inline int find_bcells(GwState &gw, std::string &rtid, GwDomain &gdom, Parallel &par, int nPoly, realArr &xPoly, realArr &yPoly)
	{
		int foundInSubdom;			 // rank ID of the current subdomain if BC cells are found here; -1 otherwise
		std::vector<int> tmpbcells;	 // temporary list of boundary-cell global indices (with halo)
		std::vector<int> tmpgcells;	 // temporary list of ghost-cell global indices (with halo)
		std::vector<int> subdomains; // list of MPI ranks that contain at least one BC cell

		// Loop over all cells in this subdomain (including halo) to find BC cells
		for (int kk = 0; kk < gdom.nz; kk++)
		{
			for (int jj = 0; jj < gdom.ny; jj++)
			{
				for (int ii = 0; ii < gdom.nx; ii++)
				{
					// Compute the global index (with halo) of the current cell.
					// gdom.hc is the halo width in cells.
					int iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;
					foundInSubdom = -1;

					// Compute the physical (x, y) coordinates of the cell center.
					// par.i_beg and par.j_beg are the global starting indices of this
					// subdomain in the x and y directions, respectively.
					// Units: [m].  gdom.dx is the grid spacing [m].
					real xCoord = gdom.xll + (par.i_beg + ii + 0.5) * gdom.dx;
					real yCoord = gdom.yll + gdom.ny_glob * gdom.dx - (par.j_beg + jj + 0.5) * gdom.dx;

					// Check if the cell center lies inside the boundary polygon
					if (geometry::isInsidePoly(nPoly, xPoly, yPoly, xCoord, yCoord))
					{
						// Keep RT boundary-cell selection consistent with GwBC:
						// direction 5 selects kk == nz-1 (bottom/lower face),
						// direction 6 selects kk == 0 (surface/upper face).
						if (direction == 5)
						{
							if (kk == gdom.nz - 1)
							{
								tmpbcells.push_back(iGlob);
								tmpgcells.push_back(iGlob + gdom.nxhc * gdom.nyhc);
							}
						}
						else if (direction == 6)
						{
							if (kk == 0)
							{
								tmpbcells.push_back(iGlob);
								tmpgcells.push_back(iGlob - gdom.nxhc * gdom.nyhc);
							}
						}
						// For lateral boundaries, all vertical layers are included.
						// NOTE: In the future, this should be customized to allow only
						// certain vertical layers to be included.
						else
						{
							tmpbcells.push_back(iGlob);
							if (direction == XPLUS)
							{
								// Right boundary: ghost cell is one column to the right (iGlob + 1).
								tmpgcells.push_back(iGlob + 1);
							}
							else if (direction == XMINUS)
							{
								// Left boundary: ghost cell is one column to the left (iGlob - 1).
								tmpgcells.push_back(iGlob - 1);
							}
							else if (direction == YPLUS)
							{
								// Back boundary: ghost cell is one row above (offset by nxhc).
								tmpgcells.push_back(iGlob + gdom.nxhc);
							}
							else if (direction == YMINUS)
							{
								// Front boundary: ghost cell is one row below (offset by -nxhc).
								tmpgcells.push_back(iGlob - gdom.nxhc);
							}
						}
					}
				}
			}
		}

		ncellsBC = int(tmpbcells.size());
		if (ncellsBC > 0)
			foundInSubdom = par.myrank; // Tag this subdomain as containing BC cells

		// Sum the number of BC cells across all MPI ranks to get the global total.
		int ncells_all;
		SERGHEI_MPI_CHECK(MPI_Allreduce(&ncellsBC, &ncells_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));

		// Gather which ranks have BC cells.
		// subdoms[i] >= 0 means rank i participates in this BC.
		int *subdoms;
		subdoms = (int *)malloc(par.nranks * sizeof(int));
		SERGHEI_MPI_CHECK(MPI_Allgather(&foundInSubdom, 1, MPI_INT, subdoms, 1, MPI_INT, MPI_COMM_WORLD));
#if SERGHEI_DEBUG_BOUNDARY
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "GW ncellsBC " << ncells_all << std::endl;
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "GW BC subdomains: ";
		for (int i = 0; i < par.nranks; i++)
		{
			std::cout << GGD << " ";
			if (subdoms[i] == par.myrank)
				std::cout << RED;
			std::cout << subdoms[i] << "\t" << RESET;
		}
		std::cout << std::endl;
#endif

		// Build a consolidated list of participating ranks (those with subdoms[i] >= 0).
		for (int i = 0; i < par.nranks; i++)
		{
			if (subdoms[i] >= 0)
			{
				subdomains.push_back(subdoms[i]);
			}
		}
#if SERGHEI_DEBUG_BOUNDARY
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "GW Consolidated BC subdomains = ";
		for (int i = 0; i < subdomains.size(); i++)
		{
			std::cout << GGD << " ";
			std::cout << subdomains[i] << "\t";
		}
		std::cout << std::endl;
#endif
		// Create an MPI communicator that includes only the participating ranks,
		// so that collective operations on this BC only involve relevant processes.
		MPI_Group group, subgroup;
		SERGHEI_MPI_CHECK(MPI_Comm_group(MPI_COMM_WORLD, &group));
		SERGHEI_MPI_CHECK(MPI_Group_incl(group, subdomains.size(), subdomains.data(), &subgroup));
		SERGHEI_MPI_CHECK(MPI_Comm_create(MPI_COMM_WORLD, subgroup, &comm));

		// If boundary cells were found, copy the temporary lists into Kokkos
		// device arrays (or host arrays for non-CUDA builds).  If no cells
		// were found across the entire domain, report an error and return failure.
		if (ncells_all > 0)
		{
			bcells = intArr("bcells", ncellsBC);
			gcells = intArr("gcells", ncellsBC);
			rtgwcells = intArr("rtgwcells", ncellsBC);

			if (rtbctype == SUB_RT_BC_SWE)
			{
				swgw_type = intArr("swgw_type", ncellsBC);
			}
#ifdef __NVCC__
			// Asynchronously copy cell-index data from host to GPU device
			cudaMemcpyAsync(bcells.data(), tmpbcells.data(), ncellsBC * sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpyAsync(gcells.data(), tmpgcells.data(), ncellsBC * sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpyAsync(rtgwcells.data(), tmpgcells.data(), ncellsBC * sizeof(int), cudaMemcpyHostToDevice);
			cudaDeviceSynchronize();
#else
			// On CPU, simply copy the host data into the Kokkos view
			std::memcpy(bcells.data(), tmpbcells.data(), ncellsBC * sizeof(int));
			std::memcpy(gcells.data(), tmpgcells.data(), ncellsBC * sizeof(int));
			std::memcpy(rtgwcells.data(), tmpgcells.data(), ncellsBC * sizeof(int));
#endif
		}
		else
		{
			if (par.masterproc)
			{
				std::cerr << RERROR << "No boundary cells found for subsurface reaction_transport boundary with id '" << rtid << "'" << std::endl;
			}
			return 0;
		}
		return 1;
	}

	/* ==================================================================
	 * Method: applyConcentrationBC
	 * ==================================================================
	 * Applies concentration-type boundary conditions for a specific
	 * chemical species by setting the ghost-cell concentration values.
	 * Legacy advective mass-flux accounting is only updated when
	 * calculate_fluxes is true.  The groundwater transport driver normally
	 * leaves this false and fills the mass rates after the linear solve via
	 * boundary-type-specific diagnostics.
	 *
	 * This method is called once per species per time step.
	 *
	 * Parameters:
	 *   rt          - (in/out) Reactive transport state for groundwater.
	 *                On output, ghost-cell concentrations (rt.c) are set
	 *                for the target species at all boundary cells.
	 *   gw          - (in)     Groundwater flow state (head, flux, etc.).
	 *   gdom        - (in)     Groundwater domain geometry.
	 *   par         - (in)     MPI parallel decomposition information.
	 *   target_spec - (in)     Zero-based index of the chemical species
	 *                to apply the BC for.  Default: 0 (first species).
	 *   calculate_fluxes - (in) If true, compute the legacy v*C*A flux
	 *                diagnostic.  Default: false.
	 *
	 * Returns:
	 *   void.  Results are stored in rt.c (concentrations) and in the
	 *   member variables QMassTot, QMassInflow, QMassOutflow [mg/s] when
	 *   calculate_fluxes is true.
	 *
	 * BC type behavior:
	 *   - Cauchy / FD: zero-gradient condition (ghost = interior).
	 *   - Dirichlet_CONST: ghost cell set to a constant or spatially
	 *     varying concentration from bcvals or spec_bcvals.
	 *   - Dirichlet_T: ghost cell set to a time-interpolated value.
	 *   - SWE: ghost cell concentration determined by the direction of
	 *     surface-water / groundwater exchange flux (infiltration vs.
	 *     exfiltration).
	 * ================================================================== */
	inline void applyConcentrationBC(RTStateGW &rt, GwState &gw, GwDomain &gdom, Parallel &par, int target_spec = 0, bool calculate_fluxes = false)
	{
		// ConQss is shared by all boundary objects and species.  It must be
		// cleared by the transport driver before a complete boundary refresh,
		// not here: this routine is called once for every boundary, so clearing
		// on target_spec == 0 would erase the species-0 SWE flux when the next
		// (for example, bottom free-drainage) boundary is processed.
		if (spec_mass_in_rate.size() != rt.n_mass)
		{
			spec_mass_in_rate.assign(rt.n_mass, 0.0);
			spec_mass_out_rate.assign(rt.n_mass, 0.0);
		}
		if (calculate_fluxes)
		{
			// Legacy advective-flux diagnostics.  The main GW transport
			// driver now uses boundary-type-specific diagnostics after the
			// linear solve.
			QMassTot = 0;
			QMassInflow = 0;
			QMassOutflow = 0;
			if (target_spec >= 0 && target_spec < rt.n_mass)
			{
				spec_mass_in_rate[target_spec] = 0.0;
				spec_mass_out_rate[target_spec] = 0.0;
			}
		}

		// Determine whether this subdomain actually lies on the global
		// domain boundary face specified by 'direction'.  Only subdomains
		// on the domain edge should apply the BC.
		bool onBoundary = 0;
		if (direction == 2 && par.px == 0) // Left face: first column of processes
		{
			onBoundary = 1;
		}
		else if (direction == 1 && par.px == par.nproc_x - 1) // Right face: last column of processes
		{
			onBoundary = 1;
		}
		else if (direction == 4 && par.py == 0) // Front face: first row of processes
		{
			onBoundary = 1;
		}
		else if (direction == 3 && par.py == par.nproc_y - 1) // Back face: last row of processes
		{
			onBoundary = 1;
		}
		else if (direction == 5 || direction == 6) // Top/bottom: always a boundary (Z faces)
		{
			onBoundary = 1;
		}

		if (ncellsBC > 0 && onBoundary == 1)
		{
			// Determine the effective boundary type for this species.
			// If multi-species BC data exists, use the per-species type;
			// otherwise fall back to the global default (rtbctype).
			int current_bctype = rtbctype;
			realArr current_bcvals = bcvals;
			TimeSeries current_ts = ts;
			if (target_spec < spec_bctype.size())
			{
				current_bctype = spec_bctype[target_spec];
			}
			if (target_spec < spec_bcvals.size() && spec_bcvals[target_spec].size() > 0)
			{
				current_bcvals = spec_bcvals[target_spec];
			}
			if (target_spec < spec_ts.size() && spec_ts[target_spec].np > 0)
			{
				current_ts = spec_ts[target_spec];
			}

			real Conbc = 0.0;
			// For time-varying Dirichlet BC, interpolate the boundary
			// concentration from the species-specific time series.
			if (current_bctype == SUB_RT_BC_Dirichlet_T || current_bctype == SUB_RT_BC_Cauchy_T)
			{
				// Use the per-species time series if available
				if (current_ts.np > 0)
				{
					Conbc = interpolateLinear(current_ts, gdom.etime);
				}
				else
				{
					if (ncellsBC > 0)
					{
						std::cerr << RERROR "GWRTM: Boundary time-series BC: No timeseries data for species " << target_spec << std::endl;
					}
					Conbc = 0.0;
				}
			}

			/* ----------------------------------------------------------
			 * Cauchy BCs: zero-gradient (Neumann) condition.
			 * ---------------------------------------------------------- */
			if (current_bctype == SUB_RT_BC_Cauchy_CONST || current_bctype == SUB_RT_BC_Cauchy_T)
			{
				Kokkos::parallel_for("rt_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc) {
					int iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc];
					rt.c(target_spec, rtiGhost, 1) = rt.c(target_spec, iGlob, 1); });
			}

			/* ----------------------------------------------------------
			 * Free-drainage BCs are outflow-only.  Use zero gradient only
			 * where the local Darcy velocity leaves the modeled domain.  If
			 * the hydraulic solution produces local backflow, prescribe the
			 * external/background concentration instead of recycling the
			 * interior concentration through the ghost cell.
			 * ---------------------------------------------------------- */
			else if (current_bctype == SUB_RT_BC_FD)
			{
				real fd_inflow_conc = 0.0;
				if (current_bcvals.size() > 0)
					fd_inflow_conc = current_bcvals(0);

				Kokkos::parallel_for("rt_bc_fd_ghost", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc) {
					int iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc];
					const real FLOW_EPSILON = 1.0e-12;
					real vel = 0.0;
					bool inflow = false;

					if (direction == 1)
					{
						vel = rt.aveVB(iGlob, 0);
						inflow = vel < -FLOW_EPSILON;
					}
					else if (direction == 2)
					{
						vel = rt.aveVB(iGlob - 1, 0);
						inflow = vel > FLOW_EPSILON;
					}
					else if (direction == 3)
					{
						vel = rt.aveVB(iGlob, 1);
						inflow = vel < -FLOW_EPSILON;
					}
					else if (direction == 4)
					{
						vel = rt.aveVB(iGlob - gdom.nxhc, 1);
						inflow = vel > FLOW_EPSILON;
					}
					else if (direction == 5)
					{
						vel = rt.aveVB(iGlob, 2);
						inflow = vel < -FLOW_EPSILON;
					}
					else if (direction == 6)
					{
						vel = rt.aveVB(iGlob - gdom.nxhc * gdom.nyhc, 2);
						inflow = vel > FLOW_EPSILON;
					}

					rt.c(target_spec, rtiGhost, 1) = inflow ? fd_inflow_conc : rt.c(target_spec, iGlob, 1);
				});
			}

			/* ----------------------------------------------------------
			 * Dirichlet constant BC: ghost-cell concentration is set to a
			 * prescribed value.  If multi-species spatial data exists, use
			 * the per-species array; otherwise use the global bcvals array.
			 * ---------------------------------------------------------- */
			else if (current_bctype == SUB_RT_BC_Dirichlet_CONST)
			{
				// Use per-species boundary values if available
				if (current_bcvals.size() > 1)
				{
					Kokkos::parallel_for("rt_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc) {
					                 int iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc];
					                 rt.c(target_spec, rtiGhost, 1) = current_bcvals(ibc); });
				}
				else
				{
					Kokkos::parallel_for("rt_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc) {
					                 int iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc];
					                 rt.c(target_spec, rtiGhost, 1) = current_bcvals(0); });
				}
			}

			/* ----------------------------------------------------------
			 * Dirichlet time-series BC: ghost-cell concentration is set to
			 * the time-interpolated value Conbc computed above.
			 * ---------------------------------------------------------- */
			else if (current_bctype == SUB_RT_BC_Dirichlet_T)
			{
				Kokkos::parallel_for("rt_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc) {
					                 int iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc];
					                 rt.c(target_spec, rtiGhost, 1) = Conbc; });
			}

			/* ----------------------------------------------------------
			 * Surface-water / groundwater exchange BC (SUB_RT_BC_SWE).
			 *
			 * The SWE interface is a pure advective exchange boundary.
			 * Groundwater dispersion remains active inside the domain, but
			 * no dispersive flux is allowed across the surface interface.
			 * The matrix receives the explicit interface source/sink in
			 * applyRTMatBC(); the ghost value is therefore kept equal to the
			 * interior value and is never used to inject surface dispersion.
			 *
			 * ConQss keeps the historical groundwater sign convention:
			 *   ConQss = -J_SW_to_GW,
			 * where J_SW_to_GW is positive from surface water to groundwater.
			 * This makes ConQss positive for groundwater -> surface inflow.
			 * ---------------------------------------------------------- */
			else if (current_bctype == SUB_RT_BC_SWE)
			{

				if (direction == 6)
				{
					Kokkos::parallel_for("rt_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc) {
							int iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc];
							int ii, jj, kk, iGlobSW;
							gdom.unpackIndicesHalo(iGlob, kk, jj, ii);

						// Map from subsurface halo index to surface grid index
						iGlobSW = jj * gdom.nxhc + ii;				// with halo
						int iGlobSW2 = (jj - 1) * gdom.nx + ii - 1; // without halo

						// Small flow-velocity threshold to avoid numerical noise
						// causing spurious boundary-direction reversal. [m/s]
						const real FLOW_EPSILON = 1.0e-10;

							// The ghost is not a concentration source for this boundary.
							// Keeping it equal to the interior also enforces zero normal
							// dispersive flux if any external stencil entry remains.
							rt.c(target_spec, rtiGhost, 1) = rt.c(target_spec, iGlob, 1);

							real qss = gw.qss(iGlobSW2);
							if (qss != qss || fabs(qss) <= FLOW_EPSILON)
							{
								rt.ConQss(target_spec, iGlobSW2) = 0.0;
							}
							else
							{
								// During infiltration, surface water is the upstream donor.
								real Csw = rt.csw(target_spec, iGlobSW);
								if (gdom.isRain == 1 && gdom.rainRate(iGlobSW) > 0.0 && Csw == 0.0)
									Csw = rt.RainCon_arr(target_spec);

								real q_in = fmax(-qss, 0.0);
								real q_out = fmax(qss, 0.0);
								real J_sw_to_gw = q_in * Csw - q_out * rt.c(target_spec, iGlob, 1);

								// ConQss is positive for mass entering surface water.
								rt.ConQss(target_spec, iGlobSW2) = -J_sw_to_gw;
							}
						});
				}
				else
				{
					if (par.masterproc)
					{
						std::cerr << RERROR "RTBC direction must be 6 for RT-SW-GW exchange boundary! " << "\n";
					}
				}
			}

			/* ----------------------------------------------------------
			 * Mass balance calculation: compute the solute mass flux for
			 * the current species across the boundary.  The flux is:
			 *   Q_mass = v * C * A
			 * where v is the Darcy velocity (aveVB) [m/s],
			 *       C is the concentration [mg/L],
			 *       A is the cross-sectional area [m^2].
			 * The conversion factor 1000 changes mg/L to mg/m^3.
			 * The total mass flux QMassTot has units [mg/s].
			 *
			 * The sign convention distinguishes inflow from outflow
			 * based on the boundary face direction:
			 *   - Right (X+): positive velocity = outflow
			 *   - Left  (X-): positive velocity = inflow
			 *   - Back  (Y+): positive velocity = outflow
			 *   - Front (Y-): positive velocity = inflow
			 *   - Direction 5 (bottom): positive aveVB_z = outflow
			 *   - Direction 6 (top):    positive aveVB_z = inflow
			 *
			 * Note: rt.aveVB is transformed from the groundwater flow
			 * convention as aveVB = -gw.q_new.  Therefore the vertical
			 * sign test below must follow aveVB, not the raw GwBC q sign.
			 * ---------------------------------------------------------- */

			if (calculate_fluxes)
			{
				if (direction == 1) // Right boundary (X maximum)
				{
					// Use the X-component of Darcy velocity (aveVB index 0).
					// Cross-section area for an X-face = dy * dz.  [m^2]
					int bcval_size = current_bcvals.size();
					Kokkos::parallel_reduce("reducexRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
								int iGlob = bcells[ibc], rtiGhost = rtgwcells[ibc];
								real vel = rt.aveVB(iGlob, 0);
								real c_in = rt.c(target_spec, rtiGhost, 1);
								if (current_bctype == SUB_RT_BC_Cauchy_CONST || current_bctype == SUB_RT_BC_Dirichlet_CONST)
									c_in = (bcval_size > 1) ? current_bcvals(ibc) : current_bcvals(0);
								else if (current_bctype == SUB_RT_BC_Cauchy_T || current_bctype == SUB_RT_BC_Dirichlet_T)
									c_in = Conbc;
								real c_flux = (vel < 0.0) ? c_in : rt.c(target_spec, iGlob, 1);
								tmp += vel * c_flux * gdom.dy * gdom.dz(iGlob) * 1000.0; }, Kokkos::Sum<real>(QMassTot));
					if (QMassTot > 0)
					{
						QMassOutflow = QMassTot;
					}
					else
					{
						QMassInflow = -QMassTot;
					}
				}
				else if (direction == 2) // Left boundary (X minimum)
				{
					int bcval_size = current_bcvals.size();
					Kokkos::parallel_reduce("reducexRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
								int iGlob = bcells[ibc], rtiGhost = rtgwcells[ibc];
								real vel = rt.aveVB(iGlob - 1, 0);
								real c_in = rt.c(target_spec, rtiGhost, 1);
								if (current_bctype == SUB_RT_BC_Cauchy_CONST || current_bctype == SUB_RT_BC_Dirichlet_CONST)
									c_in = (bcval_size > 1) ? current_bcvals(ibc) : current_bcvals(0);
								else if (current_bctype == SUB_RT_BC_Cauchy_T || current_bctype == SUB_RT_BC_Dirichlet_T)
									c_in = Conbc;
								real c_flux = (vel > 0.0) ? c_in : rt.c(target_spec, iGlob, 1);
								tmp += vel * c_flux * gdom.dy * gdom.dz(iGlob) * 1000.0; }, Kokkos::Sum<real>(QMassTot));
					if (QMassTot > 0)
					{
						QMassInflow = QMassTot;
					}
					else
					{
						QMassOutflow = -QMassTot;
					}
				}

				// --- Y direction boundaries ---

				if (direction == 3) // Back boundary (Y maximum)
				{
					// Use the Y-component of Darcy velocity (aveVB index 1).
					// Cross-section area for a Y-face = dx * dz.  [m^2]
					int bcval_size = current_bcvals.size();
					Kokkos::parallel_reduce("reduceyRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
								int iGlob = bcells[ibc], rtiGhost = rtgwcells[ibc];
								real vel = rt.aveVB(iGlob, 1);
								real c_in = rt.c(target_spec, rtiGhost, 1);
								if (current_bctype == SUB_RT_BC_Cauchy_CONST || current_bctype == SUB_RT_BC_Dirichlet_CONST)
									c_in = (bcval_size > 1) ? current_bcvals(ibc) : current_bcvals(0);
								else if (current_bctype == SUB_RT_BC_Cauchy_T || current_bctype == SUB_RT_BC_Dirichlet_T)
									c_in = Conbc;
								real c_flux = (vel < 0.0) ? c_in : rt.c(target_spec, iGlob, 1);
								tmp += vel * c_flux * gdom.dx * gdom.dz(iGlob) * 1000.0; }, Kokkos::Sum<real>(QMassTot));
					if (QMassTot > 0)
					{
						QMassOutflow = QMassTot;
					}
					else
					{
						QMassInflow = -QMassTot;
					}
				}
				else if (direction == 4) // Front boundary (Y minimum)
				{
					int bcval_size = current_bcvals.size();
					Kokkos::parallel_reduce("reduceyRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
								int iGlob = bcells[ibc], rtiGhost = rtgwcells[ibc];
								real vel = rt.aveVB(iGlob - gdom.nxhc, 1);
								real c_in = rt.c(target_spec, rtiGhost, 1);
								if (current_bctype == SUB_RT_BC_Cauchy_CONST || current_bctype == SUB_RT_BC_Dirichlet_CONST)
									c_in = (bcval_size > 1) ? current_bcvals(ibc) : current_bcvals(0);
								else if (current_bctype == SUB_RT_BC_Cauchy_T || current_bctype == SUB_RT_BC_Dirichlet_T)
									c_in = Conbc;
								real c_flux = (vel > 0.0) ? c_in : rt.c(target_spec, iGlob, 1);
								tmp += vel * c_flux * gdom.dx * gdom.dz(iGlob) * 1000.0; }, Kokkos::Sum<real>(QMassTot));
					if (QMassTot > 0)
					{
						QMassInflow = QMassTot;
					}
					else
					{
						QMassOutflow = -QMassTot;
					}
				}

				else if (direction == 5)
				{
					// Use the Z-component of Darcy velocity (aveVB index 2).
					// Cross-section area for a Z-face = dx * dy.  [m^2]
					int bcval_size = current_bcvals.size();
					Kokkos::parallel_reduce("reducezRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
								int iGlob = bcells[ibc], rtiGhost = rtgwcells[ibc];
								real vel = rt.aveVB(iGlob, 2);
								real c_in = rt.c(target_spec, rtiGhost, 1);
								if (current_bctype == SUB_RT_BC_Cauchy_CONST || current_bctype == SUB_RT_BC_Dirichlet_CONST)
									c_in = (bcval_size > 1) ? current_bcvals(ibc) : current_bcvals(0);
								else if (current_bctype == SUB_RT_BC_Cauchy_T || current_bctype == SUB_RT_BC_Dirichlet_T)
									c_in = Conbc;
								real c_flux = (vel < 0.0) ? c_in : rt.c(target_spec, iGlob, 1);
								tmp += vel * c_flux * gdom.dx * gdom.dy * 1000.0; }, Kokkos::Sum<real>(QMassTot));
					if (QMassTot > 0)
					{
						QMassOutflow = QMassTot;
					}
					else
					{
						QMassInflow = -QMassTot;
					}
				}
				else if (direction == 6)
				{
					int bcval_size = current_bcvals.size();
					Kokkos::parallel_reduce("reducezRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
								int iGlob = bcells[ibc], rtiGhost = rtgwcells[ibc];
								real vel = rt.aveVB(iGlob - gdom.nxhc * gdom.nyhc, 2);
								real c_in = rt.c(target_spec, rtiGhost, 1);
								if (current_bctype == SUB_RT_BC_Cauchy_CONST || current_bctype == SUB_RT_BC_Dirichlet_CONST)
									c_in = (bcval_size > 1) ? current_bcvals(ibc) : current_bcvals(0);
								else if (current_bctype == SUB_RT_BC_Cauchy_T || current_bctype == SUB_RT_BC_Dirichlet_T)
									c_in = Conbc;
								real c_flux = (vel > 0.0) ? c_in : rt.c(target_spec, iGlob, 1);
								tmp += vel * c_flux * gdom.dx * gdom.dy * 1000.0; }, Kokkos::Sum<real>(QMassTot));
					if (QMassTot > 0)
					{
						QMassInflow = QMassTot;
					}
					else
					{
						QMassOutflow = -QMassTot;
					}
				}
				spec_mass_in_rate[target_spec] = QMassInflow;
				spec_mass_out_rate[target_spec] = QMassOutflow;
			}
		}
	}

	/* ==================================================================
	 * Method: applyRTMatBC
	 * ==================================================================
	 * Applies boundary conditions to the matrix coefficients (RTcoef)
	 * of the linear system that is solved for the implicit reactive
	 * transport step.  This modifies the coefficients of the discrete
	 * transport equation at boundary cells to enforce the prescribed
	 * boundary condition type.
	 *
	 * The coefficient array RTcoef stores the 7-point stencil plus RHS:
	 *   Index 0: Center cell coefficient (diagonal, LHS)
	 *   Index 1: +X neighbor coefficient (i+1)
	 *   Index 2: -X neighbor coefficient (i-1)
	 *   Index 3: +Y neighbor coefficient (j+1)
	 *   Index 4: -Y neighbor coefficient (j-1)
	 *   Index 5: +Z neighbor coefficient (k+1)
	 *   Index 6: -Z neighbor coefficient (k-1)
	 *   Index 7: Right-hand side (RHS) source/sink term
	 *
	 * Parameters:
	 *   rt          - (in/out) Reactive transport state.  On output,
	 *                RTcoef for boundary cells is modified.
	 *   gw          - (in)     Groundwater flow state (head, Darcy flux).
	 *   gdom        - (in)     Groundwater domain geometry.
	 *   par         - (in)     MPI parallel decomposition information.
	 *   target_spec - (in)     Zero-based index of the chemical species.
	 *                Default: 0 (first species).
	 *
	 * Returns:
	 *   void.  Results are stored in rt.RTcoef.
	 *
	 * BC type behavior:
	 *   - Dirichlet_CONST / Dirichlet_T:
	 *     Enforces a prescribed concentration by exact row replacement:
	 *     all physical stencil entries are cleared, the diagonal is set
	 *     to one, and the RHS is set to C_boundary.
	 *     For the left boundary (direction 2), a variable water-table
	 *     adjustment is applied: cells above the water table use a
	 *     zero-gradient outflow, while saturated cells with inflow use
	 *     the Dirichlet penalty.
	 *
	 *   - SWE (surface-water exchange):
	 *     Adds source/sink terms based on the exchange flux:
	 *       Infiltration (qss < 0): mass source added to RHS.
	 *       Exfiltration (qss > 0): mass sink added to LHS diagonal.
	 *
	 *   - Cauchy_CONST / Cauchy_T:
	 *     Adds flux-proportional source (inflow) or sink (outflow)
	 *     terms to the matrix, handling all six directions (X, Y, Z)
	 *     with correct sign conventions.
	 *
	 *   - FD (free drainage):
	 *     Conservative zero-gradient condition: keeps the assembled
	 *     physical equation, merges exterior ghost-cell couplings into
	 *     the corresponding interior stencil entries, and treats local
	 *     backflow with an external/background concentration.
	 * ================================================================== */
	inline void applyRTMatBC(RTStateGW &rt, GwState &gw, GwDomain &gdom, Parallel &par, int target_spec = 0)
	{

		// Determine whether this subdomain is on the global boundary face.
		// Same logic as in applyConcentrationBC().
		bool onBoundary = 0;
		if (direction == 2 && par.px == 0)
		{
			onBoundary = 1;
		}
		else if (direction == 1 && par.px == par.nproc_x - 1)
		{
			onBoundary = 1;
		}
		else if (direction == 4 && par.py == 0)
		{
			onBoundary = 1;
		}
		else if (direction == 3 && par.py == par.nproc_y - 1)
		{
			onBoundary = 1;
		}
		else if (direction == 5 || direction == 6)
		{
			onBoundary = 1;
		}

		if (ncellsBC > 0 && onBoundary == 1)
		{
			// Resolve the effective boundary type, values, and time series
			// for the target species (multi-species support).
			int current_bctype = rtbctype;
			realArr current_bcvals = bcvals;
			TimeSeries current_ts = ts;

			// Override with per-species data if available
			if (target_spec < spec_bctype.size())
			{
				current_bctype = spec_bctype[target_spec];
			}
			if (target_spec < spec_bcvals.size() && spec_bcvals[target_spec].size() > 0)
			{
				current_bcvals = spec_bcvals[target_spec];
			}
			// Override with per-species time series if available
			if (target_spec < spec_ts.size() && spec_ts[target_spec].np > 0)
			{
				current_ts = spec_ts[target_spec];
			}
			else
			{
				// Only warn if the BC type actually requires time-series data.
				// SWE boundaries do not need a time series.
				if (current_bctype == SUB_RT_BC_Cauchy_T || current_bctype == SUB_RT_BC_Dirichlet_T)
				{
					if (ncellsBC > 0)
					{
						std::cerr << RERROR "GWRTM: Boundary type " << current_bctype << ": No timeseries data for species " << target_spec << std::endl;
					}
				}
			}

			// For time-varying concentration-type BCs, interpolate the
			// boundary value from the time series at the current time.
			real Cauchybc = 0.0;
			if (current_bctype == SUB_RT_BC_Cauchy_T || current_bctype == SUB_RT_BC_Dirichlet_T)
			{
				if (current_ts.np > 0)
				{
					Cauchybc = interpolateLinear(current_ts, gdom.etime);
				}
				else
				{
					Cauchybc = 0.0;
				}
			}

			switch (current_bctype)
			{
			/* --------------------------------------------------------
			 * No-flow/no-transport boundary:
			 * apply constant normal extension to the complete exterior
			 * stencil.  Hydraulic no-flow implies no advective contribution,
			 * while merging the normal and cross-dispersion coefficients
			 * guarantees zero dispersive flux on all six face directions.
			 * -------------------------------------------------------- */
			case SUB_RT_BC_NOFLOW:
				Kokkos::parallel_for("rt_bc_noflow", ncellsBC,
					KOKKOS_CLASS_LAMBDA(int ibc) {
						int ii, jj, kk;
						const int iGlob = bcells[ibc];
						gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
						const int idom = (kk - gdom.hc) * gdom.nx * gdom.ny
						                + (jj - gdom.hc) * gdom.nx + ii - gdom.hc;
						applyConservativeFDElimination(rt, idom, direction);
					});
				break;

			/* --------------------------------------------------------
			 * Dirichlet BC (constant or time-varying):
			 * Enforced by exact row replacement.  The boundary cell is
			 * constrained to the prescribed ghost-cell concentration.
			 * This preserves the historical SUB_RT_BC_Dirichlet behavior used
			 * by the GW benchmark cases, where bctype=2 means a fixed
			 * concentration boundary rather than a flux-only Cauchy inlet.
			 * -------------------------------------------------------- */
			case SUB_RT_BC_Dirichlet_CONST:
			case SUB_RT_BC_Dirichlet_T:
				Kokkos::parallel_for("rt_bc_Dirichlet_exact", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc) {
	                        int ii, jj, kk, idom, iGlob = bcells[ibc], rtiGhost = rtgwcells[ibc];
	                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
	                        // Convert from halo-indexed global ID to domain-indexed (no halo) ID
	                        idom = (kk-gdom.hc)*gdom.nx*gdom.ny + (jj-gdom.hc)*gdom.nx + ii - gdom.hc;

							const real c_bc = rt.c(target_spec, rtiGhost, 1);
							for (int m = 0; m < 27; ++m)
								rt.RTcoef(idom, m) = 0.0;
							rt.RTcoef(idom, 0) = 1.0;
							rt.RTcoef(idom, 7) = c_bc;
						});
				break;

			/* --------------------------------------------------------
			 * Surface-water exchange BC (SUB_RT_BC_SWE): pure advection.
			 * This routine removes the top external stencil and adds the
			 * complete interface flux exactly once:
			 *
			 *   qss < 0: RHS += dt * (-qss) * Csw / dz
			 *   qss > 0: LHS += dt * qss / dz
			 *   |qss| ~ 0: no solute exchange
			 *
			 * No SWE normal dispersive flux is included.
			 * -------------------------------------------------------- */
			case SUB_RT_BC_SWE:

#if SERGHEI_SURFACE_TRANSPORT
				if (direction == 6)
				{
					Kokkos::parallel_for("rtgw_bc_swe_check", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc) {
	    int ii, jj, kk, idom, iGlob = bcells[ibc];
	    gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
	    idom = (kk - gdom.hc) * gdom.nx * gdom.ny + (jj - gdom.hc) * gdom.nx + ii - gdom.hc;
	    int iGlobSW2 = (jj - 1) * gdom.nx + ii - 1;

	    // 1. Read the exchange flux and check for NaN.
	    //    gw.qss is the surface-groundwater exchange rate [m/s].
	    real flux_swgw = gw.qss(iGlobSW2);

	    if (flux_swgw != flux_swgw) { // NaN check (NaN != NaN is true)
	        printf("ERROR: NaN detected in gw.qss at iGlob=%d! \n", iGlob);
	        flux_swgw = 0.0; // Force to zero to prevent propagation
	    }

	    // 2. Read the cell thickness and validate. [m]
	    real dz_val = gdom.dz(iGlob);
	    if (dz_val <= 1e-12) {
	        printf("ERROR: dz too small or zero at iGlob=%d! dz=%e\n", iGlob, dz_val);
	        dz_val = 1.0; // Prevent division by zero
	    }

							// Remove the exterior top-face contribution from the
							// diagonal and horizontal cross-dispersion corrections.
							// This is deliberately applied only to active SWE cells.
							const int plane = gdom.nxhc * gdom.nyhc;
							const int iGlobExterior = iGlob - plane;
							const real exterior_diag =
								rt.c_difzz(iGlobExterior) / dz_val
								- (1.0 - rt.Up_Weighting_zm(iGlob)) *
								  rt.c_advzz(iGlobExterior) / dz_val;
							rt.RTcoef(idom, 0) -= exterior_diag;

							const real lower_zx =
								(1.0 - rt.wz(iGlobExterior)) * rt.c_difzx(iGlobExterior) /
								(2.0 * gdom.dx * dz_val);
							const real lower_zy =
								(1.0 - rt.wz(iGlobExterior)) * rt.c_difzy(iGlobExterior) /
								(2.0 * gdom.dy * dz_val);
							rt.RTcoef(idom, 1) -= lower_zx;
							rt.RTcoef(idom, 2) += lower_zx;
							rt.RTcoef(idom, 3) -= lower_zy;
							rt.RTcoef(idom, 4) += lower_zy;

							// Do not let MPI_RHS_Update re-introduce a ghost-cell
							// concentration as a dispersive or advective source.
							clearRTCoef(rt, idom, 6);
							clearRTCoef(rt, idom, 17);
							clearRTCoef(rt, idom, 19);
							clearRTCoef(rt, idom, 23);
							clearRTCoef(rt, idom, 25);

							const real FLOW_EPSILON = 1.0e-10;
							if (flux_swgw < -FLOW_EPSILON)
							{
								// rt.csw is allocated on the surface grid including halo
								// cells.  Use the same halo index as the exchange audit.
								int iGlobSW = jj * gdom.nxhc + ii;
								real Csw = rt.csw(target_spec, iGlobSW);
								// Use the same rain fallback as applyConcentrationBC().
								if (gdom.isRain == 1 && gdom.rainRate(iGlobSW) > 0.0 && Csw == 0.0)
									Csw = rt.RainCon_arr(target_spec);

								// Infiltration source: |qss| * Csw.
								rt.RTcoef(idom, 7) +=
									rt.dt * (-flux_swgw) * Csw / dz_val;
							}
							else if (flux_swgw > FLOW_EPSILON)
							{
								// Exfiltration sink: qss * Cgw, implicit in Cgw.
								real sink_term = rt.dt * flux_swgw / dz_val;
								if (sink_term > 1e10 || sink_term != sink_term)
								{
									printf("WARNING: Unstable SWE sink at iGlob=%d. Flux=%e, dt=%e, dz=%e\n",
									       iGlob, flux_swgw, rt.dt, dz_val);
								}
								rt.RTcoef(idom, 0) += sink_term;
							}
						});
				}
				else
				{
					if (par.masterproc)
					{
						std::cerr << RERROR "RTGWBC direction must be 6 for rt-sw-gw boundary! " << "\n";
					}
				}
#endif
				break;

			/* --------------------------------------------------------
			 * Cauchy BC (constant or time-varying):
			 * Adds advective flux-proportional source or sink terms to
			 * the matrix.  The total flux across the boundary is:
			 *   F = v * C  (advective component)
			 * where v is the Darcy velocity normal to the face [m/s]
			 * and C is the concentration [mg/L].
			 *
			 * The sign convention for each direction determines whether
			 * positive velocity means inflow or outflow:
			 *
			 *   Direction | Inflow condition | Outflow condition
			 *   ----------|------------------|-------------------
			 *   1 (X+)    | v_x < 0          | v_x > 0
			 *   2 (X-)    | v_x > 0          | v_x < 0
			 *   3 (Y+)    | v_y < 0          | v_y > 0
			 *   4 (Y-)    | v_y > 0          | v_y < 0
			 *   5 (Z-)    | v_z > 0          | v_z < 0
			 *   6 (Z+)    | v_z < 0          | v_z > 0
			 *
			 * Inflow:  RHS (index 7) += dt * |v| * C_bc / cell_size
			 *          (known source, uses prescribed boundary concentration)
			 *
			 * Outflow: LHS diagonal (index 0) += dt * |v| / cell_size
			 *          (unknown sink, uses cell concentration - implicit)
			 * -------------------------------------------------------- */
			case SUB_RT_BC_Cauchy_CONST:
			case SUB_RT_BC_Cauchy_T:
			{
				int bcval_size = current_bcvals.size();
				Kokkos::parallel_for("rt_bc_Cauchy", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc) {
					int ii, jj, kk, idom, iGlob = bcells[ibc];
					gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
					idom = (kk - gdom.hc) * gdom.nx * gdom.ny + (jj - gdom.hc) * gdom.nx + ii - gdom.hc;

					// 1. Determine the boundary concentration value:
					//    constant/spatially varying (from array) or time-interpolated.
					real val_bc = (current_bctype == SUB_RT_BC_Cauchy_CONST)
					                 ? ((bcval_size > 1) ? current_bcvals(ibc) : current_bcvals(0))
					                 : Cauchybc;

					// 2. Get the cell thickness for flux discretization [m]
					real dz_val = gdom.dz(iGlob);

					// 3. Unit conversion factor (kept at 1.0 for SI consistency)
					real factor = 1.0;

					// =========================================================
					// Z direction boundaries (Direction 5: Bottom, Direction 6: Top)
					// =========================================================
					if (direction == 5 || direction == 6)
					{
						// Z-component of Darcy velocity [m/s].  Positive = upward.
						const int face = (direction == 5)
						                   ? iGlob
						                   : iGlob - gdom.nxhc * gdom.nyhc;
						real q_z = gw.q(face, 2);

						// --- Top boundary (Direction 6) ---
						if (direction == 6)
						{
							// Infiltration: water flows downward (q_z < 0).
							// Solute enters the system: known source -> RHS.
							// Uses water content (wc) as a saturation factor.
							if (q_z < 0)
							{
								// Mass source: |q_z| * C_boundary * factor / dz.
								// Multiplied by dt and water content for time integration.
								// Do NOT modify the LHS (index 0): this is a source term only.
								rt.RTcoef(idom, 7) += gw.wc(iGlob, 1) * rt.dt * (-q_z) * val_bc * factor / dz_val;
							}
							// Exfiltration / evaporation: water flows upward (q_z > 0).
							// Solute leaves the system: unknown sink -> LHS diagonal.
							// The outflow concentration equals the cell concentration.
							else if (q_z > 0)
							{
								// Implicit sink: q_z * dt / dz added to diagonal.
								// Do NOT modify the RHS (index 7): boundary concentration is irrelevant.
								rt.RTcoef(idom, 0) += rt.dt * (q_z)*factor / dz_val;
							}
						}
						// --- Bottom boundary (Direction 5) ---
						else if (direction == 5)
						{
							// Recharge from below: water flows upward (q_z > 0).
							// Known source -> RHS.
							if (q_z > 0)
							{
								rt.RTcoef(idom, 7) += rt.dt * (q_z)*val_bc * factor / dz_val;
							}
							// Leakage downward: water flows downward (q_z < 0).
							// Unknown sink -> LHS diagonal.
							else if (q_z < 0)
							{
								rt.RTcoef(idom, 0) += rt.dt * (-q_z) * factor / dz_val;
							}
						}
					}

					// =========================================================
					// X direction boundaries (Direction 1: Right, Direction 2: Left)
					// =========================================================
					else if (direction == 1 || direction == 2)
					{
						const int face = (direction == 1) ? iGlob : iGlob - 1;
						real q_x = gw.q(face, 0); // X-component of Darcy velocity [m/s].  Positive = rightward.
						real dx_val = gdom.dx;      // Cell width in X direction [m]

						// Right boundary (Direction 1): inflow if q < 0, outflow if q > 0
						if (direction == 1)
						{
							if (q_x < 0)
							{ // Inflow: known source -> RHS
								rt.RTcoef(idom, 7) += rt.dt * (-q_x) * val_bc * factor / dx_val;
							}
							else if (q_x > 0)
							{ // Outflow: unknown sink -> LHS diagonal
								rt.RTcoef(idom, 0) += rt.dt * (q_x)*factor / dx_val;
							}
						}
						// Left boundary (Direction 2): inflow if q > 0, outflow if q < 0
						else if (direction == 2)
						{
							if (q_x > 0)
							{ // Inflow: known source -> RHS
								rt.RTcoef(idom, 7) += rt.dt * (q_x)*val_bc * factor / dx_val;
							}
							else if (q_x < 0)
							{ // Outflow: unknown sink -> LHS diagonal
								rt.RTcoef(idom, 0) += rt.dt * (-q_x) * factor / dx_val;
							}
						}
					}

					// =========================================================
					// Y direction boundaries (Direction 3: Back/Ymax, Direction 4: Front/Ymin)
					// =========================================================
					else if (direction == 3 || direction == 4)
					{
						const int face = (direction == 3) ? iGlob : iGlob - gdom.nxhc;
						real q_y = gw.q(face, 1); // Y-component of Darcy velocity [m/s]
						real dy_val = gdom.dy;      // Cell width in Y direction [m]

						// Back boundary (Direction 3, Ymax): inflow if q < 0, outflow if q > 0
						if (direction == 3)
						{
							if (q_y < 0)
							{
								rt.RTcoef(idom, 7) += rt.dt * (-q_y) * val_bc * factor / dy_val;
							}
							else if (q_y > 0)
							{
								rt.RTcoef(idom, 0) += rt.dt * (q_y)*factor / dy_val;
							}
						}
						// Front boundary (Direction 4, Ymin): inflow if q > 0, outflow if q < 0
						else if (direction == 4)
						{
							if (q_y > 0)
							{
								rt.RTcoef(idom, 7) += rt.dt * (q_y)*val_bc * factor / dy_val;
							}
							else if (q_y < 0)
							{
								rt.RTcoef(idom, 0) += rt.dt * (-q_y) * factor / dy_val;
							}
						}
					} });
				break;
			}
			break;

			/* --------------------------------------------------------
			 * Conservative free-drainage BC (SUB_RT_BC_FD):
			 *
			 * Keep the assembled physical conservation equation. Enforce the
			 * zero-gradient condition by merging exterior normal ghost
			 * couplings into the corresponding interior stencil entries,
			 * instead of replacing the whole row with a concentration-copy
			 * constraint.
			 * -------------------------------------------------------- */
			case SUB_RT_BC_FD:
			{
				int bcval_size = current_bcvals.size();
				real fd_inflow_conc = 0.0;
				if (bcval_size > 0)
					fd_inflow_conc = current_bcvals(0);

				Kokkos::parallel_for("rt_bc_fd", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc) {
				    int ii, jj, kk, idom, iGlob = bcells[ibc];
				    gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
				    idom = (kk-gdom.hc)*gdom.nx*gdom.ny + (jj-gdom.hc)*gdom.nx + ii - gdom.hc;

					const int plane = gdom.nxhc * gdom.nyhc;
					real vn_out = 0.0;
					if (direction == XPLUS)
						vn_out = rt.aveVB(iGlob, 0);
					else if (direction == XMINUS)
						vn_out = -rt.aveVB(iGlob - 1, 0);
					else if (direction == YPLUS)
						vn_out = rt.aveVB(iGlob, 1);
					else if (direction == YMINUS)
						vn_out = -rt.aveVB(iGlob - gdom.nxhc, 1);
					else if (direction == ZPLUS)
						vn_out = rt.aveVB(iGlob, 2);
					else if (direction == ZMINUS)
						vn_out = -rt.aveVB(iGlob - plane, 2);

					const real c_external =
						(bcval_size > 1) ? current_bcvals(ibc) : fd_inflow_conc;

					if (vn_out >= 0.0)
						applyConservativeFDElimination(rt, idom, direction);
					else
						eliminateFDInflowGhost(rt, idom, direction, c_external);
				});
				break;
			}
			}
		}
	}
};

/* ==========================================================================
 * Section 4: RTSubsurfaceBoundaries Container Class
 * ==========================================================================
 * Holds all reactive-transport boundary condition objects for the
 * subsurface domain.  Each element of rtgwbc corresponds to one BC
 * region (one polygon on one face of the domain).
 *
 * Thread safety: This class should NOT be invoked from within a parallel
 * region because it contains std::string members (not thread-safe).
 * ========================================================================== */
class RTSubsurfaceBoundaries
{
	// This class should not be invoked from a parallel region as it contains strings
public:
	/* Boundary type name strings, indexed by the SUB_RT_BC_* integer values.
	 * Used for I/O and diagnostics.
	 *   Index 0: unused (BC types start at 1)
	 *   Index 1: NOFLOW         - No-flow / no-transport boundary
	 *   Index 2: CONST_H        - Constant head (hydraulic, not used in RT)
	 *   Index 3: CONST_Q        - Constant flux (hydraulic, not used in RT)
	 *   Index 4: CONST_WT       - Constant water table (hydraulic)
	 *   Index 5: H_TIMESERIES   - Time-series head
	 *   Index 6: Q_TIMESERIES   - Time-series flux
	 *   Index 7: WT_TIMESERIES  - Time-series water table
	 *   Index 8: SWEXCHANGE     - Surface-water exchange
	 *   Index 9: FREE_DRAINAGE  - Free drainage (zero-gradient outflow)
	 */
	std::string BoundaryTypes[9] = {"NOFLOW", "CONST_H", "CONST_Q", "CONST_WT", "H_TIMESERIES", "Q_TIMESERIES", "WT_TIMESERIES", "SWEXCHANGE", "FREE_DRAINAGE"};

	/* String identifiers for each boundary condition, matching the 'id'
	 * field in the input file.  Used for lookup and error reporting. */
	std::vector<std::string> id;

	/* Vector of RTBCGW objects, one per boundary condition region. */
	std::vector<RTBCGW> rtgwbc;
};
#endif

#endif
