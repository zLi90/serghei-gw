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
 *   ZPLUS  = top     face (maximum Z, land surface or uppermost layer)
 *   ZMINUS = bottom  face (minimum Z, bedrock or lowest layer)
 * ========================================================================== */

#define XPLUS 1	 /* Right  face (X maximum) */
#define XMINUS 2 /* Left   face (X minimum) */
#define YPLUS 3	 /* Back   face (Y maximum) */
#define YMINUS 4 /* Front  face (Y minimum) */
#define ZPLUS 5	 /* Top    face (Z maximum, near land surface) */
#define ZMINUS 6 /* Bottom face (Z minimum, near bedrock) */

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
	real QMassTot = 0, QMassInflow = 0, QMassOutflow = 0, QMassAbsflow = 0;

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
						// For the top boundary (ZPLUS), only the uppermost layer (kk == nz-1)
						// is selected.  The ghost cell is one layer above in the Z direction.
						if (direction == ZPLUS)
						{
							if (kk == gdom.nz - 1)
							{
								tmpbcells.push_back(iGlob);
								// Ghost cell for top boundary: one layer above iGlob in the
								// halo region (offset by nxhc * nyhc in global indexing).
								tmpgcells.push_back(iGlob + gdom.nxhc * gdom.nyhc);
							}
						}
						// For the bottom boundary (ZMINUS), only the lowest layer (kk == 0)
						// is selected.  The ghost cell is one layer below in the Z direction.
						else if (direction == ZMINUS)
						{
							if (kk == 0)
							{
								tmpbcells.push_back(iGlob);
								// Ghost cell for bottom boundary: one layer below iGlob.
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
		MPI_Allreduce(&ncellsBC, &ncells_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		// Gather which ranks have BC cells.
		// subdoms[i] >= 0 means rank i participates in this BC.
		int *subdoms;
		subdoms = (int *)malloc(par.nranks * sizeof(int));
		MPI_Allgather(&foundInSubdom, 1, MPI_INT, subdoms, 1, MPI_INT, MPI_COMM_WORLD);
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
		MPI_Comm_group(MPI_COMM_WORLD, &group);
		MPI_Group_incl(group, subdomains.size(), subdomains.data(), &subgroup);
		MPI_Comm_create(MPI_COMM_WORLD, subgroup, &comm);

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
	 * Also computes the total mass flux across the boundary for mass
	 * balance accounting.
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
	 *
	 * Returns:
	 *   void.  Results are stored in rt.c (concentrations) and in the
	 *   member variables QMassTot, QMassInflow, QMassOutflow [mg/s].
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
	inline void applyConcentrationBC(RTStateGW &rt, GwState &gw, GwDomain &gdom, Parallel &par, int target_spec = 0)
	{
		// Reset mass flux accumulators for this time step
		QMassTot = 0;
		QMassInflow = 0;
		QMassOutflow = 0;
		QMassAbsflow = 0;

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
			if (target_spec < spec_bctype.size())
			{
				current_bctype = spec_bctype[target_spec];
			}

			real Conbc;
			// For time-varying Dirichlet BC, interpolate the boundary
			// concentration from the species-specific time series.
			if (current_bctype == SUB_RT_BC_Dirichlet_T)
			{
				// Use the per-species time series if available
				if (target_spec < spec_ts.size() && spec_ts[target_spec].np > 0)
				{
					Conbc = interpolateLinear(spec_ts[target_spec], gdom.etime);
				}
				else
				{
					if (ncellsBC > 0)
					{
						std::cerr << RERROR "GWRTM: Boundary Dirichlet_T: No timeseries data for species " << target_spec << std::endl;
					}
					Conbc = 0.0;
				}
			}

			/* ----------------------------------------------------------
			 * Cauchy and Free-Drainage BCs: zero-gradient (Neumann) condition.
			 * The ghost-cell concentration is set equal to the interior-cell
			 * concentration, which enforces a zero concentration gradient at
			 * the boundary.  This prevents diffusive flux across the boundary
			 * while allowing advective outflow.
			 * ---------------------------------------------------------- */
			if (current_bctype == SUB_RT_BC_Cauchy_CONST || current_bctype == SUB_RT_BC_Cauchy_T || current_bctype == SUB_RT_BC_FD)
			{
				Kokkos::parallel_for("rt_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc) {
			                 int iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc];
				  rt.c(target_spec, rtiGhost, 1) = rt.c(target_spec, iGlob, 1); });
			}

			/* ----------------------------------------------------------
			 * Dirichlet constant BC: ghost-cell concentration is set to a
			 * prescribed value.  If multi-species spatial data exists, use
			 * the per-species array; otherwise use the global bcvals array.
			 * ---------------------------------------------------------- */
			else if (current_bctype == SUB_RT_BC_Dirichlet_CONST)
			{
				// Use per-species boundary values if available
				if (target_spec < spec_bcvals.size() && spec_bcvals[target_spec].size() > 0)
				{
					Kokkos::parallel_for("rt_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc) {
					                 int iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc];
					                 rt.c(target_spec, rtiGhost, 1) = spec_bcvals[target_spec](ibc); });
				}
				else
				{
					Kokkos::parallel_for("rt_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc) {
					                 int iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc];
					                 rt.c(target_spec, rtiGhost, 1) = bcvals(ibc); });
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
			 * This BC couples the subsurface transport model with the
			 * overland (surface) flow model.  The ghost-cell concentration
			 * is determined by the direction of the exchange flux:
			 *
			 *   Infiltration (surface -> subsurface, qss < 0):
			 *     Ghost cell receives the surface water concentration.
			 *     If the surface is dry and it is raining, the rain
			 *     concentration is used instead.
			 *
			 *   Exfiltration / stagnation (subsurface -> surface, qss >= 0):
			 *     Zero-gradient condition: ghost cell = interior cell,
			 *     preventing artificial dispersion from the surface into
			 *     the groundwater.
			 *
			 * Only valid for direction == ZMINUS (top boundary).
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

						// Determine the exchange direction based on the surface-
						// groundwater exchange flux gw.qss (positive = exfiltration,
						// negative = infiltration).  Units: [m/s].
						if (gw.qss(iGlobSW2) < -FLOW_EPSILON) // Significant infiltration (surface -> subsurface)
						{
							// Only read the surface concentration during genuine infiltration.
							// If the surface is dry but it is raining, use the rain concentration.
							if (gdom.isRain == 1 && rt.csw(target_spec, iGlobSW) == 0) {
								rt.c(target_spec, rtiGhost, 1) = rt.RainCon_arr(target_spec);
							} else {
								rt.c(target_spec, rtiGhost, 1) = rt.csw(target_spec, iGlobSW);
							}
						}
						else // Exfiltration (subsurface -> surface) or stagnation (qss ~ 0)
						{
							// Groundwater flows toward the surface, or velocity is negligible.
							// Apply zero-gradient condition (ghost = interior) to prevent
							// high surface concentrations from artificially dispersing
							// into the groundwater through the dispersive flux term.
							rt.c(target_spec, rtiGhost, 1) = rt.c(target_spec, iGlob, 1);
						}

						// Compute the solute mass exchange flux across the surface-
						// groundwater interface: Q_solute = qss * C_ghost.  [mg/(m^2*s)]
						rt.ConQss(target_spec, iGlobSW2) = gw.qss(iGlobSW2) * rt.c(target_spec, rtiGhost, 1); });
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
			 * Mass balance calculation: compute the total solute mass
			 * flux across the boundary for all species.  The flux is:
			 *   Q_mass = v * C * A
			 * where v is the Darcy velocity (aveVB) [m/s],
			 *       C is the concentration [mg/L or kg/m^3],
			 *       A is the cross-sectional area [m^2].
			 * The total mass flux QMassTot has units [mg/s or kg/s].
			 *
			 * The sign convention distinguishes inflow from outflow
			 * based on the boundary face direction:
			 *   - Right (X+): positive velocity = outflow
			 *   - Left  (X-): positive velocity = inflow
			 *   - Back  (Y+): positive velocity = outflow
			 *   - Front (Y-): positive velocity = inflow
			 *   - Top   (Z+): positive velocity = inflow (upward into domain)
			 *   - Bottom(Z-): positive velocity = outflow (downward out)
			 * ---------------------------------------------------------- */

			if (direction == 1) // Right boundary (X maximum)
			{
				// Use the X-component of Darcy velocity (aveVB index 0).
				// Cross-section area for an X-face = dy * dz.  [m^2]
				Kokkos::parallel_reduce("reducexRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
							int iGlob = bcells[ibc];
							const int iGhost = gcells[ibc];
							const int kkOff = gdom.nxhc * gdom.nyhc;
							const real area = gdom.dy * gdom.dz(iGlob);
							const real cCell = rt.c(target_spec, iGlob, 1);
							const real cGhost = rt.c(target_spec, iGhost, 1);
							const real dcdx = (cGhost - cCell) / gdom.dx;
							const real dcdy = (rt.c(target_spec, iGlob + gdom.nxhc, 1) - rt.c(target_spec, iGlob - gdom.nxhc, 1)) / (2.0 * gdom.dy);
							const real dcdz = (rt.c(target_spec, iGlob + kkOff, 1) - rt.c(target_spec, iGlob - kkOff, 1)) / (2.0 * gdom.dz(iGlob));
							// Per-species total boundary mass flux (advective + dispersive/diffusive).
							const real qAdv = rt.aveVB(iGlob, 0) * cCell * area / 1000.0; // preserve legacy x+ scaling
							const real qDisp = -(rt.dcal(target_spec, iGlob, 0) * dcdx +
												 rt.dcal(target_spec, iGlob, 3) * dcdy +
												 rt.dcal(target_spec, iGlob, 4) * dcdz) * area / 1000.0;
							tmp += qAdv + qDisp;
							}, Kokkos::Sum<real>(QMassTot));
				Kokkos::parallel_reduce("reducexRT_abs", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
							int iGlob = bcells[ibc];
							const int iGhost = gcells[ibc];
							const int kkOff = gdom.nxhc * gdom.nyhc;
							const real area = gdom.dy * gdom.dz(iGlob);
							const real cCell = rt.c(target_spec, iGlob, 1);
							const real cGhost = rt.c(target_spec, iGhost, 1);
							const real dcdx = (cGhost - cCell) / gdom.dx;
							const real dcdy = (rt.c(target_spec, iGlob + gdom.nxhc, 1) - rt.c(target_spec, iGlob - gdom.nxhc, 1)) / (2.0 * gdom.dy);
							const real dcdz = (rt.c(target_spec, iGlob + kkOff, 1) - rt.c(target_spec, iGlob - kkOff, 1)) / (2.0 * gdom.dz(iGlob));
							const real qAdv = rt.aveVB(iGlob, 0) * cCell * area / 1000.0;
							const real qDisp = -(rt.dcal(target_spec, iGlob, 0) * dcdx +
												 rt.dcal(target_spec, iGlob, 3) * dcdy +
												 rt.dcal(target_spec, iGlob, 4) * dcdz) * area / 1000.0;
							real qcell = qAdv + qDisp;
							tmp += myfabs(qcell); }, Kokkos::Sum<real>(QMassAbsflow));
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
				Kokkos::parallel_reduce("reducexRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
							int iGlob = bcells[ibc];
							const int iGhost = gcells[ibc];
							const int kkOff = gdom.nxhc * gdom.nyhc;
							const real area = gdom.dy * gdom.dz(iGlob);
							const real cCell = rt.c(target_spec, iGlob, 1);
							const real cGhost = rt.c(target_spec, iGhost, 1);
							const real dcdx = (cCell - cGhost) / gdom.dx;
							const real dcdy = (rt.c(target_spec, iGlob + gdom.nxhc, 1) - rt.c(target_spec, iGlob - gdom.nxhc, 1)) / (2.0 * gdom.dy);
							const real dcdz = (rt.c(target_spec, iGlob + kkOff, 1) - rt.c(target_spec, iGlob - kkOff, 1)) / (2.0 * gdom.dz(iGlob));
							const real qAdv = rt.aveVB(iGlob, 0) * cCell * area;
							const real qDisp = -(rt.dcal(target_spec, iGlob, 0) * dcdx +
												 rt.dcal(target_spec, iGlob, 3) * dcdy +
												 rt.dcal(target_spec, iGlob, 4) * dcdz) * area;
							tmp += qAdv + qDisp;
							}, Kokkos::Sum<real>(QMassTot));
				Kokkos::parallel_reduce("reducexRT_abs", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
							int iGlob = bcells[ibc];
							const int iGhost = gcells[ibc];
							const int kkOff = gdom.nxhc * gdom.nyhc;
							const real area = gdom.dy * gdom.dz(iGlob);
							const real cCell = rt.c(target_spec, iGlob, 1);
							const real cGhost = rt.c(target_spec, iGhost, 1);
							const real dcdx = (cCell - cGhost) / gdom.dx;
							const real dcdy = (rt.c(target_spec, iGlob + gdom.nxhc, 1) - rt.c(target_spec, iGlob - gdom.nxhc, 1)) / (2.0 * gdom.dy);
							const real dcdz = (rt.c(target_spec, iGlob + kkOff, 1) - rt.c(target_spec, iGlob - kkOff, 1)) / (2.0 * gdom.dz(iGlob));
							const real qAdv = rt.aveVB(iGlob, 0) * cCell * area;
							const real qDisp = -(rt.dcal(target_spec, iGlob, 0) * dcdx +
												 rt.dcal(target_spec, iGlob, 3) * dcdy +
												 rt.dcal(target_spec, iGlob, 4) * dcdz) * area;
							real qcell = qAdv + qDisp;
							tmp += myfabs(qcell); }, Kokkos::Sum<real>(QMassAbsflow));
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
				Kokkos::parallel_reduce("reduceyRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
							int iGlob = bcells[ibc];
							const int iGhost = gcells[ibc];
							const int kkOff = gdom.nxhc * gdom.nyhc;
							const real area = gdom.dx * gdom.dz(iGlob);
							const real cCell = rt.c(target_spec, iGlob, 1);
							const real cGhost = rt.c(target_spec, iGhost, 1);
							const real dcdx = (rt.c(target_spec, iGlob + 1, 1) - rt.c(target_spec, iGlob - 1, 1)) / (2.0 * gdom.dx);
							const real dcdy = (cGhost - cCell) / gdom.dy;
							const real dcdz = (rt.c(target_spec, iGlob + kkOff, 1) - rt.c(target_spec, iGlob - kkOff, 1)) / (2.0 * gdom.dz(iGlob));
							const real qAdv = rt.aveVB(iGlob, 1) * cCell * area;
							const real qDisp = -(rt.dcal(target_spec, iGlob, 5) * dcdx +
												 rt.dcal(target_spec, iGlob, 1) * dcdy +
												 rt.dcal(target_spec, iGlob, 6) * dcdz) * area;
							tmp += qAdv + qDisp;
							}, Kokkos::Sum<real>(QMassTot));
				Kokkos::parallel_reduce("reduceyRT_abs", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
							int iGlob = bcells[ibc];
							const int iGhost = gcells[ibc];
							const int kkOff = gdom.nxhc * gdom.nyhc;
							const real area = gdom.dx * gdom.dz(iGlob);
							const real cCell = rt.c(target_spec, iGlob, 1);
							const real cGhost = rt.c(target_spec, iGhost, 1);
							const real dcdx = (rt.c(target_spec, iGlob + 1, 1) - rt.c(target_spec, iGlob - 1, 1)) / (2.0 * gdom.dx);
							const real dcdy = (cGhost - cCell) / gdom.dy;
							const real dcdz = (rt.c(target_spec, iGlob + kkOff, 1) - rt.c(target_spec, iGlob - kkOff, 1)) / (2.0 * gdom.dz(iGlob));
							const real qAdv = rt.aveVB(iGlob, 1) * cCell * area;
							const real qDisp = -(rt.dcal(target_spec, iGlob, 5) * dcdx +
												 rt.dcal(target_spec, iGlob, 1) * dcdy +
												 rt.dcal(target_spec, iGlob, 6) * dcdz) * area;
							real qcell = qAdv + qDisp;
							tmp += myfabs(qcell); }, Kokkos::Sum<real>(QMassAbsflow));
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
				Kokkos::parallel_reduce("reduceyRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
							int iGlob = bcells[ibc];
							const int iGhost = gcells[ibc];
							const int kkOff = gdom.nxhc * gdom.nyhc;
							const real area = gdom.dx * gdom.dz(iGlob);
							const real cCell = rt.c(target_spec, iGlob, 1);
							const real cGhost = rt.c(target_spec, iGhost, 1);
							const real dcdx = (rt.c(target_spec, iGlob + 1, 1) - rt.c(target_spec, iGlob - 1, 1)) / (2.0 * gdom.dx);
							const real dcdy = (cCell - cGhost) / gdom.dy;
							const real dcdz = (rt.c(target_spec, iGlob + kkOff, 1) - rt.c(target_spec, iGlob - kkOff, 1)) / (2.0 * gdom.dz(iGlob));
							const real qAdv = rt.aveVB(iGlob, 1) * cCell * area;
							const real qDisp = -(rt.dcal(target_spec, iGlob, 5) * dcdx +
												 rt.dcal(target_spec, iGlob, 1) * dcdy +
												 rt.dcal(target_spec, iGlob, 6) * dcdz) * area;
							tmp += qAdv + qDisp;
							}, Kokkos::Sum<real>(QMassTot));
				Kokkos::parallel_reduce("reduceyRT_abs", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
							int iGlob = bcells[ibc];
							const int iGhost = gcells[ibc];
							const int kkOff = gdom.nxhc * gdom.nyhc;
							const real area = gdom.dx * gdom.dz(iGlob);
							const real cCell = rt.c(target_spec, iGlob, 1);
							const real cGhost = rt.c(target_spec, iGhost, 1);
							const real dcdx = (rt.c(target_spec, iGlob + 1, 1) - rt.c(target_spec, iGlob - 1, 1)) / (2.0 * gdom.dx);
							const real dcdy = (cCell - cGhost) / gdom.dy;
							const real dcdz = (rt.c(target_spec, iGlob + kkOff, 1) - rt.c(target_spec, iGlob - kkOff, 1)) / (2.0 * gdom.dz(iGlob));
							const real qAdv = rt.aveVB(iGlob, 1) * cCell * area;
							const real qDisp = -(rt.dcal(target_spec, iGlob, 5) * dcdx +
												 rt.dcal(target_spec, iGlob, 1) * dcdy +
												 rt.dcal(target_spec, iGlob, 6) * dcdz) * area;
							real qcell = qAdv + qDisp;
							tmp += myfabs(qcell); }, Kokkos::Sum<real>(QMassAbsflow));
				if (QMassTot > 0)
				{
					QMassInflow = QMassTot;
				}
				else
				{
					QMassOutflow = -QMassTot;
				}
			}

			else if (direction == 5) // Bottom boundary (Z minimum)
			{
				// Use the Z-component of Darcy velocity (aveVB index 2).
				// Cross-section area for a Z-face = dx * dy.  [m^2]
				Kokkos::parallel_reduce("reducezRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
							int iGlob = bcells[ibc];
							const int iGhost = gcells[ibc];
							const real area = gdom.dx * gdom.dy;
							const real cCell = rt.c(target_spec, iGlob, 1);
							const real cGhost = rt.c(target_spec, iGhost, 1);
							const real dcdx = (rt.c(target_spec, iGlob + 1, 1) - rt.c(target_spec, iGlob - 1, 1)) / (2.0 * gdom.dx);
							const real dcdy = (rt.c(target_spec, iGlob + gdom.nxhc, 1) - rt.c(target_spec, iGlob - gdom.nxhc, 1)) / (2.0 * gdom.dy);
							const real dcdz = (cGhost - cCell) / gdom.dz(iGlob);
							const real qAdv = rt.aveVB(iGlob, 2) * cCell * area;
							const real qDisp = -(rt.dcal(target_spec, iGlob, 7) * dcdx +
												 rt.dcal(target_spec, iGlob, 8) * dcdy +
												 rt.dcal(target_spec, iGlob, 2) * dcdz) * area;
							tmp += qAdv + qDisp;
							}, Kokkos::Sum<real>(QMassTot));
				Kokkos::parallel_reduce("reducezRT_abs", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
							int iGlob = bcells[ibc];
							const int iGhost = gcells[ibc];
							const real area = gdom.dx * gdom.dy;
							const real cCell = rt.c(target_spec, iGlob, 1);
							const real cGhost = rt.c(target_spec, iGhost, 1);
							const real dcdx = (rt.c(target_spec, iGlob + 1, 1) - rt.c(target_spec, iGlob - 1, 1)) / (2.0 * gdom.dx);
							const real dcdy = (rt.c(target_spec, iGlob + gdom.nxhc, 1) - rt.c(target_spec, iGlob - gdom.nxhc, 1)) / (2.0 * gdom.dy);
							const real dcdz = (cGhost - cCell) / gdom.dz(iGlob);
							const real qAdv = rt.aveVB(iGlob, 2) * cCell * area;
							const real qDisp = -(rt.dcal(target_spec, iGlob, 7) * dcdx +
												 rt.dcal(target_spec, iGlob, 8) * dcdy +
												 rt.dcal(target_spec, iGlob, 2) * dcdz) * area;
							real qcell = qAdv + qDisp;
							tmp += myfabs(qcell); }, Kokkos::Sum<real>(QMassAbsflow));
				if (QMassTot > 0)
				{
					QMassOutflow = QMassTot;
				}
				else
				{
					QMassInflow = -QMassTot;
				}
			}
			else if (direction == 6) // Top boundary (Z maximum, land surface)
			{
				Kokkos::parallel_reduce("reducezRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
							int iGlob = bcells[ibc];
							const int iGhost = gcells[ibc];
							const real area = gdom.dx * gdom.dy;
							const real cCell = rt.c(target_spec, iGlob, 1);
							const real cGhost = rt.c(target_spec, iGhost, 1);
							const real dcdx = (rt.c(target_spec, iGlob + 1, 1) - rt.c(target_spec, iGlob - 1, 1)) / (2.0 * gdom.dx);
							const real dcdy = (rt.c(target_spec, iGlob + gdom.nxhc, 1) - rt.c(target_spec, iGlob - gdom.nxhc, 1)) / (2.0 * gdom.dy);
							const real dcdz = (cCell - cGhost) / gdom.dz(iGlob);
							const real qAdv = rt.aveVB(iGlob, 2) * cCell * area;
							const real qDisp = -(rt.dcal(target_spec, iGlob, 7) * dcdx +
												 rt.dcal(target_spec, iGlob, 8) * dcdy +
												 rt.dcal(target_spec, iGlob, 2) * dcdz) * area;
							tmp += qAdv + qDisp;
							}, Kokkos::Sum<real>(QMassTot));
				Kokkos::parallel_reduce("reducezRT_abs", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
							int iGlob = bcells[ibc];
							const int iGhost = gcells[ibc];
							const real area = gdom.dx * gdom.dy;
							const real cCell = rt.c(target_spec, iGlob, 1);
							const real cGhost = rt.c(target_spec, iGhost, 1);
							const real dcdx = (rt.c(target_spec, iGlob + 1, 1) - rt.c(target_spec, iGlob - 1, 1)) / (2.0 * gdom.dx);
							const real dcdy = (rt.c(target_spec, iGlob + gdom.nxhc, 1) - rt.c(target_spec, iGlob - gdom.nxhc, 1)) / (2.0 * gdom.dy);
							const real dcdz = (cCell - cGhost) / gdom.dz(iGlob);
							const real qAdv = rt.aveVB(iGlob, 2) * cCell * area;
							const real qDisp = -(rt.dcal(target_spec, iGlob, 7) * dcdx +
												 rt.dcal(target_spec, iGlob, 8) * dcdy +
												 rt.dcal(target_spec, iGlob, 2) * dcdz) * area;
							real qcell = qAdv + qDisp;
							tmp += myfabs(qcell); }, Kokkos::Sum<real>(QMassAbsflow));
				if (QMassTot > 0)
				{
					QMassInflow = QMassTot;
				}
				else
				{
					QMassOutflow = -QMassTot;
				}
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
	 *     Enforces a prescribed concentration using a large penalty
	 *     (big-number) method: the diagonal coefficient is multiplied
	 *     by ~1e50 and the RHS is set to 1e50 * diag * C_boundary.
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
	 *     Zero-gradient condition: resets all coefficients and enforces
	 *     C_cell = C_neighbor in the outflow direction (upwind).
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

			// For time-varying Cauchy BC, interpolate the boundary value
			// from the time series at the current simulation time.
			real Cauchybc;
			if (current_bctype == SUB_RT_BC_Cauchy_T)
			{
				if (target_spec < spec_ts.size() && spec_ts[target_spec].np > 0)
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
			 * Dirichlet BC (constant or time-varying):
			 * Enforced via the "big number" (penalty) method.
			 * The diagonal coefficient is scaled by ~1e50, and the
			 * RHS is set to 1e50 * diag * C_boundary.  This makes
			 * the solution at the boundary cell approximately equal
			 * to the prescribed concentration C_boundary.
			 * -------------------------------------------------------- */
			case SUB_RT_BC_Dirichlet_CONST:
			case SUB_RT_BC_Dirichlet_T:
				Kokkos::parallel_for("rt_bc_Dirichlet_const", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc) {
	                        int ii, jj, kk, idom,  iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc];
	                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
	                        // Convert from halo-indexed global ID to domain-indexed (no halo) ID
	                        idom = (kk-gdom.hc)*gdom.nx*gdom.ny + (jj-gdom.hc)*gdom.nx + ii - gdom.hc;

							if (direction == 1) {  // Right boundary (X+)
								// Penalty method: ghost cell is at iGlob+1 (right neighbor)
								rt.RTcoef(idom,7) = 1e50 * rt.RTcoef(idom,0) * rt.c(target_spec, iGlob+1,1);
								rt.RTcoef(idom,0) = 1e50 * rt.RTcoef(idom,0);
							}
							else if (direction == 2) {  // Left boundary (X-)
								// Variable water-table Dirichlet BC:
								// If the cell is within the saturated zone (head >= 0),
								// the BC behavior depends on the flow direction:
								//   - Inflow (q <= 0): apply Dirichlet penalty to enforce
								//     the prescribed concentration from the boundary.
								//   - Outflow (q > 0): zero-gradient condition (free outflow
								//     of solute from the aquifer).
								// If the cell is in the unsaturated zone (head < 0),
								// a seepage-face condition is applied: solute freely
								// drains out (zero-gradient).

								if (gw.h(iGlob, 1) >= 0) { // Saturated zone: cell is within the aquifer
									if (gw.q(iGlob, 0) <= 0) { // Solute flows into the aquifer (inflow)
										rt.RTcoef(idom,7) = 1e50 * rt.RTcoef(idom,0) * rt.c(target_spec, iGlob-1,1);
										rt.RTcoef(idom,0) = 1e50 * rt.RTcoef(idom,0);
									}
									else { // Aquifer head > river head: solute freely flows out
										for (int i = 0; i < 30; i++) {
											rt.RTcoef(idom, i) = 0.0;
										}
										rt.RTcoef(idom, 0) = 1.0;
										rt.RTcoef(idom, 1) = -1; // C_{i+1}: zero-gradient (C = C_{i+1})
									}
								}
								else { // Unsaturated zone: seepage face, solute freely drains out
									for (int i = 0; i < 30; i++) {
										rt.RTcoef(idom, i) = 0.0;
									}
									rt.RTcoef(idom, 0) = 1.0;
									rt.RTcoef(idom, 1) = -1; // C_{i+1}: zero-gradient
								}
							}
							else if (direction == 3) {  // Back boundary (Y+)
								rt.RTcoef(idom,7) = 1e50 * rt.RTcoef(idom,0) * rt.c(target_spec, iGlob+gdom.nxhc,1);
								rt.RTcoef(idom,0) = 1e50 * rt.RTcoef(idom,0);
							}
							else if (direction == 4) {  // Front boundary (Y-)
								rt.RTcoef(idom,7) = 1e50 * rt.RTcoef(idom,0) * rt.c(target_spec, iGlob-gdom.nxhc,1);
								rt.RTcoef(idom,0) = 1e50 * rt.RTcoef(idom,0);
							}
							else if (direction == 5) {  // Bottom boundary (Z-)
								rt.RTcoef(idom,7) = 1e50 * rt.RTcoef(idom,0) * rt.c(target_spec, iGlob+gdom.nxhc*gdom.nyhc,1);
								rt.RTcoef(idom,0) = 1e50 * rt.RTcoef(idom,0);
							}
							else if (direction == 6) {  // Top boundary (Z+)
								rt.RTcoef(idom,7) = 1e50 * rt.RTcoef(idom,0) * rt.c(target_spec, iGlob-gdom.nxhc*gdom.nyhc,1);
								rt.RTcoef(idom,0) = 1e50 * rt.RTcoef(idom,0);
							} });
				break;

			/* --------------------------------------------------------
			 * Surface-water exchange BC (SUB_RT_BC_SWE):
			 * Modifies the matrix coefficients to account for mass
			 * exchange between surface water and groundwater.
			 *
			 *   Infiltration (qss < 0): mass enters the subsurface.
			 *     The solute flux is added to the RHS (index 7) as:
			 *       RHS += dt * (-ConQss) / dz
			 *     where ConQss = qss * C_ghost [mg/(m^2*s)].
			 *
			 *   Exfiltration (qss > 0): mass leaves the subsurface.
			 *     The sink term is added to the LHS diagonal (index 0):
			 *       diag += dt * qss / dz
			 *     This is treated implicitly for numerical stability.
			 *
			 *   NaN checks are performed to prevent numerical blow-up
			 *   from corrupted flux values.
			 *
			 * Only valid for direction == ZMINUS (top boundary).
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

	    // 3. Apply the exchange terms to the matrix coefficients
	    if (flux_swgw <= 0.0) // Infiltration (surface -> subsurface): source term
	    {
	        // Read the solute mass exchange flux and validate.
	        // ConQss = qss * C_ghost [mg/(m^2*s)]
	        real mass_flux = rt.ConQss(target_spec, iGlobSW2);
	        if (mass_flux != mass_flux) { // NaN check
	             printf("ERROR: NaN in ConQss at iGlobSW2=%d\n", iGlobSW2);
	             mass_flux = 0.0;
	        }
	        // Add the mass source to the RHS (index 7).
	        // The negative sign converts infiltration flux direction.
	        rt.RTcoef(idom, 7) += rt.dt * (-1.0 * mass_flux) / dz_val;
	    }
	    else // Exfiltration (subsurface -> surface): sink term
	    {
	        // Add the mass sink to the LHS diagonal (index 0) implicitly.
	        // sink_term = dt * qss / dz  [dimensionless coefficient increment]
	        real sink_term = rt.dt * (flux_swgw) / dz_val;

	        // 4. Check for abnormally large coefficients (potential instability)
	        if (sink_term > 1e10 || sink_term != sink_term) {
	             printf("WARNING: Unstable sink term at iGlob=%d. Flux=%e, dt=%e, dz=%e\n",
	                    iGlob, flux_swgw, rt.dt, dz_val);
	        }

	        rt.RTcoef(idom, 0) += sink_term;
	    } });
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
				Kokkos::parallel_for("rt_bc_Cauchy", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc) {
					int ii, jj, kk, idom, iGlob = bcells[ibc];
					gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
					idom = (kk - gdom.hc) * gdom.nx * gdom.ny + (jj - gdom.hc) * gdom.nx + ii - gdom.hc;

					// 1. Determine the boundary concentration value:
					//    constant (from array) or time-interpolated.
					real val_bc = (current_bctype == SUB_RT_BC_Cauchy_CONST) ? current_bcvals(0) : Cauchybc;

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
						real q_z = gw.q(iGlob, 2);

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
						real q_x = gw.q(iGlob, 0); // X-component of Darcy velocity [m/s].  Positive = rightward.
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
						real q_y = gw.q(iGlob, 1); // Y-component of Darcy velocity [m/s]
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
			 * Free drainage BC (SUB_RT_BC_FD):
			 * Zero-gradient (upwind) outflow condition.
			 *
			 * All 30 matrix coefficients are zeroed, then a simple
			 * equation is set up: C_cell = C_upwind_neighbor, i.e.,
			 *   RTcoef(idom, 0) = 1.0  (diagonal)
			 *   RTcoef(idom, neighbor_idx) = -1.0  (upwind neighbor)
			 *
			 * The upwind neighbor index depends on the boundary face:
			 *   Direction 1 (X+): neighbor is index 2 (C_{i-1}, left neighbor)
			 *   Direction 2 (X-): neighbor is index 1 (C_{i+1}, right neighbor)
			 *   Direction 3 (Y+): neighbor is index 4 (C_{j-1}, front neighbor)
			 *   Direction 4 (Y-): neighbor is index 3 (C_{j+1}, back neighbor)
			 *   Direction 5 (Z-): neighbor is index 6 (C_{k-1}, lower neighbor)
			 *   Direction 6 (Z+): neighbor is index 5 (C_{k+1}, upper neighbor)
			 * -------------------------------------------------------- */
			case SUB_RT_BC_FD:
				Kokkos::parallel_for("rt_bc_fd", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc) {
				    int ii, jj, kk, idom, iGlobSW, iGlob = bcells[ibc];
				    gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
				    idom = (kk-gdom.hc)*gdom.nx*gdom.ny + (jj-gdom.hc)*gdom.nx + ii - gdom.hc;
				    iGlobSW = jj*gdom.nxhc + ii;

				    // Determine upwind weighting factors based on the local
				    // flow velocity direction in each coordinate.  These are
				    // currently computed but not directly used in the matrix
				    // modification below (reserved for future upwind schemes).
				    double Up_Weighting_x, Up_Weighting_y, Up_Weighting_z;
				    if (rt.aveVB(iGlob, 0) > 0)
				    {
				    	Up_Weighting_x = rt.Up_Weighting_vplus;
				    }
				    else
				    {
				    	Up_Weighting_x = rt.Up_Weighting_vminus;
				    }

				    // Y-direction velocity check
				    if (rt.aveVB(iGlob, 1) > 0)
				    {
				    	Up_Weighting_y = rt.Up_Weighting_vplus;
				    }
				    else
				    {
				    	Up_Weighting_y = rt.Up_Weighting_vminus;
				    }

				    // Z-direction velocity check
				    if (rt.aveVB(iGlob, 2) > 0)
				    {
				    	Up_Weighting_z = rt.Up_Weighting_vplus;
				    }
				    else
				    {
				    	Up_Weighting_z = rt.Up_Weighting_vminus;
				    }

				    // Apply zero-gradient condition for the appropriate face.
				    // Reset all 30 coefficients, then set diagonal = 1 and
				    // the upwind neighbor coefficient = -1.
					if (direction == 1) {   // Right face (X+): C = C_{i-1}
						for (int i = 0; i < 30; i++) {
							rt.RTcoef(idom, i) = 0.0;
						}
						rt.RTcoef(idom, 0) = 1.0;
						rt.RTcoef(idom, 2) = -1; // C_{i-1}: left neighbor
					}
					else if (direction == 2) { // Left face (X-): C = C_{i+1}
						for (int i = 0; i < 30; i++) {
							rt.RTcoef(idom, i) = 0.0;
						}
						rt.RTcoef(idom, 0) = 1.0;
						rt.RTcoef(idom, 1) = -1; // C_{i+1}: right neighbor
					}
					else if (direction == 3) { // Back face (Y+): C = C_{j-1}
						for (int i = 0; i < 30; i++) {
							rt.RTcoef(idom, i) = 0.0;
						}
						rt.RTcoef(idom, 0) = 1.0;
						rt.RTcoef(idom, 4) = -1; // C_{j-1}: front neighbor
					}
					else if (direction == 4) { // Front face (Y-): C = C_{j+1}
						for (int i = 0; i < 30; i++) {
							rt.RTcoef(idom, i) = 0.0;
						}
						rt.RTcoef(idom, 0) = 1.0;
						rt.RTcoef(idom, 3) = -1; // C_{j+1}: back neighbor
					}
					else if (direction == 5) { // Bottom face (Z-): C = C_{k-1}
						for (int i = 0; i < 30; i++) {
							rt.RTcoef(idom, i) = 0.0;
						}
						rt.RTcoef(idom, 0) = 1.0;
						rt.RTcoef(idom, 6) = -1; // C_{k-1}: lower neighbor
					}
					else if (direction == 6) { // Top face (Z+): C = C_{k+1}
						rt.RTcoef(idom, 6) = 0.0;
						for (int i = 0; i < 30; i++) {
							rt.RTcoef(idom, i) = 0.0;
						}
						rt.RTcoef(idom, 0) = 1.0;
						rt.RTcoef(idom, 5) = -1;  // C_{k+1}: upper neighbor
					} });
				break;
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
