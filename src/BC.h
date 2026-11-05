/* -*- mode: c++ -*- */

#ifndef _BC_H_
#define _BC_H_

#include "define.h"
#include "Indexing.h"
#if CROP_GROWTH_MODEL
#include "./cropsrc/Wofost72.h"
#endif
#if CROP_GROWTH_MODEL && SERGHEI_RE_MODEL
#include "SourceSinkCrop.h" // 根系层土壤水触发：单柱 helper rootZoneMoistureSums
#endif

// outer boundary direction definitions
#define SERGHEI_BC_OUTER_REFLECTIVE_CELL 9999
#if SERGHEI_SURFACE_TRANSPORT
#include "RTStateSW.h"
// ==========================================
// Surface solute transport boundary conditions
// ==========================================
#define RTSW_BC_NONE 1		// no intervention
#define RTSW_BC_CON_CONST 2 // constant concentration (steady inflow)
#define RTSW_BC_ZEROGRAD 3	// zero concentration gradient (free outflow)
#define RTSW_BC_CON_T 4		// concentration time series (dynamic inflow)
#define RTSW_BC_TRIGGERED_IRRIGATION 5	// triggered irrigation: Dirichlet while watering, zero-gradient otherwise
#define RTSW_BC_TRIGGERED_DRAINAGE 6	// triggered drainage: zero-gradient outflow (explicit type)
#endif
// these definitions are meant for hydraulics
#define SWE_BC_PERIODIC 1
#define SWE_BC_REFLECTIVE 2
#define SWE_BC_TRANSMISSIVE 3
#define SWE_BC_CRITICAL 5
#define SWE_BC_H_CONST 6
#define SWE_BC_Q_CONST 7
#define SWE_BC_WSE_CONST 8
#define SWE_BC_FREE_OUTFLOW 9
#define SWE_BC_HZ_T_INLET 10
#define SWE_BC_HZ_T_OUTLET 11
#define SWE_BC_Q_T 12
#define SWE_BC_TRIGGERED_IRRIGATION 13
#define SWE_BC_TRIGGERED_DRAINAGE 14

KOKKOS_INLINE_FUNCTION real criticalDepth(real hu, real hv, real Fr)
{
	return (cbrt((hu * hu + hv * hv) / (GRAV * Fr * Fr)));
}

class ExtBC
{
	// this class is safe to invoke in a parallel region
public:
	int ncellsBC = 0;	   // number of bcells
	int ncellsTrigger = 0; // number of trigger-observation cells
	intArr bcells;		   // array of indexes of boundary cells
	intArr triggerCells;   // array of indexes used to compute trigger averages
	real normalx, normaly; // direction set by user for inflow/outflow
	int location;		   // 1->west, 2->north, 3->east, 4-> south
	int bctype = 0;
	int hydro_bctype = 0;  // hydrodynamic BC type copied into RT BC objects
	int hydroTriggerActive = 0; // runtime on/off state synced from the paired hydro BC each step
	int isInDomain;
	realArr bcvals;
	real outflowDischarge;
	real outflowAccumulated = 0;
	real inflowDischarge;
	real inflowAccumulated = 0;
	real adjustedVolume = 0;
		int irrigationActive = 0;
		int rtIrrigationFallbackWarned = 0;
		int irrigationStrategyWarned = 0;
		real irrigationMeanH = 0.0;
		real irrigationMeanGwDepth = -1.0;
		real irrigationAppliedDischarge = 0.0;
		std::string irrigationStrategyFile;
		int useDVSForIrrigation = 0;
		int nIrrigationStages = 0;
		int currentIrrigationStage = -1;
		real cropStartTime = 0.0;
		real irrigationDAT = 0.0;
		real irrigationMeanDVS = 0.0;
		realArr irrigationDatMin;
		realArr irrigationDatMax;
		realArr irrigationDvsMin;
		realArr irrigationDvsMax;
		realArr irrigationHMin;
		realArr irrigationHMax;
		realArr irrigationGwDepthTrigger;
		realArr irrigationQ;
		intArr irrigationEnabled;
		intArr drainageEnabled;
		intArr irrigationTrigType;  // 每阶段触发类型：0=沿用 bcvals(0) 全局模式；1=地表水深；2=根系层土壤水；3=任一；4=同时
		realArr irrigationThetaTrig; // 触发阈值（根系层相对饱和度分数，如 0.70）；-1=未配置
		realArr irrigationThetaStop; // 提前停阈值（相对饱和度）；-1=不启用（仅按 h_max 停）
		real irrigationMeanTheta = -1.0; // 触发单元上的根系层平均相对饱和度（θ kernel 每步刷新）
		int irrigationThetaWarned = 0;
#if CROP_GROWTH_MODEL && SERGHEI_RE_MODEL
		GwState *gw = nullptr;	// 一次性 wiring（serghei.h 初始化，先例 GwIntegrator::initialize），土壤水触发读 wc/vgTable/soilID
		GwDomain *gdom = nullptr; // 同上，读 dz/depth 与网格维度
#endif
		int drainageActive = 0;
		int rtDrainageFallbackWarned = 0;
		int drainageStrategyWarned = 0;
		real drainageMeanH = 0.0;
		real drainageAppliedDischarge = 0.0;
		std::string drainageStrategyFile;
		int useDVSForDrainage = 0;
		int nDrainageStages = 0;
		int currentDrainageStage = -1;
		real drainageDAT = 0.0;
		real drainageMeanDVS = 0.0;
		realArr drainageDatMin;
		realArr drainageDatMax;
		realArr drainageDvsMin;
		realArr drainageDvsMax;
		realArr drainageOnH;
		realArr drainageOffH;
		intArr drainageMode;
		realArr drainageQ;
		intArr drainageStageEnabled;

#if SERGHEI_SUSPENDED_SEDIMENT
	real outflowSolidDischarge;
	real outflowSolidAccumulated = 0;
	real inflowSolidDischarge;
	real inflowSolidAccumulated = 0;
	real adjustedSolidVolume = 0;
#endif

	TimeSeries hydrograph;
	real netQ, netVol;

	real hzMin = 1E6; // lowest water surface in boundary cross section
	real zMin = 1E6;  // lowest bed elevation in boundary cross section
	real zMax = -1E6;
	int nzMin = 0;

	MPI_Comm comm; // communicator for ranks associated to the BC

	// Multi-species boundary condition support
	std::vector<int> spec_bctype;		  // boundary type per species
	std::vector<realArr> spec_bcvals;	  // spatially distributed boundary values per species
	std::vector<std::string> spec_bcfile; // boundary file name per species
	std::vector<TimeSeries> spec_ts;	  // time series per species
	std::vector<real> spec_bcval_const;	  // constant boundary value per species
	std::vector<int> has_bcfile;		  // whether a boundary file is present per species
	std::vector<real> spec_mass_in_rate;  // solute boundary influx rate [mg/s]
	std::vector<real> spec_mass_out_rate; // solute boundary outflux rate [mg/s]

public:
	inline int copyHostIntVectorToView(const std::vector<int> &tmpcells, intArr &cells)
	{
		int ncells = int(tmpcells.size());
		if (ncells <= 0)
			return ncells;
		cells = intArr("cells", ncells);
#if defined(KOKKOS_ENABLE_CUDA)
		int err = cudaMemcpyAsync(cells.data(), tmpcells.data(), ncells * sizeof(int),
								  cudaMemcpyHostToDevice);
		err = cudaDeviceSynchronize();
#elif defined(KOKKOS_ENABLE_HIP)
		int err = hipMemcpyAsync(cells.data(), tmpcells.data(), ncells * sizeof(int),
								 hipMemcpyHostToDevice);
		err = hipDeviceSynchronize();
#elif defined(KOKKOS_ENABLE_SYCL)
		sycl::queue q{sycl::gpu_selector_v};
		q.memcpy(cells.data(), tmpcells.data(), ncells * sizeof(int));
		q.wait();
#else
		std::memcpy(cells.data(), tmpcells.data(), ncells * sizeof(int));
#endif
		return ncells;
	}

	inline int find_bcells(State &state, std::string &id, const Domain &dom, Parallel &par, int nPoly, realArr &xPoly, realArr &yPoly)
	{
#if SERGHEI_DEBUG_BOUNDARY
		std::cout << GGD << "find_bcells called for boundary id =" << RED << id << RESET << std::endl;
#endif
		Kokkos::Timer timermpi;
		int foundInSubdom;			 // to keep track of which subdomains are associated to this boundary
		std::vector<int> tmpbcells;	 // array of indexes of boundary cells
		std::vector<int> subdomains; // keeps track of which subdomains are associated to the BC

		for (int iGlob = 0; iGlob < dom.nCell; iGlob++)
		{
			int i, j;
			dom.unpackIndices(iGlob, j, i);
			int ii = dom.getHaloExtension(i, j);
			foundInSubdom = -1;
			if (!state.isnodata(ii))
			{
				if ((j == 0 && dom.iN) || (j == dom.ny - 1 && dom.iS) || (i == 0 && dom.iW) || (i == dom.nx - 1 && dom.iE) ||
					state.isnodata(ii + 1) || state.isnodata(ii - 1) ||
					state.isnodata(ii - (dom.nx + 2 * dom.hc)) || state.isnodata(ii + (dom.nx + 2 * dom.hc)))
				{
					// boundary domain || nodata neighbours
					real xCoord = dom.xll + (par.i_beg + i + 0.5) * dom.dxConst;
					real yCoord = dom.yll + dom.ny_glob * dom.dxConst - (par.j_beg + j + 0.5) * dom.dxConst;
					if (geometry::isInsidePoly(nPoly, xPoly, yPoly, xCoord, yCoord))
					{
						tmpbcells.push_back(ii);
					}
				}
			}
		}

		ncellsBC = int(tmpbcells.size());
		if (ncellsBC > 0)
			foundInSubdom = par.myrank; // if at least one cell in this subdomain (rank) is in the BC, tag as found

		int ncells_all;

		timermpi.reset();
		SERGHEI_MPI_CHECK(MPI_Allreduce(&ncellsBC, &ncells_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));
		if (par.nranks > 1)
			dom.timers.swe.bc.mpi += timermpi.seconds();

		int *subdoms;
		subdoms = (int *)malloc(par.nranks * sizeof(int));
		timermpi.reset();
		SERGHEI_MPI_CHECK(MPI_Allgather(&foundInSubdom, 1, MPI_INT, subdoms, 1, MPI_INT, MPI_COMM_WORLD));
		if (par.nranks > 1)
			dom.timers.swe.bc.mpi += timermpi.seconds();

#if SERGHEI_DEBUG_BOUNDARY
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << " rank " << RED << par.myrank << RESET << std::endl;
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << " ncellsBC " << RED << ncells_all << RESET << std::endl;
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << " BC subdomains: ";
		for (int i = 0; i < par.nranks; i++)
		{
			std::cout << " ";
			if (subdoms[i] == par.myrank)
				std::cout << RED;
			std::cout << subdoms[i] << "\t" << RESET;
		}
		std::cout << std::endl;
#endif

		for (int i = 0; i < par.nranks; i++)
		{
			if (subdoms[i] >= 0)
			{
				subdomains.push_back(subdoms[i]);
			}
		}
#if SERGHEI_DEBUG_BOUNDARY
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "Consolidated " << RED << subdomains.size() << RESET << " BC subdomains = ";
		for (int i = 0; i < subdomains.size(); i++)
		{
			std::cout << " ";
			std::cout << subdomains[i] << "\t";
		}
		std::cout << std::endl;
#endif
		timermpi.reset();
		MPI_Group group, subgroup;
		SERGHEI_MPI_CHECK(MPI_Comm_group(MPI_COMM_WORLD, &group));
		SERGHEI_MPI_CHECK(MPI_Group_incl(group, subdomains.size(), subdomains.data(), &subgroup));
		SERGHEI_MPI_CHECK(MPI_Comm_create(MPI_COMM_WORLD, subgroup, &comm));
		if (par.nranks > 1)
			dom.timers.swe.bc.mpi += timermpi.seconds();

		// int err;
		//  we need the total boundary cells detected by all subdomain to launch an error otherwise
		if (ncells_all > 0)
		{
			bcells = intArr("bcells", ncellsBC);
#if defined(KOKKOS_ENABLE_CUDA)
			int err = cudaMemcpyAsync(bcells.data(), tmpbcells.data(), ncellsBC * sizeof(int),
									  cudaMemcpyHostToDevice);
			err = cudaDeviceSynchronize();
#elif defined(KOKKOS_ENABLE_HIP)
			int err = hipMemcpyAsync(bcells.data(), tmpbcells.data(), ncellsBC * sizeof(int),
									 hipMemcpyHostToDevice);
			err = hipDeviceSynchronize();
#elif defined(KOKKOS_ENABLE_SYCL)
			sycl::queue q{sycl::gpu_selector_v};
			q.memcpy(bcells.data(), tmpbcells.data(), ncellsBC * sizeof(int));
			q.wait();
#else
			std::memcpy(bcells.data(), tmpbcells.data(), ncellsBC * sizeof(int));
#endif
		}
		else
		{
			if (par.masterproc)
				std::cerr << RERROR << "No boundary cells found for external boundary with id '" << id << "'" << std::endl;
			return 0;
		}
		return 1;
	}

	inline int find_trigger_cells(State &state, std::string &id, const Domain &dom, Parallel &par, int nPoly, realArr &xPoly, realArr &yPoly)
	{
#if SERGHEI_DEBUG_BOUNDARY
		std::cout << GGD << "find_trigger_cells called for boundary id =" << RED << id << RESET << std::endl;
#endif
		std::vector<int> tmpcells;

		if (nPoly == 0)
		{
			for (int iGlob = 0; iGlob < dom.nCell; iGlob++)
			{
				int i, j;
				dom.unpackIndices(iGlob, j, i);
				int ii = dom.getHaloExtension(i, j);
				if (!state.isnodata(ii))
					tmpcells.push_back(ii);
			}
		}
		else if (nPoly == 1)
		{
			real xObs = xPoly(0);
			real yObs = yPoly(0);
			real localMinDist2 = 1.0e30;
			int localCell = -1;
			for (int iGlob = 0; iGlob < dom.nCell; iGlob++)
			{
				int i, j;
				dom.unpackIndices(iGlob, j, i);
				int ii = dom.getHaloExtension(i, j);
				if (!state.isnodata(ii))
				{
					real xCoord = dom.xll + (par.i_beg + i + 0.5) * dom.dxConst;
					real yCoord = dom.yll + dom.ny_glob * dom.dxConst - (par.j_beg + j + 0.5) * dom.dxConst;
					real dx = xCoord - xObs;
					real dy = yCoord - yObs;
					real dist2 = dx * dx + dy * dy;
					if (dist2 < localMinDist2)
					{
						localMinDist2 = dist2;
						localCell = ii;
					}
				}
			}
			real globalMinDist2 = 1.0e30;
			SERGHEI_MPI_CHECK(MPI_Allreduce(&localMinDist2, &globalMinDist2, 1, SERGHEI_MPI_REAL, MPI_MIN, MPI_COMM_WORLD));
			if (localCell >= 0 && localMinDist2 <= globalMinDist2 + TOL12)
				tmpcells.push_back(localCell);
		}
		else if (nPoly >= 3)
		{
			for (int iGlob = 0; iGlob < dom.nCell; iGlob++)
			{
				int i, j;
				dom.unpackIndices(iGlob, j, i);
				int ii = dom.getHaloExtension(i, j);
				if (!state.isnodata(ii))
				{
					real xCoord = dom.xll + (par.i_beg + i + 0.5) * dom.dxConst;
					real yCoord = dom.yll + dom.ny_glob * dom.dxConst - (par.j_beg + j + 0.5) * dom.dxConst;
					if (geometry::isInsidePoly(nPoly, xPoly, yPoly, xCoord, yCoord))
						tmpcells.push_back(ii);
				}
			}
		}
		else
		{
			if (par.masterproc)
				std::cerr << RERROR << "Invalid trigger polygon for boundary id '" << id << "'" << std::endl;
			return 0;
		}

		ncellsTrigger = copyHostIntVectorToView(tmpcells, triggerCells);
		int ncells_all = 0;
		SERGHEI_MPI_CHECK(MPI_Allreduce(&ncellsTrigger, &ncells_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));
		if (ncells_all <= 0 && par.masterproc)
		{
			std::cerr << YEXC << "No trigger observation cells found for boundary id '" << id << "'. Irrigation trigger will stay off." << std::endl;
		}
		return 1;
	}

	void inline getMinBedElevation(Domain const &dom, State &state)
	{
		Kokkos::Timer timermpi;
		Kokkos::parallel_reduce("swe_bc_z_min", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &zMin) {
			int ii = bcells[iGlob];
				real z = state.z(ii);
				zMin = min(zMin,z); }, Kokkos::Min<real>(zMin));
		real zMin_all;
		timermpi.reset();
		SERGHEI_MPI_CHECK_COMM(MPI_Allreduce(&zMin, &zMin_all, 1, SERGHEI_MPI_REAL, MPI_MIN, comm), comm);
		if (dom.nsubdom > 1)
			dom.timers.swe.bc.mpi += timermpi.seconds();
		zMin = zMin_all;

		real zMax = -1E6;

		Kokkos::parallel_reduce("swe_bc_z_max", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &zMax) {
		  int ii = bcells[iGlob];
			real z = state.z(ii);
			zMax = max(zMax,z); }, Kokkos::Max<real>(zMax));
		real zMax_all;
		timermpi.reset();
		SERGHEI_MPI_CHECK_COMM(MPI_Allreduce(&zMax, &zMax_all, 1, SERGHEI_MPI_REAL, MPI_MAX, comm), comm);
		if (dom.nsubdom > 1)
			dom.timers.swe.bc.mpi += timermpi.seconds();
		zMax = zMax_all;

		nzMin = 0;
		Kokkos::parallel_reduce("swe_bc_z_min_count", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, int &nzMin) {
		  int ii = bcells[iGlob];
			if(state.z(ii) == zMin) nzMin++; }, Kokkos::Sum<int>(nzMin));
		int nzMin_all;
		timermpi.reset();
		SERGHEI_MPI_CHECK_COMM(MPI_Allreduce(&nzMin, &nzMin_all, 1, MPI_INT, MPI_SUM, comm), comm);
		if (dom.nsubdom > 1)
			dom.timers.swe.bc.mpi += timermpi.seconds();
		nzMin = nzMin_all;

#if SERGHEI_DEBUG_BOUNDARY
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "zMin = " << zMin << "\tnzMin = " << nzMin << std::endl;
#endif
	}

	void inline flattenWaterSurface(State &state, Domain &dom)
	{
#if SERGHEI_DEBUG_BOUNDARY
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
#endif
		Kokkos::Timer timermpi;
		hzMin = zMin = 1E6;
		real extraVol = 0;
		real extraArea = 0;

		// find the lowest water surface in the cross section
		Kokkos::parallel_reduce("swe_bc_flatten_hz", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &hzMin) {
		  int ii = bcells[iGlob];
			real h = state.h(ii);
			real z = state.z(ii);
			if(h > state.hmin) hzMin = min(hzMin,h+z); }, Kokkos::Min<real>(hzMin));

		real hzMin_all;
		timermpi.reset();
		SERGHEI_MPI_CHECK_COMM(MPI_Allreduce(&hzMin, &hzMin_all, 1, SERGHEI_MPI_REAL, MPI_MIN, comm), comm);
		if (dom.nsubdom > 1)
			dom.timers.swe.bc.mpi += timermpi.seconds();
		hzMin = hzMin_all;

		// find the volume above the lowest water surface
		real area = dom.cellArea(); // WARNING assumes uniform mesh
		Kokkos::parallel_reduce("swe_bc_vol_minhz", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &volume) {
			int ii = bcells[iGlob];
			real h = state.h(ii);
			real z = state.z(ii);
			if(h > state.hmin){
				if(z <= hzMin){
					volume += area * (h + z - hzMin);
				}else{
					volume += area * h;
				}
			} }, Kokkos::Sum<real>(extraVol));
		real extraVol_all;

		timermpi.reset();
		SERGHEI_MPI_CHECK_COMM(MPI_Allreduce(&extraVol, &extraVol_all, 1, SERGHEI_MPI_REAL, MPI_SUM, comm), comm);
		if (dom.nsubdom > 1)
			dom.timers.swe.bc.mpi += timermpi.seconds();

		extraVol = extraVol_all;

		if (extraVol > ZERO)
		{
			// compute the surface area with water level higher than the cross-sectional minimum
			Kokkos::parallel_reduce("swe_bc_area_minhz", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &sumarea) {
				int ii = bcells[iGlob];
				real h = state.h(ii);
				real z = state.z(ii);
				if(h > state.hmin && z < hzMin && hzMin-z >0){
					sumarea += area;
				} }, Kokkos::Sum<real>(extraArea));
			real extraArea_all;
			timermpi.reset();
			SERGHEI_MPI_CHECK_COMM(MPI_Allreduce(&extraArea, &extraArea_all, 1, SERGHEI_MPI_REAL, MPI_SUM, comm), comm);
			if (dom.nsubdom > 1)
				dom.timers.swe.bc.mpi += timermpi.seconds();
			extraArea = extraArea_all;

			// compute and assign homogenised water surface
			real hz = hzMin + extraVol / extraArea;
			Kokkos::parallel_for("swe_bc_h_minhz", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob) {
				int ii = bcells[iGlob];
				real h = state.h(ii);
				real z = state.z(ii);
				if(h > state.hmin){
					state.h(ii) = 0.;
					if(z <= hzMin) state.h(ii) = max(hz - z, (real) 0.);
				} });
		}
	}

	void inline distributeDischarge(State &state, Domain const &dom, real Q)
	{
#if SERGHEI_DEBUG_BOUNDARY
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
#endif
		Kokkos::Timer timermpi;
		real hsum;
		// 1. Total water depth over the boundary, counting wet cells only
		Kokkos::parallel_reduce("swe_bc_h_weight", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &hsum) {
			int ii = bcells[iGlob];
			if (state.h(ii) > state.hmin) {
				hsum += state.h(ii);
			} }, Kokkos::Sum<real>(hsum));

		real hsum_all;
		timermpi.reset();
		SERGHEI_MPI_CHECK_COMM(MPI_Allreduce(&hsum, &hsum_all, 1, SERGHEI_MPI_REAL, MPI_SUM, comm), comm);
		if (dom.nsubdom > 1)
			dom.timers.swe.bc.mpi += timermpi.seconds();
		hsum = hsum_all;

		// 2. Total number of boundary cells, used to derive the unit discharge
		int ncells_all;
		int local_ncells = ncellsBC;
		timermpi.reset();
		SERGHEI_MPI_CHECK_COMM(MPI_Allreduce(&local_ncells, &ncells_all, 1, MPI_INT, MPI_SUM, comm), comm);
		if (dom.nsubdom > 1)
			dom.timers.swe.bc.mpi += timermpi.seconds();

		real ds = dom.dx();
		real total_width = ncells_all * ds;
		real q_unit = 0.0;
		if (total_width > 0)
		{
			q_unit = Q / total_width; // unit discharge
		}

		Kokkos::parallel_for("swe_bc_distribute", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob) {
			int ii = bcells[iGlob];
			real h = state.h(ii);
			real z = state.z(ii);
			real weight = 0;

			// Case A: boundary already wet (distribute the discharge proportionally to depth for a smooth water surface)
			if (hsum > 1e-6 && h > state.hmin)
			{
				weight = h / hsum;
				state.hu(ii) = Q * weight / ds * normalx;
				state.hv(ii) = Q * weight / ds * normaly;
			}
			// Case B: dry boundary -> raise the water depth to the critical depth
			else
			{
				// Critical depth hc = (q^2 / g)^(1/3)
				real hc_val = cbrt((q_unit * q_unit) / GRAV);
				if (hc_val < state.hmin)
					hc_val = state.hmin; // keep it above the minimum depth

				// Force the depth to the critical value so momentum is not zeroed
				state.h(ii) = hc_val;
				h = hc_val; // update the local copy used for DEBUG printing

				real hu_val = fabs(q_unit);
				state.hu(ii) = hu_val * normalx;
				state.hv(ii) = hu_val * normaly;
			}

#if SERGHEI_DEBUG_BOUNDARY > 1
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "Q = " << Q << "\tii = " << ii << "\tz = " << z << "\t h = " << h << "\tweight = " << weight << "\tds = " << ds << "\t(hu,hv) = " << state.hu(ii) << " " << state.hv(ii) << std::endl;
#endif
		});
	}

	void inline distributeDrainageDischarge(State &state, Domain const &dom, real Q)
	{
#if SERGHEI_DEBUG_BOUNDARY
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
#endif
		Kokkos::Timer timermpi;
		real hsum = 0.0;
		int wetCount = 0;
		Kokkos::parallel_reduce("swe_bc_drainage_h_weight", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &sumH, int &countWet) {
			int ii = bcells[iGlob];
			if (state.h(ii) > state.hmin)
			{
				sumH += state.h(ii);
				countWet++;
			}
		}, Kokkos::Sum<real>(hsum), Kokkos::Sum<int>(wetCount));

		real hsum_all = 0.0;
		int wetCount_all = 0;
		timermpi.reset();
		SERGHEI_MPI_CHECK_COMM(MPI_Allreduce(&hsum, &hsum_all, 1, SERGHEI_MPI_REAL, MPI_SUM, comm), comm);
		SERGHEI_MPI_CHECK_COMM(MPI_Allreduce(&wetCount, &wetCount_all, 1, MPI_INT, MPI_SUM, comm), comm);
		if (dom.nsubdom > 1)
			dom.timers.swe.bc.mpi += timermpi.seconds();

		hsum = hsum_all;
		wetCount = wetCount_all;
		real ds = dom.dx();
		real q_unit = 0.0;
		if (wetCount > 0)
			q_unit = Q / (wetCount * ds);

		Kokkos::parallel_for("swe_bc_drainage_distribute", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob) {
			int ii = bcells[iGlob];
			real h = state.h(ii);
			if (h > state.hmin && hsum > 1e-6)
			{
				real weight = h / hsum;
				state.hu(ii) = Q * weight / ds * normalx;
				state.hv(ii) = Q * weight / ds * normaly;
			}
			else if (h > state.hmin && wetCount > 0)
			{
				state.hu(ii) = q_unit * normalx;
				state.hv(ii) = q_unit * normaly;
			}
			else
			{
				state.hu(ii) = 0.0;
				state.hv(ii) = 0.0;
			}
		});
	}

		inline int findIrrigationStageByDAT(real dat)
		{
			for (int i = 0; i < nIrrigationStages; i++)
			{
				if (dat >= irrigationDatMin(i) && dat <= irrigationDatMax(i))
					return i;
			}
			return -1;
		}

		inline int findIrrigationStageByDVS(real dvs)
		{
			for (int i = 0; i < nIrrigationStages; i++)
			{
				if (dvs >= irrigationDvsMin(i) && dvs <= irrigationDvsMax(i))
					return i;
			}
			return -1;
		}

		inline int selectIrrigationStage(Domain &dom)
		{
			if (nIrrigationStages <= 0)
			{
				if (!irrigationStrategyWarned)
				{
					std::cerr << YEXC << "Triggered irrigation boundary has no irrigation strategy table; irrigation is disabled." << std::endl;
					irrigationStrategyWarned = 1;
				}
				return -1;
			}

			if (useDVSForIrrigation == 1)
			{
				if (!irrigationStrategyWarned)
				{
					std::cerr << YEXC << "Triggered irrigation use_DVS=1 was requested, but DVS data are not connected to ExtBC yet. Irrigation is disabled." << std::endl;
					irrigationStrategyWarned = 1;
				}
				return -1;
			}

			irrigationDAT = floor((dom.etime - cropStartTime) / 86400.0) + 1.0;
			int stage = findIrrigationStageByDAT(irrigationDAT);
			if (stage < 0 && !irrigationStrategyWarned)
			{
				std::cerr << YEXC << "Triggered irrigation has no strategy row for DAT=" << irrigationDAT << "; irrigation is disabled." << std::endl;
				irrigationStrategyWarned = 1;
			}
			return stage;
		}

#if CROP_GROWTH_MODEL
		inline int selectIrrigationStage(Domain &dom, Wofost72 &wofost, State &state)
		{
			if (useDVSForIrrigation == 0)
				return selectIrrigationStage(dom);

			if (nIrrigationStages <= 0)
			{
				if (!irrigationStrategyWarned)
				{
					std::cerr << YEXC << "Triggered irrigation boundary has no irrigation strategy table; irrigation is disabled." << std::endl;
					irrigationStrategyWarned = 1;
				}
				return -1;
			}

			real localDvsSum = 0.0;
			int localCount = 0;
			auto cropDVS = wofost.s.DVS;
			if (ncellsTrigger > 0)
			{
				Kokkos::parallel_reduce("swe_bc_irrigation_trigger_dvs", ncellsTrigger, KOKKOS_CLASS_LAMBDA(int iTrig, real &dvsSum, int &count) {
					int ii = triggerCells(iTrig);
					int iLocal = ii % (dom.nx + 2 * dom.hc) - dom.hc;
					int jLocal = ii / (dom.nx + 2 * dom.hc) - dom.hc;
					if (iLocal >= 0 && iLocal < dom.nx && jLocal >= 0 && jLocal < dom.ny && !state.isnodata(ii))
					{
						int iDom = jLocal * dom.nx + iLocal;
						dvsSum += cropDVS(iDom);
						count++;
					}
				}, Kokkos::Sum<real>(localDvsSum), Kokkos::Sum<int>(localCount));
			}

			real dvsSumAll = 0.0;
			int countAll = 0;
			Kokkos::Timer timermpi;
			SERGHEI_MPI_CHECK(MPI_Allreduce(&localDvsSum, &dvsSumAll, 1, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD));
			SERGHEI_MPI_CHECK(MPI_Allreduce(&localCount, &countAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));
			if (dom.nsubdom > 1)
				dom.timers.swe.bc.mpi += timermpi.seconds();

			if (countAll <= 0)
			{
				if (!irrigationStrategyWarned)
				{
					std::cerr << YEXC << "Triggered irrigation use_DVS=1 has no valid crop DVS cells in trigger_polygon; irrigation is disabled." << std::endl;
					irrigationStrategyWarned = 1;
				}
				return -1;
			}

			irrigationMeanDVS = dvsSumAll / countAll;
			int stage = findIrrigationStageByDVS(irrigationMeanDVS);
			if (stage < 0 && !irrigationStrategyWarned)
			{
				std::cerr << YEXC << "Triggered irrigation has no strategy row for DVS=" << irrigationMeanDVS << "; irrigation is disabled." << std::endl;
				irrigationStrategyWarned = 1;
			}
			return stage;
		}
#endif

		inline int checkIrrigationTrigger(State &state, Domain &dom, const int stage, const realArr &RD = realArr())
		{
			real localHsum = 0.0;
			int localCount = 0;

		if (ncellsTrigger > 0)
		{
			Kokkos::parallel_reduce("swe_bc_irrigation_trigger_h", ncellsTrigger, KOKKOS_CLASS_LAMBDA(int iGlob, real &hsum, int &count) {
				int ii = triggerCells[iGlob];
				if (!state.isnodata(ii))
				{
					hsum += state.h(ii);
					count++;
				} }, Kokkos::Sum<real>(localHsum), Kokkos::Sum<int>(localCount));
		}

		real hsumAll = 0.0;
		int countAll = 0;
		Kokkos::Timer timermpi;
		SERGHEI_MPI_CHECK(MPI_Allreduce(&localHsum, &hsumAll, 1, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD));
		SERGHEI_MPI_CHECK(MPI_Allreduce(&localCount, &countAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));
		if (dom.nsubdom > 1)
			dom.timers.swe.bc.mpi += timermpi.seconds();

		if (countAll <= 0)
		{
			irrigationActive = 0;
			irrigationMeanH = 0.0;
			return 0;
		}

			irrigationMeanH = hsumAll / countAll;

			real hTrigger = irrigationHMin(stage);
			real hStop = irrigationHMax(stage);
			if (hStop < hTrigger)
				hStop = hTrigger;

			int mode = 1;
			if (bcvals.extent(0) > 0)
				mode = int(bcvals(0) + 0.5);
			if (irrigationTrigType.extent(0) > stage && irrigationTrigType(stage) != 0)
				mode = irrigationTrigType(stage);
			real thetaTrig = -1.0;
			real thetaStop = -1.0;
			if (irrigationThetaTrig.extent(0) > stage)
				thetaTrig = irrigationThetaTrig(stage);
			if (irrigationThetaStop.extent(0) > stage)
				thetaStop = irrigationThetaStop(stage);

			// 根系层土壤水触发值：触发单元上的平均相对饱和度
			// S = Σ(wc·dz)/Σ(θs·dz)，与 "θ_rz < frac·θsat,rz" 严格等价。
			// 每个水动力步用当前 gw.wc 新鲜计算（wofost.ets.RZSM 每天只刷新一次，太陈旧；
			// 异步耦合下 θ 天然滞后一个 GW 子循环间隔）。
#if CROP_GROWTH_MODEL && SERGHEI_RE_MODEL
			irrigationMeanTheta = -1.0;
			if (mode >= 2 && mode <= 4 && gw != nullptr && gdom != nullptr && RD.extent(0) > 0)
			{
				real localWcDz = 0.0;
				real localThsDz = 0.0;
				int localThetaCount = 0;
				{
					// view 本地拷贝后按值捕获（勿按引用捕获宿主对象，CUDA 后端不安全）
					auto wc = gw->wc;
					auto vgTab = gw->vgTable;
					auto sID = gw->soilID;
					auto dzArr = gdom->dz;
					auto depthArr = gdom->depth;
					const int ghc = gdom->hc;
					const int gnxhc = gdom->nxhc;
					const int gnyhc = gdom->nyhc;
					const int gnz = gdom->nz;
					Kokkos::parallel_reduce("swe_bc_irrigation_trigger_theta", ncellsTrigger, KOKKOS_CLASS_LAMBDA(int iTrig, real &sumWc, real &sumThs, int &cnt) {
						int ii = triggerCells(iTrig);
						int iLocal = ii % (dom.nx + 2 * dom.hc) - dom.hc;
						int jLocal = ii / (dom.nx + 2 * dom.hc) - dom.hc;
						if (iLocal >= 0 && iLocal < dom.nx && jLocal >= 0 && jLocal < dom.ny && !state.isnodata(ii))
						{
							int iDom = jLocal * dom.nx + iLocal;
							real sWc = 0.0;
							real sThs = 0.0;
							SourceSinkCrop::rootZoneMoistureSums(wc, vgTab, sID, dzArr, depthArr,
																ghc, gnxhc, gnyhc, gnz, iLocal, jLocal,
																RD(iDom) * 0.01, sWc, sThs);
							if (sThs > 0.0)
							{
								sumWc += sWc;
								sumThs += sThs;
								cnt++;
							}
						} }, Kokkos::Sum<real>(localWcDz), Kokkos::Sum<real>(localThsDz), Kokkos::Sum<int>(localThetaCount));
				}

				real wcDzAll = 0.0;
				real thsDzAll = 0.0;
				int thetaCountAll = 0;
				timermpi.reset();
				SERGHEI_MPI_CHECK(MPI_Allreduce(&localWcDz, &wcDzAll, 1, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD));
				SERGHEI_MPI_CHECK(MPI_Allreduce(&localThsDz, &thsDzAll, 1, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD));
				SERGHEI_MPI_CHECK(MPI_Allreduce(&localThetaCount, &thetaCountAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));
				if (dom.nsubdom > 1)
					dom.timers.swe.bc.mpi += timermpi.seconds();

				if (thetaCountAll > 0 && thsDzAll > TOL12)
					irrigationMeanTheta = wcDzAll / thsDzAll;
			}

			if (mode >= 2 && mode <= 4 && irrigationMeanTheta < 0.0)
			{
				if (!irrigationThetaWarned)
				{
					std::cerr << YEXC << "Triggered irrigation soil-moisture trigger requested, but root-zone saturation is not computable (no crop/GW data on the trigger cells); irrigation stays off." << std::endl;
					irrigationThetaWarned = 1;
				}
				irrigationActive = 0;
				return 0;
			}
#endif

			// 滞回真值表：
			// 开启（inactive）：mode1 meanH<hTrigger；mode2 S<θtrig；mode3 任一；mode4 同时。
			// 继续（active，四种模式统一）：meanH<hStop 且 θ 未达提前停阈值。
			// hStop 恒为硬上限（防"田面已满但土壤干"时无界灌水/振荡补灌）；
			// thetaStop<=0 表示提前停条件透明（勿做字面比较，否则恒真导致立即关闭）。
			int swOpen = (irrigationMeanH < hTrigger) ? 1 : 0;
#if CROP_GROWTH_MODEL && SERGHEI_RE_MODEL
			int thetaOpen = (irrigationMeanTheta < thetaTrig) ? 1 : 0;
#else
			int thetaOpen = 0;
			(void)thetaTrig;
#endif
			int openCond = 0;
			if (mode == 1)
				openCond = swOpen;
			else if (mode == 2)
				openCond = thetaOpen;
			else if (mode == 3)
				openCond = (swOpen || thetaOpen) ? 1 : 0;
			else if (mode == 4)
				openCond = (swOpen && thetaOpen) ? 1 : 0;
			else
			{
				std::cerr << YEXC << "Invalid triggered irrigation mode " << mode << "; using surface-water trigger only." << std::endl;
				openCond = swOpen;
			}

			int contCond = (irrigationMeanH < hStop) ? 1 : 0;
#if CROP_GROWTH_MODEL && SERGHEI_RE_MODEL
			if (contCond && thetaStop > 0.0 && irrigationMeanTheta >= thetaStop)
				contCond = 0;
#else
			(void)thetaStop;
#endif

			int active = ((irrigationActive ? contCond : openCond) && irrigationMeanH < hStop) ? 1 : 0;

#if SERGHEI_DEBUG_BOUNDARY > 1
			std::cout << GGD << "checkIrrigationTrigger: stage=" << stage << " mode=" << mode
					  << " meanH=" << irrigationMeanH << " S=" << irrigationMeanTheta << " active=" << active << std::endl;
#endif

		irrigationActive = active;
		return active;
	}

	inline void closeTriggeredIrrigation(State &state)
	{
		if (ncellsBC <= 0)
			return;
		Kokkos::parallel_for("swe_bc_irrigation_off", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob) {
			int ii = bcells[iGlob];
			state.hu(ii) = 0.0;
			state.hv(ii) = 0.0;
		});
	}

		inline real applyTriggeredIrrigationStage(State &state, Domain &dom, const int stage, const realArr &RD = realArr())
		{
			real extraMass = 0.0;
			irrigationAppliedDischarge = 0.0;
			currentIrrigationStage = stage;
			if (stage < 0 || irrigationEnabled(stage) == 0)
			{
				irrigationActive = 0;
				closeTriggeredIrrigation(state);
				return 0.0;
			}

			int active = checkIrrigationTrigger(state, dom, stage, RD);
			if (!active)
			{
				closeTriggeredIrrigation(state);
				return 0.0;
			}

			real Qirr = irrigationQ(stage);
			irrigationAppliedDischarge = Qirr;
			if (ncellsBC <= 0)
				return 0.0;

			flattenWaterSurface(state, dom);
			distributeDischarge(state, dom, Qirr);

			return extraMass;
		}

		inline real applyTriggeredIrrigation(State &state, Domain &dom)
		{
#if SERGHEI_DEBUG_BOUNDARY
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "SWE_BC_TRIGGERED_IRRIGATION" << std::endl;
#endif
			int stage = selectIrrigationStage(dom);
			return applyTriggeredIrrigationStage(state, dom, stage);
		}

#if CROP_GROWTH_MODEL
		inline real applyTriggeredIrrigation(State &state, Domain &dom, Wofost72 &wofost)
		{
#if SERGHEI_DEBUG_BOUNDARY
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "SWE_BC_TRIGGERED_IRRIGATION" << std::endl;
#endif
			int stage = selectIrrigationStage(dom, wofost, state);
			return applyTriggeredIrrigationStage(state, dom, stage, wofost.s.RD);
		}
#endif

		inline int findDrainageStageByDAT(real dat)
		{
			for (int i = 0; i < nDrainageStages; i++)
			{
				if (dat >= drainageDatMin(i) && dat <= drainageDatMax(i))
					return i;
			}
			return -1;
		}

		inline int findDrainageStageByDVS(real dvs)
		{
			for (int i = 0; i < nDrainageStages; i++)
			{
				if (dvs >= drainageDvsMin(i) && dvs <= drainageDvsMax(i))
					return i;
			}
			return -1;
		}

		inline int selectDrainageStage(Domain &dom)
		{
			if (nDrainageStages <= 0)
			{
				if (!drainageStrategyWarned)
				{
					std::cerr << YEXC << "Triggered drainage boundary has no drainage strategy table; drainage is disabled." << std::endl;
					drainageStrategyWarned = 1;
				}
				return -1;
			}

			if (useDVSForDrainage == 1)
			{
				if (!drainageStrategyWarned)
				{
					std::cerr << YEXC << "Triggered drainage use_DVS=1 was requested, but DVS data are not connected to this ExtBC call. Drainage is disabled." << std::endl;
					drainageStrategyWarned = 1;
				}
				return -1;
			}

			drainageDAT = floor((dom.etime - cropStartTime) / 86400.0) + 1.0;
			int stage = findDrainageStageByDAT(drainageDAT);
			if (stage < 0 && !drainageStrategyWarned)
			{
				std::cerr << YEXC << "Triggered drainage has no strategy row for DAT=" << drainageDAT << "; drainage is disabled." << std::endl;
				drainageStrategyWarned = 1;
			}
			return stage;
		}

#if CROP_GROWTH_MODEL
		inline int selectDrainageStage(Domain &dom, Wofost72 &wofost, State &state)
		{
			if (useDVSForDrainage == 0)
				return selectDrainageStage(dom);

			if (nDrainageStages <= 0)
			{
				if (!drainageStrategyWarned)
				{
					std::cerr << YEXC << "Triggered drainage boundary has no drainage strategy table; drainage is disabled." << std::endl;
					drainageStrategyWarned = 1;
				}
				return -1;
			}

			real localDvsSum = 0.0;
			int localCount = 0;
			auto cropDVS = wofost.s.DVS;
			if (ncellsTrigger > 0)
			{
				Kokkos::parallel_reduce("swe_bc_drainage_trigger_dvs", ncellsTrigger, KOKKOS_CLASS_LAMBDA(int iTrig, real &dvsSum, int &count) {
					int ii = triggerCells(iTrig);
					int iLocal = ii % (dom.nx + 2 * dom.hc) - dom.hc;
					int jLocal = ii / (dom.nx + 2 * dom.hc) - dom.hc;
					if (iLocal >= 0 && iLocal < dom.nx && jLocal >= 0 && jLocal < dom.ny && !state.isnodata(ii))
					{
						int iDom = jLocal * dom.nx + iLocal;
						dvsSum += cropDVS(iDom);
						count++;
					}
				}, Kokkos::Sum<real>(localDvsSum), Kokkos::Sum<int>(localCount));
			}

			real dvsSumAll = 0.0;
			int countAll = 0;
			Kokkos::Timer timermpi;
			SERGHEI_MPI_CHECK(MPI_Allreduce(&localDvsSum, &dvsSumAll, 1, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD));
			SERGHEI_MPI_CHECK(MPI_Allreduce(&localCount, &countAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));
			if (dom.nsubdom > 1)
				dom.timers.swe.bc.mpi += timermpi.seconds();

			if (countAll <= 0)
			{
				if (!drainageStrategyWarned)
				{
					std::cerr << YEXC << "Triggered drainage use_DVS=1 has no valid crop DVS cells in trigger_polygon; drainage is disabled." << std::endl;
					drainageStrategyWarned = 1;
				}
				return -1;
			}

			drainageMeanDVS = dvsSumAll / countAll;
			int stage = findDrainageStageByDVS(drainageMeanDVS);
			if (stage < 0 && !drainageStrategyWarned)
			{
				std::cerr << YEXC << "Triggered drainage has no strategy row for DVS=" << drainageMeanDVS << "; drainage is disabled." << std::endl;
				drainageStrategyWarned = 1;
			}
			return stage;
		}
#endif

		inline int checkDrainageTrigger(State &state, Domain &dom, const int stage)
		{
			real localHsum = 0.0;
			int localCount = 0;

			if (ncellsTrigger > 0)
			{
				Kokkos::parallel_reduce("swe_bc_drainage_trigger_h", ncellsTrigger, KOKKOS_CLASS_LAMBDA(int iGlob, real &hsum, int &count) {
					int ii = triggerCells[iGlob];
					if (!state.isnodata(ii))
					{
						hsum += state.h(ii);
						count++;
					}
				}, Kokkos::Sum<real>(localHsum), Kokkos::Sum<int>(localCount));
			}

			real hsumAll = 0.0;
			int countAll = 0;
			Kokkos::Timer timermpi;
			SERGHEI_MPI_CHECK(MPI_Allreduce(&localHsum, &hsumAll, 1, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD));
			SERGHEI_MPI_CHECK(MPI_Allreduce(&localCount, &countAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));
			if (dom.nsubdom > 1)
				dom.timers.swe.bc.mpi += timermpi.seconds();

			if (countAll <= 0)
			{
				drainageActive = 0;
				drainageMeanH = 0.0;
				return 0;
			}

			drainageMeanH = hsumAll / countAll;
			real hOn = drainageOnH(stage);
			real hOff = drainageOffH(stage);
			if (hOn < hOff)
				hOn = hOff;

			int active = 0;
			if (drainageActive)
				active = (drainageMeanH > hOff) ? 1 : 0;
			else
				active = (drainageMeanH > hOn) ? 1 : 0;

			drainageActive = active;
			return active;
		}

		inline void closeTriggeredDrainage(State &state)
		{
			if (ncellsBC <= 0)
				return;
			Kokkos::parallel_for("swe_bc_drainage_off", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob) {
				int ii = bcells[iGlob];
				state.hu(ii) = 0.0;
				state.hv(ii) = 0.0;
			});
		}

		inline real applyTriggeredDrainageStage(State &state, Domain &dom, const int stage)
		{
			drainageAppliedDischarge = 0.0;
			currentDrainageStage = stage;
			if (stage < 0 || drainageStageEnabled(stage) == 0)
			{
				drainageActive = 0;
				closeTriggeredDrainage(state);
				return 0.0;
			}

			int active = checkDrainageTrigger(state, dom, stage);
			if (!active)
			{
				closeTriggeredDrainage(state);
				return 0.0;
			}

			int mode = drainageMode(stage);
			if (mode == 1)
			{
				Kokkos::parallel_for("swe_bc_triggered_drainage_critical", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob) {
					int ii = bcells[iGlob];
					real h = state.h(ii);
					if (h >= state.hmin)
					{
						real vel = sqrt(GRAV * h);
						state.hu(ii) = vel * h * normalx;
						state.hv(ii) = vel * h * normaly;
					}
					else
					{
						state.hu(ii) = 0.0;
						state.hv(ii) = 0.0;
					}
				});
			}
			else if (mode == 2)
			{
				real Qdrain = drainageQ(stage);
				drainageAppliedDischarge = Qdrain;
				if (ncellsBC <= 0)
					return 0.0;
				distributeDrainageDischarge(state, dom, Qdrain);
			}
			else
			{
				if (!drainageStrategyWarned)
				{
					std::cerr << YEXC << "Invalid triggered drainage mode " << mode << "; drainage is disabled." << std::endl;
					drainageStrategyWarned = 1;
				}
				drainageActive = 0;
				closeTriggeredDrainage(state);
			}
			return 0.0;
		}

		inline real applyTriggeredDrainage(State &state, Domain &dom)
		{
#if SERGHEI_DEBUG_BOUNDARY
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "SWE_BC_TRIGGERED_DRAINAGE" << std::endl;
#endif
			int stage = selectDrainageStage(dom);
			return applyTriggeredDrainageStage(state, dom, stage);
		}

#if CROP_GROWTH_MODEL
		inline real applyTriggeredDrainage(State &state, Domain &dom, Wofost72 &wofost)
		{
#if SERGHEI_DEBUG_BOUNDARY
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "SWE_BC_TRIGGERED_DRAINAGE" << std::endl;
#endif
			int stage = selectDrainageStage(dom, wofost, state);
			return applyTriggeredDrainageStage(state, dom, stage);
		}
#endif

#if SERGHEI_SEDIMENT_TRANSPORT
	void inline distributeBedChange(State &state, Domain const &dom)
	{
#if SERGHEI_DEBUG_BOUNDARY
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
#endif
		real DZsum;
		real Asum;
		Kokkos::parallel_reduce("sed_bc_sumVol", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &var1, real &var2) {
			int ii = bcells[iGlob];
			real h = state.h(ii);

			state.z(ii) += state.sediment.bedExchangeVol(ii);

			var1 += state.sediment.bedExchangeVol(ii);
			if(h>TOL12){
				var2 += dom.cellArea();
			} }, Kokkos::Sum<real>(DZsum), Kokkos::Sum<real>(Asum));

		real DZsum_all, Asum_all;
		SERGHEI_MPI_CHECK_COMM(MPI_Allreduce(&DZsum, &DZsum_all, 1, SERGHEI_MPI_REAL, MPI_SUM, comm), comm);
		SERGHEI_MPI_CHECK_COMM(MPI_Allreduce(&Asum, &Asum_all, 1, SERGHEI_MPI_REAL, MPI_SUM, comm), comm);
		DZsum = DZsum_all;
		Asum = Asum_all;

		if (Asum > TOL12)
		{
			real excVol = DZsum / Asum;

			Kokkos::parallel_for("swe_bc_shareVol", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob) {
				int ii = bcells[iGlob];
				real h = state.h(ii);
				if (h > TOL12)
				{
					state.z(ii) -= excVol * dom.cellArea();
				}
#if SERGHEI_DEBUG_BOUNDARY > 1
				std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "Q = " << Q << "\tii = " << ii << "\tz = " << z << "\t h = " << h << "\tweight = " << weight << "\tds = " << ds << "\t(hu,hv) = " << state.hu(ii) << " " << state.hv(ii) << std::endl;
#endif
			});
		}
	}
#endif

		inline void apply(State &state, Domain &dom)
		{
#if SERGHEI_DEBUG_WORKFLOW
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
#endif
		Kokkos::fence("bc.apply-end");
		Kokkos::Timer timer;

		real extraMass = 0.0;

#if SERGHEI_SEDIMENT_TRANSPORT
		if (ncellsBC > 0)
		{
			switch (bctype)
			{
			// outlet boundary conditions
			case SWE_BC_CRITICAL:
			case SWE_BC_H_CONST:
			case SWE_BC_WSE_CONST:
			case SWE_BC_FREE_OUTFLOW:
			case SWE_BC_HZ_T_OUTLET:
			case SWE_BC_TRIGGERED_DRAINAGE:

				distributeBedChange(state, dom);

				// do nothing already done in outletScalarFlux in ScalarTransport.h
				break;

			case SWE_BC_Q_CONST:
			case SWE_BC_HZ_T_INLET:
			case SWE_BC_Q_T:
			case SWE_BC_TRIGGERED_IRRIGATION:
				// to be defined for inlet boundaries with solutes (need an input file)
				break;
			}
		}
#endif

		if (bctype == SWE_BC_TRIGGERED_IRRIGATION)
		{
			extraMass += applyTriggeredIrrigation(state, dom);
		}
		else if (bctype == SWE_BC_TRIGGERED_DRAINAGE)
		{
			extraMass += applyTriggeredDrainage(state, dom);
		}
		else if (ncellsBC > 0)
		{
			switch (bctype)
			{
			default:
				std::cerr << RERROR "Boundary type: " << bctype << " not recognised." << std::endl;
				std::cerr << RERROR "No boundary condition applied." << std::endl;
				exit(EXIT_FAILURE);
				break;

			case SWE_BC_CRITICAL: // critical flow boundary condition
#if SERGHEI_DEBUG_BOUNDARY
				std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "SWE_BC_CRITICAL" << std::endl;
#endif
				Kokkos::parallel_for("swe_bc_critical", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob) {
		      	int ii=bcells[iGlob];
		        real h=state.h(ii);
		        if( h>=state.hmin) {
            	real hu=state.hu(ii);
		          real hv=state.hv(ii);
							//----------------------------
							real modQ=sqrt(hu*hu+hv*hv);
							real vel=modQ/h;
							if(vel/sqrt(GRAV*h)<1.){
								vel=1.0*sqrt(GRAV*h); //Froude 1.0 (critical)
							}
		          hu=vel*h*normalx;
		          hv=vel*h*normaly;
		          state.hu(ii)=hu;
		          state.hv(ii)=hv;
		        } });
				break;

			case SWE_BC_WSE_CONST: // constant free surface elevation
#if SERGHEI_DEBUG_BOUNDARY
				std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "SWE_BC_WSE_CONST" << std::endl;
#endif
				Kokkos::parallel_reduce("swe_bc_swe_const", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &sumM) {
            int ii = bcells[iGlob];
            real h = state.h(ii);
						if(h > TOL12){
							real hu= state.hu(ii);
							real hv= state.hv(ii);
							real z = state.z(ii);
							state.h(ii) = max(bcvals(0) - z, (real) 0.0); // enforce water depth positivity
							sumM += (state.h(ii)-h)*dom.cellArea();
							//orientation wrt to the outflow normal direction
							real modQ=sqrt(hu*hu+hv*hv);
							state.hu(ii)=normalx*modQ;
							state.hv(ii)=normaly*modQ;
						} }, Kokkos::Sum<real>(extraMass));
				break;

			case SWE_BC_H_CONST: // constant depth boundary condition
#if SERGHEI_DEBUG_BOUNDARY
				std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "SWE_BC_H_CONST" << std::endl;
#endif
				Kokkos::parallel_reduce("swe_bc_h_const", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &sumM) {
	          int ii = bcells[iGlob];
			      real h = state.h(ii);
						if(h>TOL12){
							real hu= state.hu(ii);
							real hv= state.hv(ii);
							state.h(ii) = bcvals(0);
							sumM += (state.h(ii)-h)*dom.cellArea();
							//orientation wrt to the outflow normal direction
							real modQ=sqrt(hu*hu+hv*hv);
							state.hu(ii)=normalx*modQ;
							state.hv(ii)=normaly*modQ;
						} }, Kokkos::Sum<real>(extraMass));
				break;

			case SWE_BC_Q_CONST: // constant inflow discharge boundary condition
#if SERGHEI_DEBUG_BOUNDARY
				std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "SWE_BC_Q_CONST" << std::endl;
#endif
				Kokkos::parallel_reduce("swe_bc_q_const", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &sumM) {
          int ii = bcells[iGlob];
          real h = state.h(ii);
          real qbc_x = bcvals(1);
 			    real qbc_y = bcvals(2);

			    //check if water depth is subcritical (Froude
			    //number less than 0.99). Otherwise impose
			    //boundary condition water depth
			    real inletFr=0.99;
			    //real hcr=cbrt((qbc_x*qbc_x+qbc_y*qbc_y)/(GRAV*inletFr*inletFr));
          real hcr = criticalDepth(qbc_x,qbc_y,inletFr);

			    //TODO: add this mass in case h is less than hcr
			    // If the current depth is below the critical depth, force it up to critical
			    if (hcr > h) state.h(ii) = max(h, hcr);

          sumM += (state.h(ii)-h)*dom.cellArea();

			    state.hu(ii) = qbc_x;
			    state.hv(ii) = qbc_y; }, Kokkos::Sum<real>(extraMass));
				break;

			case SWE_BC_FREE_OUTFLOW: // zero gradient or free boundary
#if SERGHEI_DEBUG_BOUNDARY
				std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "SWE_BC_FREE_OUTFLOW" << std::endl;
#endif
				Kokkos::parallel_for("swe_bc_free", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob) {
            int ii=bcells[iGlob];
				    //orientation wrt to the outflow direction
						real hu= state.hu(ii);
						real hv= state.hv(ii);
						real modQ=sqrt(hu*hu+hv*hv);
						state.hu(ii)=normalx*modQ;
						state.hv(ii)=normaly*modQ; });
				break;

			case SWE_BC_HZ_T_INLET: // stage hydrograph inlet
			{
#if SERGHEI_DEBUG_BOUNDARY
				std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "SWE_BC_HZ_T_INLET" << std::endl;
#endif
				real hzBC = interpolateLinear(hydrograph, dom.etime);
#if SERGHEI_DEBUG_BOUNDARY
				std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "time = " << dom.etime << "\th+z = " << hzBC << std::endl;
#endif
				Kokkos::parallel_reduce("swe_bc_hz_t_inlet", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &sumM) {
					int ii = bcells[iGlob];
					real h = state.h(ii);
					real hu = state.hu(ii);
					real hv = state.hv(ii);

					state.h(ii) = max(hzBC - state.z(ii), 0.0);
					sumM += (state.h(ii) - h) * dom.cellArea();

					/*
							//orientation wrt to the outflow normal direction
								real modQ=sqrt(hu*hu+hv*hv);
						state.hu(ii)=normalx*modQ;
						state.hv(ii)=normaly*modQ;
					*/ }, Kokkos::Sum<real>(extraMass));
				break;
			}

			case SWE_BC_HZ_T_OUTLET: // stage hydrograph outlet
			{
#if SERGHEI_DEBUG_BOUNDARY
				std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "SWE_BC_HZ_T_OUTLET" << std::endl;
#endif
				real hzBC = interpolateLinear(hydrograph, dom.etime);
#if SERGHEI_DEBUG_BOUNDARY
				std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "time = " << dom.etime << "\th+z = " << hzBC << std::endl;
#endif
				Kokkos::parallel_reduce("swe_bc_hz_t_outlet", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &sumM) {
					int ii = bcells[iGlob];
					real h = state.h(ii);
					state.h(ii) = max(hzBC - state.z(ii), (real)0.0);
					sumM += (state.h(ii) - h) * dom.cellArea();

					// no orientation wrt to the outflow normal direction to allow tidal wave coming into the domain
				},
										Kokkos::Sum<real>(extraMass));
				break;
			}

			case SWE_BC_Q_T: // hydrograph
			{
#if SERGHEI_DEBUG_BOUNDARY
				std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "SWE_BC_Q_T" << std::endl;
#endif
				real Q = interpolateLinear(hydrograph, dom.etime);
#if SERGHEI_DEBUG_BOUNDARY
				std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "HYDROGRAPH BC: time = " << dom.etime << "\tQ = " << Q << std::endl;
#endif
				flattenWaterSurface(state, dom);
				distributeDischarge(state, dom, Q);
				// extraMass = Q;
				break;
			}
			} // end switch
		} // endif ncellsBC

#if SERGHEI_SCALAR_TRANSPORT
		if (ncellsBC > 0)
		{
			switch (bctype)
			{
			// outlet boundary conditions
			case SWE_BC_CRITICAL:
			case SWE_BC_H_CONST:
			case SWE_BC_WSE_CONST:
			case SWE_BC_FREE_OUTFLOW:
			case SWE_BC_HZ_T_OUTLET:
			case SWE_BC_TRIGGERED_DRAINAGE:

				// do nothing already done in outletScalarFlux in ScalarTransport.h
				break;

			case SWE_BC_Q_CONST:
			case SWE_BC_HZ_T_INLET:
			case SWE_BC_Q_T:
			case SWE_BC_TRIGGERED_IRRIGATION:
				// to be defined for inlet boundaries with solutes (need an input file)
				break;
			}
		}
#endif

		adjustedVolume = extraMass; // extraMass per bc
		Kokkos::fence("bc.apply-end");
			dom.timers.swe.bc.total += timer.seconds();
		}

#if CROP_GROWTH_MODEL
		inline void apply(State &state, Domain &dom, Wofost72 &wofost)
		{
			if (bctype != SWE_BC_TRIGGERED_IRRIGATION && bctype != SWE_BC_TRIGGERED_DRAINAGE)
			{
				apply(state, dom);
				return;
			}

#if SERGHEI_DEBUG_WORKFLOW
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
#endif
			Kokkos::fence("bc.apply-crop-begin");
			Kokkos::Timer timer;
			real extraMass = 0.0;
			if (bctype == SWE_BC_TRIGGERED_IRRIGATION)
				extraMass = applyTriggeredIrrigation(state, dom, wofost);
			else
				extraMass = applyTriggeredDrainage(state, dom, wofost);
			adjustedVolume = extraMass;
			Kokkos::fence("bc.apply-crop-end");
			dom.timers.swe.bc.total += timer.seconds();
		}
#endif

		inline void integrate(State &state, Domain &dom)
	{
#if SERGHEI_DEBUG_WORKFLOW
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
#endif
		Kokkos::fence("bc.integrate-begin");
		Kokkos::Timer timer;
		real outDischarge = 0.0;
		real inDischarge = 0.0;
		real totalDischarge = 0.0;

#if SERGHEI_SUSPENDED_SEDIMENT
		real outSolidDischarge = 0.0;
		real inSolidDischarge = 0.0;
#endif
		real totalSolidDischarge = 0.0;

		if (ncellsBC > 0)
		{
			// discharge integration
			Kokkos::parallel_reduce("reduceDischargeBC", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &sumD, real &sumSD) {
				int ii = bcells[iGlob];
				real qbound;
		    if( state.h(ii)>=state.hmin) {
					//the integration is done over all boundary walls according to the outflow direction
					qbound= (state.hu(ii)*sgn(normalx) + state.hv(ii)*sgn(normaly));
					sumD +=  qbound * dom.dx();
#if SERGHEI_SUSPENDED_SEDIMENT
						for(int iphi=state.sediment.iphised; iphi<(state.sediment.iphised+state.sediment.nSed); iphi++){
							sumSD += qbound * state.ade.hphi(ii,iphi)/state.h(ii) * dom.dx();
						}
#endif
				} }, Kokkos::Sum<real>(totalDischarge), Kokkos::Sum<real>(totalSolidDischarge));
			Kokkos::fence();

			switch (bctype)
			{
			case SWE_BC_CRITICAL:
			case SWE_BC_H_CONST:
			case SWE_BC_WSE_CONST:
			case SWE_BC_FREE_OUTFLOW:
			case SWE_BC_HZ_T_OUTLET:
			case SWE_BC_TRIGGERED_DRAINAGE:
				outDischarge = totalDischarge;
#if SERGHEI_SUSPENDED_SEDIMENT
				outSolidDischarge = totalSolidDischarge;
#endif
				break;
			case SWE_BC_Q_CONST:
			case SWE_BC_HZ_T_INLET:
			case SWE_BC_Q_T:
			case SWE_BC_TRIGGERED_IRRIGATION:
				inDischarge = totalDischarge;
#if SERGHEI_SUSPENDED_SEDIMENT
				inSolidDischarge = totalSolidDischarge;
#endif
				break;
			default:
				std::cerr << RERROR " Boundary type: " << bctype << " not recognized." << std::endl;
				std::cerr << RERROR " No boundary condition applied." << std::endl;
				exit(EXIT_FAILURE);
			}
#if SERGHEI_DEBUG_BOUNDARY
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "totalDischarge = " << totalDischarge << std::endl;
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "inDischarge = " << inDischarge << std::endl;
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "outDischarge = " << outDischarge << std::endl;
#if SERGHEI_SUSPENDED_SEDIMENT
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "totalSolidDischarge = " << totalSolidDischarge << std::endl;
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "inSolidDischarge = " << inSolidDischarge << std::endl;
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "outSolidDischarge = " << outSolidDischarge << std::endl;
#endif
#endif
		}
		inflowDischarge = inDischarge;				 // inflowdischarge local per bc
		outflowDischarge = outDischarge;			 // outflowdischarge local per bc
		inflowAccumulated += inDischarge * dom.dt;	 // inflowAccumulated local per bc
		outflowAccumulated += outDischarge * dom.dt; // outflowaccumulated local per bc
#if SERGHEI_SUSPENDED_SEDIMENT
		inflowSolidDischarge = inSolidDischarge;			   // inflowdischarge local per bc
		outflowSolidDischarge = outSolidDischarge;			   // outflowdischarge local per bc
		inflowSolidAccumulated += inSolidDischarge * dom.dt;   // inflowAccumulated local per bc
		outflowSolidAccumulated += outSolidDischarge * dom.dt; // outflowaccumulated local per bc
#endif

		Kokkos::fence("bc.integrate-end");
		dom.timers.swe.bc.integrate += timer.seconds();
	}

	inline void reduce(Parallel const &par)
	{
		// only used to write out to file
		real Qin, Qout;

		SERGHEI_MPI_CHECK(MPI_Reduce(&inflowDischarge, &Qin, 1, SERGHEI_MPI_REAL, MPI_SUM, SERGHEI_MASTERPROC, MPI_COMM_WORLD));
		SERGHEI_MPI_CHECK(MPI_Reduce(&outflowDischarge, &Qout, 1, SERGHEI_MPI_REAL, MPI_SUM, SERGHEI_MASTERPROC, MPI_COMM_WORLD));
		netQ = Qin - Qout;

		SERGHEI_MPI_CHECK(MPI_Reduce(&inflowAccumulated, &Qin, 1, SERGHEI_MPI_REAL, MPI_SUM, SERGHEI_MASTERPROC, MPI_COMM_WORLD));
		SERGHEI_MPI_CHECK(MPI_Reduce(&outflowAccumulated, &Qout, 1, SERGHEI_MPI_REAL, MPI_SUM, SERGHEI_MASTERPROC, MPI_COMM_WORLD));
		netVol = Qin - Qout;
	}

	void setIsBound(State &state, int value) const
	{
		Kokkos::parallel_for("init_isBound", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob) {
	 		int ii = bcells[iGlob]; //extended domain index
			state.isBound(ii)=value; });
	}

	// ============================================================================
// Surface solute transport boundary condition application
// ============================================================================
#if SERGHEI_SURFACE_TRANSPORT
	inline void applyrtbc(RTStateSW &rtsw, const Domain &dom, const State &state)
	{
		if (ncellsBC <= 0)
			return;
		// Initialise the flux bookkeeping arrays
		if (spec_mass_in_rate.size() != rtsw.n_mass)
		{
			spec_mass_in_rate.resize(rtsw.n_mass, 0.0);
			spec_mass_out_rate.resize(rtsw.n_mass, 0.0);
		}

		for (int iSpec = 0; iSpec < rtsw.n_mass; iSpec++)
		{
			int current_bctype = bctype;
			bool has_rt_config = (current_bctype == RTSW_BC_NONE ||
								  current_bctype == RTSW_BC_CON_CONST ||
								  current_bctype == RTSW_BC_ZEROGRAD ||
								  current_bctype == RTSW_BC_CON_T ||
								  current_bctype == RTSW_BC_TRIGGERED_IRRIGATION ||
								  current_bctype == RTSW_BC_TRIGGERED_DRAINAGE);
			if (iSpec < spec_bctype.size() && spec_bctype[iSpec] != -999)
			{
				current_bctype = spec_bctype[iSpec];
				has_rt_config = true;
			}

			// Concentration is a scalar.  A reflective hydraulic wall therefore
			// uses even reflection (zero normal concentration gradient): both RT
			// ghost layers copy the adjacent interior value instead of being set
			// to a zero-concentration Dirichlet value.
			if (!has_rt_config && hydro_bctype == SWE_BC_REFLECTIVE)
			{
				current_bctype = RTSW_BC_ZEROGRAD;
			}
			else if (!has_rt_config && hydro_bctype == SWE_BC_TRIGGERED_IRRIGATION)
			{
				current_bctype = RTSW_BC_TRIGGERED_IRRIGATION;
				if (!rtIrrigationFallbackWarned)
				{
					std::cerr << YEXC << "SWRTM: Triggered irrigation boundary has no RT concentration BC; using gated clear-water fallback (Dirichlet 0.0 while irrigation is active, zero-gradient while inactive)." << std::endl;
					rtIrrigationFallbackWarned = 1;
				}
			}
			else if (!has_rt_config && hydro_bctype == SWE_BC_TRIGGERED_DRAINAGE)
			{
				current_bctype = RTSW_BC_TRIGGERED_DRAINAGE;
				if (!rtDrainageFallbackWarned)
				{
					std::cerr << YEXC << "SWRTM: Triggered drainage boundary has no RT concentration BC; using zero-gradient outflow fallback (type RTSW_BC_TRIGGERED_DRAINAGE)." << std::endl;
					rtDrainageFallbackWarned = 1;
				}
			}

			// Skip NONE
			if (current_bctype == RTSW_BC_NONE)
				continue;

			// 门控：触发式类型按配对水动力边界的运行时开关状态解析为有效类型。
			// 灌溉活跃 -> Dirichlet（有时间序列用 CON_T，否则 CON_CONST）；不活跃 -> 零梯度。
			// 排水类型恒为零梯度出流（有水流出才带出质量，不做浓度强制）。
			bool hasTs = (iSpec < spec_ts.size() && spec_ts[iSpec].np > 0);
			int eff_bctype = current_bctype;
			if (current_bctype == RTSW_BC_TRIGGERED_IRRIGATION)
				eff_bctype = hydroTriggerActive ? (hasTs ? RTSW_BC_CON_T : RTSW_BC_CON_CONST)
												: RTSW_BC_ZEROGRAD;
			else if (current_bctype == RTSW_BC_TRIGGERED_DRAINAGE)
				eff_bctype = RTSW_BC_ZEROGRAD;

			// 提取输入的浓度值
			real bc_val = 0.0;
			if (eff_bctype == RTSW_BC_CON_CONST)
			{
				if (iSpec < spec_bcval_const.size())
					bc_val = spec_bcval_const[iSpec];
				else if (bcvals.extent(0) > 0)
					bc_val = bcvals(0);
			}
			else if (eff_bctype == RTSW_BC_CON_T)
			{
				if (hasTs)
				{
					bc_val = interpolateLinear(spec_ts[iSpec], dom.etime);
				}
			}

			real dt = rtsw.dt;
			real dx = dom.dxConst;

			// Local reduction variables for the total in/out mass rate of each species
			real current_in_rate = 0.0;
			real current_out_rate = 0.0;

			Kokkos::parallel_reduce("rt_bc_fvm", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &sum_in, real &sum_out) {
				int ii = bcells[iGlob]; 
				int rt_ii = rtsw.hydroToRtIndex(ii, dom.hc);
				int nxhc = dom.nx + 2 * dom.hc;
				int row = ii / nxhc;
				int col = ii % nxhc;

				real hu = state.hu(ii);
				real hv = state.hv(ii);
				
				real c_int = rtsw.c(iSpec, rt_ii, 0);
				real c_latest = rtsw.c(iSpec, rt_ii, 1);
				real h_curr = state.h(ii);
				
				// --------------------------------------------------------------------
				// Minimum water depth for boundary concentration reconstruction; below it
				// the 1/h term amplifies numerical noise. Recommended 1e-4 m (0.1 mm), or
				// the maximum with the hydrodynamic wet-cell threshold state.hmin as here
				// --------------------------------------------------------------------
				real rt_hmin = state.hmin;

				// --------------------------------------------------------------------
				// Dead-band filter for tiny discharges: high-frequency numerical noise
				// oscillating around zero at the boundary would be rectified by the one-sided
				// sign test (if (hu > 0)) into spurious unidirectional flux spikes
				// --------------------------------------------------------------------
				real q_threshold = 1e-5; // unit-discharge dead-band threshold [m^2/s]
				real hu_filtered = (fabs(hu) > q_threshold) ? hu : 0.0;
				real hv_filtered = (fabs(hv) > q_threshold) ? hv : 0.0;

				real total_q_out = 0.0; // unit discharge leaving the cell across all exposed faces
				real total_q_in  = 0.0; // unit discharge entering the cell across all exposed faces

				// Identify physically exposed faces from the topology and filtered discharges
				// 1. East face (i+1)
				if ((col == dom.nx + dom.hc - 1 && dom.iE) || state.isnodata(ii+1)) {
					if (hu_filtered > 0.0) { total_q_out += hu_filtered; }       // velocity towards +x: water leaves the domain
					else                   { total_q_in  += (-hu_filtered); }    // velocity towards -x: water enters the domain
				}
				
				// 2. West face (i-1)
				if (col == dom.hc  && dom.iW || state.isnodata(ii - 1)) {
					if (hu_filtered < 0.0) { total_q_out += (-hu_filtered); }    // velocity towards -x: water leaves the domain
					else                   { total_q_in  += hu_filtered; }       // velocity towards +x: water enters the domain
				}
				
				// 3. North face (j-1)
				if (row == dom.hc && dom.iN || state.isnodata(ii - nxhc)) {
					if (hv_filtered > 0.0) { total_q_out += hv_filtered; }       // velocity towards +y: water leaves the domain
					else                   { total_q_in  += (-hv_filtered); }    // velocity towards -y: water enters the domain
				}
				
				// 4. South face (j+1)
				if (row == dom.ny + dom.hc - 1 && dom.iS || state.isnodata(ii + nxhc)) {
					if (hv_filtered < 0.0) { total_q_out += (-hv_filtered); }    // velocity towards -y: water leaves the domain
					else                   { total_q_in  += hv_filtered; }       // velocity towards +y: water enters the domain
				}

				// Mass flux statistics [mg/s]
				if (h_curr >= rt_hmin) {
					// Accumulate only above the safe depth threshold; ultra-shallow water distorts the fluxes
					sum_out += total_q_out * c_int * dx * 1000.0;
					
					real c_inflow = (eff_bctype == RTSW_BC_ZEROGRAD) ? 0.0 : bc_val;
					sum_in += total_q_in * c_inflow * dx * 1000.0;
				}

				// Ensure the underlying mass has a sensible initial value
				if (rtsw.solute_mass(iSpec, rt_ii, 1) == 0.0 && c_latest > 0.0 && h_curr >= rt_hmin) {
					rtsw.solute_mass(iSpec, rt_ii, 1) = c_latest * h_curr;
				}

				// ----------------------------------------------------
				// Branch handling: concentration back-calculation and dry-cell truncation
				// ----------------------------------------------------
				if (eff_bctype == RTSW_BC_ZEROGRAD) {
					if (h_curr >= rt_hmin) {
						if (total_q_out > 0.0) {
							// Subtract only the mass actually leaving through the exposed faces
							real mass_leaving = total_q_out * c_int * dt / dx;
							Kokkos::atomic_sub(&rtsw.solute_mass(iSpec, rt_ii, 1), mass_leaving);
						}
						// Back-calculate the concentration
						rtsw.c(iSpec, rt_ii, 1) = rtsw.solute_mass(iSpec, rt_ii, 1) / h_curr;
						if (rtsw.c(iSpec, rt_ii, 1) < 0.0) {
							rtsw.c(iSpec, rt_ii, 1) = 0.0;
							rtsw.solute_mass(iSpec, rt_ii, 1) = 0.0;
						}
					} else {
						// Near-dry cell after recession: zero concentration and mass to truncate noise
						rtsw.c(iSpec, rt_ii, 1) = 0.0;
						rtsw.solute_mass(iSpec, rt_ii, 1) = 0.0;
					}
				}
				else if (eff_bctype == RTSW_BC_CON_CONST || eff_bctype == RTSW_BC_CON_T) {
					if (h_curr >= rt_hmin) {
						// 注意：这里不把 Dirichlet 改写引起的 Δsolute_mass 计入 sum_in/out。
						// 输运核经幽灵单元（值=bc_val）的平流已把入流质量记入 total_q_in 项，
						// 再记 Δ 会重复计数（实测恰好多计一份 Q·c）。
						rtsw.c(iSpec, rt_ii, 1) = bc_val;
						rtsw.solute_mass(iSpec, rt_ii, 1) = bc_val * h_curr;
					} else {
						rtsw.c(iSpec, rt_ii, 1) = 0.0;
						rtsw.solute_mass(iSpec, rt_ii, 1) = 0.0;
					}
				}

				// ----------------------------------------------------
				// Update the ghost cells outside the domain
				// ----------------------------------------------------
				int ghost1 = -1, ghost2 = -1;
				if (row == dom.hc                     && dom.iN) ghost1 = ii - nxhc; 
				else if (row == dom.ny + dom.hc - 1   && dom.iS) ghost1 = ii + nxhc; 
				if (col == dom.hc                     && dom.iW) ghost2 = ii - 1; 
				else if (col == dom.nx + dom.hc - 1   && dom.iE) ghost2 = ii + 1; 

				auto apply_to_ghost = [&](int g_idx) {
					if (g_idx < 0 || g_idx >= dom.nCellMem) return;
					int rt_g_idx = rtsw.hydroToRtIndex(g_idx, dom.hc);
					
					if (h_curr >= rt_hmin) {
						// Wet boundary: derive the exterior ghost concentration from the boundary condition
						if (eff_bctype == RTSW_BC_ZEROGRAD) {
							rtsw.c(iSpec, rt_g_idx, 1) = rtsw.c(iSpec, rt_ii, 1);
							rtsw.c(iSpec, rt_g_idx, 0) = rtsw.c(iSpec, rt_ii, 1);
						} else {
							rtsw.c(iSpec, rt_g_idx, 1) = bc_val;
							rtsw.c(iSpec, rt_g_idx, 0) = bc_val;
						}
					} else {
						// Dry boundary: zero the exterior ghost cells so small numerical errors do not diffuse inward
						rtsw.c(iSpec, rt_g_idx, 1) = 0.0;
						rtsw.c(iSpec, rt_g_idx, 0) = 0.0;
					}
				};
				apply_to_ghost(ghost1);
				apply_to_ghost(ghost2); }, Kokkos::Sum<real>(current_in_rate), Kokkos::Sum<real>(current_out_rate));

			// Store in member variables for external diagnostics
			spec_mass_in_rate[iSpec] = current_in_rate;
			spec_mass_out_rate[iSpec] = current_out_rate;
		}
	}
#endif
};

class ExternalBoundaries
{
	// This class should not be invoked form a parallel region as it contains strings
public:
	std::string BoundaryTypes[15] = {"NONE", "PERIODIC", "REFLECTIVE", "TRANSMISSIVE", "NONE", "CRITICAL", "CONSTANT DEPTH", "CONSTANT INFLOW", "CONSTANT WSELEVATION", "FREE OUTFLOW", "STAGE HYDROGRAPH INLET", "STAGE HYDROGRAPH OUTLET", "HYDROGRAPH", "TRIGGERED IRRIGATION", "TRIGGERED DRAINAGE"};
	std::vector<std::string> id;
	std::vector<ExtBC> extbc;

	std::string RTTypes[7] = {"NONE", "NONE", "DIRICHLET_CONSTANT", "ZEROGRAD", "DIRICHLET_T", "TRIGGERED_IRRIGATION", "TRIGGERED_DRAINAGE"};
	std::vector<ExtBC> rtbc;
};
#endif
