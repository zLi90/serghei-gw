/* -*- mode: c++ -*- */

#ifndef _BC_H_
#define _BC_H_

#include "define.h"
#include "Indexing.h"

// outer boundary direction definitions
#define SERGHEI_BC_OUTER_REFLECTIVE_CELL 9999
#if SERGHEI_SURFACE_TRANSPORT
#include "RTStateSW.h"
// ==========================================
// 地表水溶质运移边界条件
// ==========================================
#define RTSW_BC_NONE 1		// 无变化 (不干预)
#define RTSW_BC_CON_CONST 2 // 浓度定值 (恒定输入)
#define RTSW_BC_ZEROGRAD 3	// 零浓度梯度 (自由流出)
#define RTSW_BC_CON_T 4		// 浓度时间序列 (动态输入)
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

KOKKOS_INLINE_FUNCTION real criticalDepth(real hu, real hv, real Fr)
{
	return (cbrt((hu * hu + hv * hv) / (GRAV * Fr * Fr)));
}

class ExtBC
{
	// this class is safe to invoke in a parallel region
public:
	int ncellsBC = 0;	   // number of bcells
	intArr bcells;		   // array of indexes of boundary cells
	real normalx, normaly; // direction set by user for inflow/outflow
	int location;		   // 1->west, 2->north, 3->east, 4-> south
	int bctype;
	int isInDomain;
	realArr bcvals;
	real outflowDischarge;
	real outflowAccumulated = 0;
	real inflowDischarge;
	real inflowAccumulated = 0;
	real adjustedVolume = 0;

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

	// 多组分边界条件支持
	std::vector<int> spec_bctype;		  // 每个组分的边界类型
	std::vector<realArr> spec_bcvals;	  // 每个组分的边界值数组（空间分布）
	std::vector<std::string> spec_bcfile; // 每个组分的边界文件名
	std::vector<TimeSeries> spec_ts;	  // 每个组分的时间序列
	std::vector<real> spec_bcval_const;	  // 每个组分的常量边界值
	std::vector<int> has_bcfile;		  // 每个组分是否有边界文件
	std::vector<real> spec_mass_in_rate;  // 溶质质量边界通量 (mg/s)
	std::vector<real> spec_mass_out_rate; // 溶质质量边界通量 (mg/s)

public:
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
		MPI_Allreduce(&ncellsBC, &ncells_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if (par.nranks > 1)
			dom.timers.swe.bc.mpi += timermpi.seconds();

		int *subdoms;
		subdoms = (int *)malloc(par.nranks * sizeof(int));
		timermpi.reset();
		MPI_Allgather(&foundInSubdom, 1, MPI_INT, subdoms, 1, MPI_INT, MPI_COMM_WORLD);
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
		MPI_Comm_group(MPI_COMM_WORLD, &group);
		MPI_Group_incl(group, subdomains.size(), subdomains.data(), &subgroup);
		MPI_Comm_create(MPI_COMM_WORLD, subgroup, &comm);
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

	void inline getMinBedElevation(Domain const &dom, State &state)
	{
		Kokkos::Timer timermpi;
		Kokkos::parallel_reduce("swe_bc_z_min", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &zMin) {
			int ii = bcells[iGlob];
				real z = state.z(ii);
				zMin = min(zMin,z); }, Kokkos::Min<real>(zMin));
		real zMin_all;
		timermpi.reset();
		MPI_Allreduce(&zMin, &zMin_all, 1, SERGHEI_MPI_REAL, MPI_MIN, comm);
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
		MPI_Allreduce(&zMax, &zMax_all, 1, SERGHEI_MPI_REAL, MPI_MAX, comm);
		if (dom.nsubdom > 1)
			dom.timers.swe.bc.mpi += timermpi.seconds();
		zMax = zMax_all;

		nzMin = 0;
		Kokkos::parallel_reduce("swe_bc_z_min_count", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, int &nzMin) {
		  int ii = bcells[iGlob];
			if(state.z(ii) == zMin) nzMin++; }, Kokkos::Sum<int>(nzMin));
		int nzMin_all;
		timermpi.reset();
		MPI_Allreduce(&nzMin, &nzMin_all, 1, MPI_INT, MPI_SUM, comm);
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
		MPI_Allreduce(&hzMin, &hzMin_all, 1, SERGHEI_MPI_REAL, MPI_MIN, comm);
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
		MPI_Allreduce(&extraVol, &extraVol_all, 1, SERGHEI_MPI_REAL, MPI_SUM, comm);
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
			MPI_Allreduce(&extraArea, &extraArea_all, 1, SERGHEI_MPI_REAL, MPI_SUM, comm);
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
		// 1. 计算边界总水深，只统计有水区域
		Kokkos::parallel_reduce("swe_bc_h_weight", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &hsum) {
			int ii = bcells[iGlob];
			if (state.h(ii) > state.hmin) {
				hsum += state.h(ii);
			} }, Kokkos::Sum<real>(hsum));

		real hsum_all;
		timermpi.reset();
		MPI_Allreduce(&hsum, &hsum_all, 1, SERGHEI_MPI_REAL, MPI_SUM, comm);
		if (dom.nsubdom > 1)
			dom.timers.swe.bc.mpi += timermpi.seconds();
		hsum = hsum_all;

		// 2. 计算边界总网格数，用于求解单宽流量
		int ncells_all;
		int local_ncells = ncellsBC;
		timermpi.reset();
		MPI_Allreduce(&local_ncells, &ncells_all, 1, MPI_INT, MPI_SUM, comm);
		if (dom.nsubdom > 1)
			dom.timers.swe.bc.mpi += timermpi.seconds();

		real ds = dom.dx();
		real total_width = ncells_all * ds;
		real q_unit = 0.0;
		if (total_width > 0)
		{
			q_unit = Q / total_width; // 单宽流量
		}

		Kokkos::parallel_for("swe_bc_distribute", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob) {
			int ii = bcells[iGlob];
			real h = state.h(ii);
			real z = state.z(ii);
			real weight = 0;

			// 情况A：边界已经有水 (按照水深占比分配流量，保持水面平滑)
			if (hsum > 1e-6 && h > state.hmin)
			{
				weight = h / hsum;
				state.hu(ii) = Q * weight / ds * normalx;
				state.hv(ii) = Q * weight / ds * normaly;
			}
			// 情况B：边界是干的 -> 强制垫高到临界水深
			else
			{
				// 计算临界水深 hc = (q^2 / g)^(1/3)
				real hc_val = cbrt((q_unit * q_unit) / GRAV);
				if (hc_val < state.hmin)
					hc_val = state.hmin; // 保证大于最小水深

				// 强制修改水深为临界水深，确保动量不被清零
				state.h(ii) = hc_val;
				h = hc_val; // 更新用于 DEBUG 打印的局部变量

				real hu_val = fabs(q_unit);
				state.hu(ii) = hu_val * normalx;
				state.hv(ii) = hu_val * normaly;
			}

#if SERGHEI_DEBUG_BOUNDARY > 1
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "Q = " << Q << "\tii = " << ii << "\tz = " << z << "\t h = " << h << "\tweight = " << weight << "\tds = " << ds << "\t(hu,hv) = " << state.hu(ii) << " " << state.hv(ii) << std::endl;
#endif
		});
	}

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
		MPI_Allreduce(&DZsum, &DZsum_all, 1, SERGHEI_MPI_REAL, MPI_SUM, comm);
		MPI_Allreduce(&Asum, &Asum_all, 1, SERGHEI_MPI_REAL, MPI_SUM, comm);
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

				distributeBedChange(state, dom);

				// do nothing already done in outletScalarFlux in ScalarTransport.h
				break;

			case SWE_BC_Q_CONST:
			case SWE_BC_HZ_T_INLET:
			case SWE_BC_Q_T:
				// to be defined for inlet boundaries with solutes (need an input file)
				break;
			}
		}
#endif

		if (ncellsBC > 0)
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
			    // 修复：如果当前水深小于临界水深，强制垫高到临界水深
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

				// do nothing already done in outletScalarFlux in ScalarTransport.h
				break;

			case SWE_BC_Q_CONST:
			case SWE_BC_HZ_T_INLET:
			case SWE_BC_Q_T:
				// to be defined for inlet boundaries with solutes (need an input file)
				break;
			}
		}
#endif

		adjustedVolume = extraMass; // extraMass per bc
		Kokkos::fence("bc.apply-end");
		dom.timers.swe.bc.total += timer.seconds();
	}

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
				outDischarge = totalDischarge;
#if SERGHEI_SUSPENDED_SEDIMENT
				outSolidDischarge = totalSolidDischarge;
#endif
				break;
			case SWE_BC_Q_CONST:
			case SWE_BC_HZ_T_INLET:
			case SWE_BC_Q_T:
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

		MPI_Reduce(&inflowDischarge, &Qin, 1, SERGHEI_MPI_REAL, MPI_SUM, SERGHEI_MASTERPROC, MPI_COMM_WORLD);
		MPI_Reduce(&outflowDischarge, &Qout, 1, SERGHEI_MPI_REAL, MPI_SUM, SERGHEI_MASTERPROC, MPI_COMM_WORLD);
		netQ = Qin - Qout;

		MPI_Reduce(&inflowAccumulated, &Qin, 1, SERGHEI_MPI_REAL, MPI_SUM, SERGHEI_MASTERPROC, MPI_COMM_WORLD);
		MPI_Reduce(&outflowAccumulated, &Qout, 1, SERGHEI_MPI_REAL, MPI_SUM, SERGHEI_MASTERPROC, MPI_COMM_WORLD);
		netVol = Qin - Qout;
	}

	void setIsBound(State &state, int value) const
	{
		Kokkos::parallel_for("init_isBound", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob) {
	 		int ii = bcells[iGlob]; //extended domain index
			state.isBound(ii)=value; });
	}

	// ============================================================================
// 地表溶质运移边界条件应用
// ============================================================================
#if SERGHEI_SURFACE_TRANSPORT
	inline void applyrtbc(RTStateSW &rtsw, const Domain &dom, const State &state)
	{
		if (ncellsBC <= 0)
			return;
		// 初始化通量记录数组
		if (spec_mass_in_rate.size() != rtsw.n_mass)
		{
			spec_mass_in_rate.resize(rtsw.n_mass, 0.0);
			spec_mass_out_rate.resize(rtsw.n_mass, 0.0);
		}

		for (int iSpec = 0; iSpec < rtsw.n_mass; iSpec++)
		{
			int current_bctype = bctype;
			if (iSpec < spec_bctype.size())
			{
				current_bctype = spec_bctype[iSpec];
			}

			// 跳过 NONE
			if (current_bctype == RTSW_BC_NONE)
				continue;

			// 提取输入的浓度值
			real bc_val = 0.0;
			if (current_bctype == RTSW_BC_CON_CONST)
			{
				if (iSpec < spec_bcval_const.size())
					bc_val = spec_bcval_const[iSpec];
				else if (bcvals.extent(0) > 0)
					bc_val = bcvals(0);
			}
			else if (current_bctype == RTSW_BC_CON_T)
			{
				if (iSpec < spec_ts.size() && spec_ts[iSpec].np > 0)
				{
					bc_val = interpolateLinear(spec_ts[iSpec], dom.etime);
				}
			}

			real dt = rtsw.dt;
			real dx = dom.dxConst;

			// 新增局部规约变量，用于统计每个物质的总输入输出质量率
			real current_in_rate = 0.0;
			real current_out_rate = 0.0;

			Kokkos::parallel_reduce("rt_bc_fvm", ncellsBC, KOKKOS_CLASS_LAMBDA(int iGlob, real &sum_in, real &sum_out) {
				int ii = bcells[iGlob]; 
				int nxhc = dom.nx + 2 * dom.hc;
				int row = ii / nxhc;
				int col = ii % nxhc;

				real hu = state.hu4rtsw(ii, 0);
				real hv = state.hv4rtsw(ii, 0);
				
				// real q_out = hu * normalx + hv * normaly;

				real c_int = rtsw.c(iSpec, ii, 0);
				real c_latest = rtsw.c(iSpec, ii, 1);
				// real h_curr = state.h4rtsw(ii, 1);
				real h_curr = state.h(ii);
				
				// 摒弃 normalx/normaly，基于真实暴露的网格面和实时速度计算通量
				real total_q_out = 0.0; // 当前网格所有暴露面【向外流出】的单宽流量总和
				real total_q_in  = 0.0; // 当前网格所有暴露面【向内流入】的单宽流量总和

// ----------------------------------------------------
				// 摒弃纯 isnodata 判断，结合矩阵拓扑进行物理边界判定
				// ----------------------------------------------------
				
				// 1. 检查东面 (East, i+1)
				if ((col == dom.nx + dom.hc - 1 && dom.iE) || state.isnodata(ii+1)) {
					if (hu > 0.0) { total_q_out += hu; }       // 速度向东(+x)，水流出计算域
					else          { total_q_in  += (-hu); }    // 速度向西(-x)，水流入计算域
				}
				
				// 2. 检查西面 (West, i-1)
				if (col == dom.hc  && dom.iW || state.isnodata(ii - 1)) {
					if (hu < 0.0) { total_q_out += (-hu); }    // 速度向西(-x)，水流出计算域
					else          { total_q_in  += hu; }       // 速度向东(+x)，水流入计算域
				}
				
				// 3. 检查北面 (North, j-1)
				if (row == dom.hc && dom.iN || state.isnodata(ii - nxhc)) {
					if (hv > 0.0) { total_q_out += hv; }       // 速度向北(+y)，水流出计算域
					else          { total_q_in  += (-hv); }    // 速度向南(-y)，水流入计算域
				}
				
				// 4. 检查南面 (South, j+1)
				if (row == dom.ny + dom.hc - 1 && dom.iS || state.isnodata(ii + nxhc)) {
					if (hv < 0.0) { total_q_out += (-hv); }    // 速度向南(-y)，水流出计算域
					else          { total_q_in  += hv; }       // 速度向北(+y)，水流入计算域
				}

				// 质量通量统计 (mg/s)
				if (h_curr > 1e-6) {
					// 真实流出的总质量率 = 暴露面总流出流量 * 内部浓度 * dx * 1000
					sum_out += total_q_out * c_int * dx * 1000.0;
					
					// 真实流入的总质量率 = 暴露面总流入流量 * 边界指定浓度 * dx * 1000
					real c_inflow = (current_bctype == RTSW_BC_ZEROGRAD) ? 0.0 : bc_val;
					sum_in += total_q_in * c_inflow * dx * 1000.0;
				}

				// 确保底层质量有初始值
				if (rtsw.solute_mass(iSpec, ii, 1) == 0.0 && c_latest > 0.0 && h_curr > 1e-6) {
					rtsw.solute_mass(iSpec, ii, 1) = c_latest * h_curr;
				}

				// ----------------------------------------------------
				// 分支处理
				// ----------------------------------------------------
				if (current_bctype == RTSW_BC_ZEROGRAD) {
					if (h_curr > 1e-6) {
						if (total_q_out > 0.0) {
							// 仅扣除真正从各个面“流出”的那部分质量
							real mass_leaving = total_q_out * c_int * dt / dx;
							Kokkos::atomic_sub(&rtsw.solute_mass(iSpec, ii, 1), mass_leaving);
						}
						// 动态反算浓度
						rtsw.c(iSpec, ii, 1) = rtsw.solute_mass(iSpec, ii, 1) / h_curr;
						if (rtsw.c(iSpec, ii, 1) < 0.0) {
							rtsw.c(iSpec, ii, 1) = 0.0;
							rtsw.solute_mass(iSpec, ii, 1) = 0.0;
						}
					} else {
						rtsw.c(iSpec, ii, 1) = 0.0;
						rtsw.solute_mass(iSpec, ii, 1) = 0.0;
					}
				}
				else if (current_bctype == RTSW_BC_CON_CONST || current_bctype == RTSW_BC_CON_T) {
					
					if (h_curr > 1e-6) {
						
						rtsw.c(iSpec, ii, 1) = bc_val;
						rtsw.solute_mass(iSpec, ii, 1) = bc_val * h_curr; 
					} else {
						rtsw.c(iSpec, ii, 1) = 0.0;
						rtsw.solute_mass(iSpec, ii, 1) = 0.0;
					}
				}

				// ----------------------------------------------------
				// 更新域外幽灵网格 (Ghost Cells)
				// ----------------------------------------------------
				int ghost1 = -1, ghost2 = -1;
				if (row == dom.hc                     && dom.iN) ghost1 = ii - nxhc; 
				else if (row == dom.ny + dom.hc - 1   && dom.iS) ghost1 = ii + nxhc; 
				if (col == dom.hc                     && dom.iW) ghost2 = ii - 1; 
				else if (col == dom.nx + dom.hc - 1   && dom.iE) ghost2 = ii + 1; 

				auto apply_to_ghost = [&](int g_idx) {
					if (g_idx < 0 || g_idx >= dom.nCellMem) return;
					
					if (current_bctype == RTSW_BC_ZEROGRAD) {
						// 零梯度流出：让域外浓度始终等于域内边缘浓度
						rtsw.c(iSpec, g_idx, 1) = rtsw.c(iSpec, ii, 1);
						rtsw.c(iSpec, g_idx, 0) = rtsw.c(iSpec, ii, 1);
					} else {
						// 浓度输入：域外维持指定浓度
						rtsw.c(iSpec, g_idx, 1) = bc_val;
						rtsw.c(iSpec, g_idx, 0) = bc_val;
					}
				};
				apply_to_ghost(ghost1);
				apply_to_ghost(ghost2); }, Kokkos::Sum<real>(current_in_rate), Kokkos::Sum<real>(current_out_rate));

			// 保存到成员变量中，供外部读取
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
	std::string BoundaryTypes[13] = {"NONE", "PERIODIC", "REFLECTIVE", "TRANSMISSIVE", "NONE", "CRITICAL", "CONSTANT DEPTH", "CONSTANT INFLOW", "CONSTANT WSELEVATION", "FREE OUTFLOW", "STAGE HYDROGRAPH INLET", "STAGE HYDROGRAPH OUTLET", "HYDROGRAPH"};
	std::vector<std::string> id;
	std::vector<ExtBC> extbc;

	std::string RTTypes[5] = {"NONE", "NONE", "DIRICHLET_CONSTANT", "ZEROGRAD", "DIRICHLET_T"};
	std::vector<ExtBC> rtbc;
};
#endif
