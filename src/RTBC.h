/*建立地下水溶质运移方程的边界条件*/
#ifndef _RTBC_H_
#define _RTBC_H_

#if SERGHEI_SUBSURFACE_TRANSPORT
#include "const.h"
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "GwDomain.h"
#include "GwMPI.h"
#include "GwState.h"
#include "Parallel.h"
#include "RTState.h"


// RT  bc types
#define SUB_RT_BC_NOFLOW 1
#define SUB_RT_BC_Dirichlet_CONST 2
#define SUB_RT_BC_Neumann_CONST 3
#define SUB_RT_BC_Cauchy_CONST 4
#define SUB_RT_BC_Dirichlet_T 5
#define SUB_RT_BC_Neumann_T 6
#define SUB_RT_BC_Cauchy_T 7
#define SUB_RT_BC_SWE 8
#define SUB_RT_BC_FD 9

// subsurface bc directions
#define XPLUS 1
#define XMINUS 2
#define YPLUS 3
#define YMINUS 4
#define ZPLUS 5
#define ZMINUS 6

class RTBC
{

// this class is safe to invoke in a parallel region
public:
	int ncellsBC = 0; //number of bcells
	int ncellsIT = 0;	// number of internal source/sink cells
	intArr bcells, gcells, rtgwcells, icells; //array of indexes of boundary cells
    intArr swgw_type; // type of surface-subsurface exchange: 0: No ponding, no sw-gw exchange, 1: Ponding with large h, 2: Ponding with small h
	int location, rtbctype, isInDomain, direction;
    realArr bcvals, bcdata;
	TimeSeries ts;
	real Qtot, Qinflow, Qoutflow;

	MPI_Comm comm;	// communicator for ranks associated to the BC

	inline int find_bcells(GwState &gw, std::string &rtid, GwDomain &gdom, Parallel &par, int nPoly, realArr &xPoly, realArr &yPoly){
		int foundInSubdom; // to keep track of which subdomains are associated to this boundary
		std::vector<int> tmpbcells; //array of indexes of boundary cells
		std::vector<int> tmpgcells; //array of indexes of ghost cells
		std::vector<int> subdomains;	// keeps track of which subdomains are associated to the BC
		// Loop over the entire domain to find bc cells
		for (int kk = 0; kk < gdom.nz; kk++) {
			for (int jj = 0; jj < gdom.ny; jj++) {
				for (int ii = 0; ii < gdom.nx; ii++) {
					int iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
					foundInSubdom = -1;
          //每个进程在x和y方向上处理的全局网格范围，即par.i_beg、par.i_end、par.j_beg和par.j_end。
		            real xCoord = gdom.xll + ( par.i_beg + ii + 0.5) * gdom.dx;
		            real yCoord = gdom.yll + gdom.ny_glob*gdom.dx - ( par.j_beg + jj + 0.5) * gdom.dx;
		            // (xCoord, yCoord) Check  if the cell is inside the polygon
                      if (geometry::isInsidePoly(nPoly, xPoly, yPoly, xCoord, yCoord)){
		                // If on top/bottom boundary, only the top/bottom layer counts
		                if (direction == ZPLUS) {
							if (kk == gdom.nz-1) {
								tmpbcells.push_back(iGlob);
                                        //push_back是C++标准库中向容器尾部添加新元素的函数。
								tmpgcells.push_back(iGlob+gdom.nxhc*gdom.nyhc);
							}
						}
		                else if (direction == ZMINUS)    {
							if (kk == 0)  {
								tmpbcells.push_back(iGlob);
								tmpgcells.push_back(iGlob-gdom.nxhc*gdom.nyhc);
								//上边界，direction == ZMINUS，ghost cell是iGlob-gdom.nxhc*gdom.nyhc
							}
						}
		                // Otherwise (lateral boundary), all cells in the vertical direction are included
		                // NOTE: In the future, this should be customized to allow only certain vertical layers to be included
		                else {
							tmpbcells.push_back(iGlob);
							if (direction == XPLUS)	{tmpgcells.push_back(iGlob+1);}
							//右边界，direction == XPLUS，ghost cell是iGlob+1
							else if (direction == XMINUS)	{tmpgcells.push_back(iGlob-1);}
							else if (direction == YPLUS)	{tmpgcells.push_back(iGlob+gdom.nxhc);}
							else if (direction == YMINUS)	{tmpgcells.push_back(iGlob-gdom.nxhc);}
						}
		            }
				}
			}
		}

		ncellsBC=int(tmpbcells.size());
		if(ncellsBC > 0) foundInSubdom = par.myrank; // if at least one cell in this subdomain (rank) is in the BC, tag as found

		int ncells_all;
        MPI_Allreduce(&ncellsBC, &ncells_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//使用了MPI的MPI_Allreduce函数，它的作用是将所有进程中的ncellsBC的值求和，并将结果存储在ncells_all中
		int *subdoms;
		subdoms = (int*) malloc(par.nranks * sizeof(int));
          //通过malloc函数动态分配的内存大小为par.nranks * sizeof(int)字节，返回的指针类型为int*。
		MPI_Allgather(&foundInSubdom,1,MPI_INT,subdoms,1,MPI_INT,MPI_COMM_WORLD);
		#if SERGHEI_DEBUG_BOUNDARY
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "GW ncellsBC " << ncells_all << std::endl;
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "GW BC subdomains: " ;
			for (int i=0; i<par.nranks; i++){
				std::cout << GGD << " ";
				if(subdoms[i]==par.myrank) std::cout << RED;
				std::cout << subdoms[i] << "\t"<< RESET ;
			}
			std::cout << std::endl;
		#endif

		for(int i=0; i<par.nranks; i++){
			if(subdoms[i] >= 0){
				subdomains.push_back(subdoms[i]);
			}
		}
		#if SERGHEI_DEBUG_BOUNDARY
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "GW Consolidated BC subdomains = ";
			for(int i=0; i<subdomains.size(); i++){
				std::cout << GGD << " " ;
				std::cout << subdomains[i] << "\t";
			}
 			std::cout << std::endl;
		#endif
		MPI_Group group, subgroup;
		MPI_Comm_group(MPI_COMM_WORLD,&group);
		MPI_Group_incl(group,subdomains.size(),subdomains.data(),&subgroup);
		MPI_Comm_create(MPI_COMM_WORLD,subgroup,&comm);

		//we need the total boundary cells detected by all subdomain to launch an error otherwise
		if(ncells_all>0){
			bcells=intArr("bcells", ncellsBC);
			gcells=intArr("gcells", ncellsBC);
			rtgwcells = intArr("rtgwcells", ncellsBC);
			
            if (rtbctype == SUB_RT_BC_SWE)   {swgw_type = intArr("swgw_type", ncellsBC);}
			#ifdef __NVCC__
				cudaMemcpyAsync( bcells.data() , tmpbcells.data() , ncellsBC*sizeof(int) , cudaMemcpyHostToDevice );
				cudaMemcpyAsync( gcells.data() , tmpgcells.data() , ncellsBC*sizeof(int) , cudaMemcpyHostToDevice );
				cudaMemcpyAsync( rtgwcells.data() , tmpgcells.data() , ncellsBC*sizeof(int) , cudaMemcpyHostToDevice );
				cudaDeviceSynchronize();
			#else
				std::memcpy(bcells.data(), tmpbcells.data(), ncellsBC*sizeof(int));
				std::memcpy(gcells.data(), tmpgcells.data(), ncellsBC*sizeof(int));
				std::memcpy(rtgwcells.data(), tmpgcells.data(), ncellsBC*sizeof(int));
			#endif
		}
		else{
			if(par.masterproc){
				std::cerr << RERROR << "No boundary cells found for subsurface boundary with id '" << rtid << "'" << std::endl;
			}
			return 0;
		}
		return 1;
	}

    // Apply subsurface boundary conditions
    inline void applyConcentrationBC(RTState &rt, GwState &gw, GwDomain &gdom, Parallel &par) {
		// Check if on global boundaries
		bool onBoundary = 0;
		if (direction == 2 && par.px == 0)	{onBoundary = 1;}
		else if (direction == 1 && par.px == par.nproc_x-1)	{onBoundary = 1;}
		else if (direction == 4 && par.py == 0)	{onBoundary = 1;}
		else if (direction == 3 && par.py == par.nproc_y-1)	{onBoundary = 1;}
		else if (direction == 5 || direction == 6)	{onBoundary = 1;}

        // Kokkos::Timer timer;
        if (ncellsBC > 0 && onBoundary == 1) {
            real Conbc;
			// interpolate if time-series boundary value is read
            if (rtbctype == SUB_RT_BC_Dirichlet_T ) {Conbc = interpolateLinear(ts, gdom.etime);}
            // zero gradient if Q BC is specified
          //   if (bctype == SUB_RT_BC_Q_CONST || bctype == SUB_RT_BC_Q_T || bctype == SUB_RT_BC_FD)   {
          //       Kokkos::parallel_for("gw_bc_h", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
          //           int iGlob = bcells[ibc], iGhost = gcells[ibc];
		// 			gw.h(iGhost,1) = gw.h(iGlob,1);
          //       });
          //   }
			//一类边界 CONST or 二类边界 CONST or 三类边界CONST
          else if (rtbctype == SUB_RT_BC_Dirichlet_CONST || rtbctype == SUB_RT_BC_Neumann_CONST || rtbctype == SUB_RT_BC_Cauchy_CONST){
				Kokkos::parallel_for("rt_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
                    int iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc]; //20240510
                         rt.c(rtiGhost,1) = bcvals(ibc);
			// std::cout << "rt.c(iGhost,1) = " << rt.c(rtiGhost,1) << std::endl;	
                });
			}
            // H Time series
            else if (rtbctype == SUB_RT_BC_Dirichlet_T)  {
                Kokkos::parallel_for("rt_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
                    int iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc]; //20240510
                    rt.c(rtiGhost,1) = Conbc;
                });
            }
            // Prescribed water table BC
          //   else if (bctype == SUB_RT_BC_WT_T)  {
          //       Kokkos::parallel_for("gw_bc_wt", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
          //           int ivg, iGlob = bcells[ibc], iGhost = gcells[ibc];
          //           real wcs, wcr, alpha, n;
          //           ivg = gw.soilID(iGlob) * gw.nVGparam;
          //           wcs = gw.vgTable(ivg+2);    wcr = gw.vgTable(ivg+3);
          //           n = gw.vgTable(ivg+4);  alpha = gw.vgTable(ivg+6);
		// 			gw.h(iGhost,1) = hbc - gdom.z(iGlob);
		// 			gw.wc(iGhost,1) = h2wc(gw.h(iGhost,1), alpha, n, wcs, wcr);
          //       });
          //   }
            // Surface-subsurface exchange
          //   else if (rtbctype == SUB_RT_BC_SWE)  {
		// 		#if SERGHEI_SWE_MODEL
		// 		if (direction == 6)	{
		// 			Kokkos::parallel_for("rt_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
	     //                int iGlob = bcells[ibc], iGhost = gcells[ibc], ivg = gw.soilID(iGlob) * NVG;
	     //                real ks = gw.vgTable(ivg);
		// 				int ii, jj, kk, iGlobSW;
		// 				gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
		// 				iGlobSW = jj*gdom.nxhc + ii;
		// 				gw.h(iGhost,1) = gw.hs(iGlobSW);
		// 				// get sw-gw exchange type
		// 				if (gw.h(iGhost,1) > 0.0)    {
		// 					real q_infilt = 2.0 * ks * (gw.h(iGlob,1) - gw.h(iGhost,1)) / gdom.dz(iGlob) - ks;
		// 					if (-q_infilt * gdom.dt <= gw.h(iGhost,1))   {swgw_type(ibc) = 1;}
		// 					else {swgw_type(ibc) = 2;}
		// 				}
		// 				else {
		// 					// exfiltration
		// 					if (gw.h(iGlob,1) > gw.h(iGhost,1) + 0.5*gdom.dz(iGlob)) {swgw_type(ibc) = 1;}
		// 					// no flow
		// 					else {swgw_type(ibc) = 0;}
		// 				}
	     //            });
		// 		}
		// 		else {
		// 			if (par.masterproc)	{std::cerr << RERROR "BC direction must be 6 for SW-GW exchange boundary! " << "\n";}
		// 		}
		// 		#endif
          //   }
        }
        // gdom.timers.gw += timer.seconds();
    }

    // Apply boundary conditions to the matrix coefficients
    inline void applyRTMatBC(RTState &rt, GwState &gw, GwDomain &gdom, Parallel &par) {
		// Check if on global boundaries
		bool onBoundary = 0;
		if (direction == 2 && par.px == 0)	{onBoundary = 1;}
		else if (direction == 1 && par.px == par.nproc_x-1)	{onBoundary = 1;}
		else if (direction == 4 && par.py == 0)	{onBoundary = 1;}
		else if (direction == 3 && par.py == par.nproc_y-1)	{onBoundary = 1;}
		else if (direction == 5 || direction == 6)	{onBoundary = 1;}

        if (ncellsBC > 0 && onBoundary == 1) {
            real Neumannbc;
            if (rtbctype == SUB_RT_BC_Neumann_T) {Neumannbc = interpolateLinear(ts, gdom.etime);}
	  	    switch (rtbctype) {
                case SUB_RT_BC_Dirichlet_CONST:    
                case SUB_RT_BC_Cauchy_CONST:   
                case SUB_RT_BC_Dirichlet_T:    
                case SUB_RT_BC_Cauchy_T:
                    Kokkos::parallel_for("rt_bc_Dirichlet_const", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
                        int ii, jj, kk, idom, iGlobSW, iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc]; //20240510
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
						iGlobSW = jj*gdom.nxhc + ii;
						if (direction == 1)	{
							// rt.RTcoef(idom,1) = rt.RTcoef(idom,1) * 2.0;//为什么要乘以2？
							rt.RTcoef(idom,7) -= rt.RTcoef(idom,1) * rt.c(iGlob+1,1);
						}
						else if (direction == 2)	{
							// rt.RTcoef(idom,2) = rt.RTcoef(idom,2) * 2.0;
							rt.RTcoef(idom,7) -= rt.RTcoef(idom,2) * rt.c(iGlob-1,1);
						}
						else if (direction == 3)	{
							// rt.RTcoef(idom,3) = rt.RTcoef(idom,3) * 2.0;
							rt.RTcoef(idom,7) -= rt.RTcoef(idom,3) * rt.c(iGlob-gdom.nxhc,1);
						}
						else if (direction == 4)	{
							// rt.RTcoef(idom,4) = rt.RTcoef(idom,4) * 2.0;
							rt.RTcoef(idom,7) -= rt.RTcoef(idom,4) * rt.c(iGlob-gdom.nxhc,1);
						}
						else if (direction == 5)	{
							// rt.RTcoef(idom,5) = rt.RTcoef(idom,5) * 2.0;
							rt.RTcoef(idom,7) -= rt.RTcoef(idom,5) * rt.c(iGlob+gdom.nxhc*gdom.nyhc,1);
						}
						else if (direction == 6)	{
							// rt.RTcoef(idom,6) = rt.RTcoef(idom,6) * 2.0;
							rt.RTcoef(idom,7) -= rt.RTcoef(idom,6) * rt.c(iGlob-gdom.nxhc*gdom.nyhc,1);
							// if (gdom.isEvap)	{
							// 	rt.RTcoef(idom,7) -= gdom.dt * gdom.evapRate(iGlobSW) / gdom.dz(iGlob);
							// }
						}
                    });
                    break;
               //  case SUB_RT_BC_SWE:
               //  	#if SERGHEI_SWE_MODEL
			// 		if (direction == 6)	{
			// 			Kokkos::parallel_for("gw_bc_swe_const", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
	          //               int ii, jj, kk, idom, iGlobSW, iGlob = bcells[ibc];
	          //               gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
	          //               idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
			// 				iGlobSW = jj*gdom.nxhc + ii;
			// 				if (swgw_type(ibc) == 2)    {
			// 					real q_infilt = gw.h(iGlob-gdom.nxhc*gdom.nyhc,1) / gdom.dt;
			// 					rt.RTcoef(idom,7) -= gdom.dt * gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) / gdom.dz(iGlob);
			// 					rt.RTcoef(idom,7) += gdom.dt * q_infilt / gdom.dz(iGlob);
			// 					rt.RTcoef(idom,6) = 0.0;
			// 				}
			// 				else if (swgw_type(ibc) == 0)	{
			// 					rt.RTcoef(idom,7) -= gdom.dt * gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) / gdom.dz(iGlob);
			// 					rt.RTcoef(idom,6) = 0.0;
			// 				}
			// 				else {
			// 					rt.RTcoef(idom,6) = rt.RTcoef(idom,6) * 2.0;
			// 					rt.RTcoef(idom,7) -= rt.RTcoef(idom,6) * gw.h(iGlob-gdom.nxhc*gdom.nyhc,1);
			// 				}
			// 				// evaporation
			// 				if (gdom.isEvap)	{
			// 					rt.RTcoef(idom,7) -= gdom.dt * gdom.evapRate(iGlobSW) / gdom.dz(iGlob);
			// 				}
	          //           });
			// 		}
               //      else {
			// 			if (par.masterproc)	{std::cerr << RERROR "BC direction must be 6 for sw-gw boundary! " << "\n";}
			// 		}
               //      #endif
               //      break;
				case SUB_RT_BC_Neumann_CONST:
					Kokkos::parallel_for("rt_bc_Neumann", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
						int ii, jj, kk, idom, iGlobSW, iGlob = bcells[ibc];
						gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
						idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
						iGlobSW = jj*gdom.nxhc + ii;
						if (direction == 1)	{
							rt.RTcoef(idom,7) += gdom.dt * bcvals(ibc);
							rt.RTcoef(idom,1) = 0.0;
						}
						else if (direction == 2)	{
							rt.RTcoef(idom,7) -= gdom.dt * bcvals(ibc) ;
							rt.RTcoef(idom,2) = 0.0;
						}
						else if (direction == 3)	{
							rt.RTcoef(idom,7) += gdom.dt * bcvals(ibc) ;
							rt.RTcoef(idom,3) = 0.0;
						}
						else if (direction == 4)	{
							rt.RTcoef(idom,7) -= gdom.dt * bcvals(ibc);
							rt.RTcoef(idom,4) = 0.0;
						}
						else if (direction == 5)	{
							rt.RTcoef(idom,5) = 0.0;
						}
						else if (direction == 6)	{
							rt.RTcoef(idom,7) -= gdom.dt * gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) ;
							rt.RTcoef(idom,7) -= gdom.dt * bcvals(ibc) ;
							// if (gdom.isRain)    {rt.RTcoef(idom,7) += gdom.dt * gdom.rainRate(iGlobSW) / gdom.dz(iGlob);}
							// if (gdom.isEvap)	{rt.RTcoef(idom,7) -= gdom.dt * gdom.evapRate(iGlobSW) / gdom.dz(iGlob);}
							rt.RTcoef(idom,6) = 0.0;
						}
					});
					break;
               //  case SUB_RT_BC_Q_T:    case SUB_RT_BC_FD:
               //      Kokkos::parallel_for("gw_bc_q_fd", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
			// 		// for (int ibc = 0; ibc < ncellsBC; ibc++)	{
               //          int ii, jj, kk, idom, iGlobSW, iGlob = bcells[ibc];
               //          gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
               //          idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
			// 			iGlobSW = jj*gdom.nxhc + ii;
			// 			if (direction == 1)	{
			// 				rt.RTcoef(idom,7) += gdom.dt * qbc / gdom.dx;
			// 				rt.RTcoef(idom,1) = 0.0;
			// 			}
			// 			else if (direction == 2)	{
			// 				rt.RTcoef(idom,7) -= gdom.dt * qbc / gdom.dx;
			// 				rt.RTcoef(idom,2) = 0.0;
			// 			}
			// 			else if (direction == 3)	{
			// 				rt.RTcoef(idom,7) += gdom.dt * qbc / gdom.dx;
			// 				rt.RTcoef(idom,3) = 0.0;
			// 			}
			// 			else if (direction == 4)	{
			// 				rt.RTcoef(idom,7) -= gdom.dt * qbc / gdom.dx;
			// 				rt.RTcoef(idom,4) = 0.0;
			// 			}
			// 			else if (direction == 5)	{
			// 				rt.RTcoef(idom,5) = 0.0;
			// 			}
			// 			else if (direction == 6)	{
			// 				rt.RTcoef(idom,7) -= gdom.dt * gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) / gdom.dz(iGlob);
			// 				rt.RTcoef(idom,7) -= gdom.dt * qbc / gdom.dz(iGlob);
			// 				// if (gdom.isRain)    {rt.RTcoef(idom,7) += gdom.dt * gdom.rainRate(iGlobSW) / gdom.dz(iGlob);}
			// 				// if (gdom.isEvap)	{rt.RTcoef(idom,7) -= gdom.dt * gdom.evapRate(iGlobSW) / gdom.dz(iGlob);}
			// 				rt.RTcoef(idom,6) = 0.0;
			// 			}
               //      });
			// 		// }
               //      break;
            }
        }
    }
};

class RTSubsurfaceBoundaries{
// This class should not be invoked form a parallel region as it contains strings
public:
  	std::string BoundaryTypes[9] = {"NOFLOW","CONST_H","CONST_Q","CONST_WT","H_TIMESERIES","Q_TIMESERIES","WT_TIMESERIES","SWEXCHANGE","FREE_DRAINAGE"};
	std::vector<std::string> id;
	std::vector<RTBC> rtgwbc;
};
#endif

#endif

