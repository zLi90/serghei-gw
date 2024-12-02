/*建立土壤温度传递方程的边界条件*/
#ifndef _HTBC_H_
#define _HTBC_H_

#if SERGHEI_SUBSURFACE_HEAT
#include "const.h"
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "GwDomain.h"
#include "GwMPI.h"
#include "GwState.h"
#include "Parallel.h"
#include "HTState.h"


// HT  bc types
#define SUB_HT_BC_NOFLOW 1
#define SUB_HT_BC_Dirichlet_CONST 2
#define SUB_HT_BC_Neumann_CONST 3
#define SUB_HT_BC_Cauchy_CONST 4
#define SUB_HT_BC_Dirichlet_T 5
#define SUB_HT_BC_Neumann_T 6
#define SUB_HT_BC_Cauchy_T 7
#define SUB_HT_BC_SWE 8
#define SUB_HT_BC_FD 9

// subsurface bc directions
#define XPLUS 1
#define XMINUS 2
#define YPLUS 3
#define YMINUS 4
#define ZPLUS 5
#define ZMINUS 6

class HTBC
{

// this class is safe to invoke in a parallel region
public:
	int ncellsBC = 0; //number of bcells
	int ncellsIT = 0;	// number of internal source/sink cells
	intArr bcells, gcells, htgwcells, icells; //array of indexes of boundary cells
    intArr swgw_type; // type of surface-subsurface exchange: 0: No ponding, no sw-gw exchange, 1: Ponding with large h, 2: Ponding with small h
	int location, htbctype, isInDomain, direction;
    realArr bcvals, bcdata;
	TimeSeries ts;
	real Qtot, Qinflow, Qoutflow;

	MPI_Comm comm;	// communicator for ranks associated to the BC

	inline int find_bcells(GwState &gw, std::string &htid, GwDomain &gdom, Parallel &par, int nPoly, realArr &xPoly, realArr &yPoly){
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
			htgwcells = intArr("htgwcells", ncellsBC);
			
            if (htbctype == SUB_HT_BC_SWE)   {swgw_type = intArr("swgw_type", ncellsBC);}
			#ifdef __NVCC__
				cudaMemcpyAsync( bcells.data() , tmpbcells.data() , ncellsBC*sizeof(int) , cudaMemcpyHostToDevice );
				cudaMemcpyAsync( gcells.data() , tmpgcells.data() , ncellsBC*sizeof(int) , cudaMemcpyHostToDevice );
				cudaMemcpyAsync( htgwcells.data() , tmpgcells.data() , ncellsBC*sizeof(int) , cudaMemcpyHostToDevice );
				cudaDeviceSynchronize();
			#else
				std::memcpy(bcells.data(), tmpbcells.data(), ncellsBC*sizeof(int));
				std::memcpy(gcells.data(), tmpgcells.data(), ncellsBC*sizeof(int));
				std::memcpy(htgwcells.data(), tmpgcells.data(), ncellsBC*sizeof(int));
			#endif
		}
		else{
			if(par.masterproc){
				std::cerr << RERROR << "No boundary cells found for subsurface heat_transport boundary with id '" << htid << "'" << std::endl;
			}
			return 0;
		}
		return 1;
	}

    // Apply subsurface boundary conditions
    inline void applyTemperatureBC(HTState &ht, GwState &gw, GwDomain &gdom, Parallel &par) {
		// Check if on global boundaries
		bool onBoundary = 0;
		if (direction == 2 && par.px == 0)	{onBoundary = 1;}
		else if (direction == 1 && par.px == par.nproc_x-1)	{onBoundary = 1;}
		else if (direction == 4 && par.py == 0)	{onBoundary = 1;}
		else if (direction == 3 && par.py == par.nproc_y-1)	{onBoundary = 1;}
		else if (direction == 5 || direction == 6)	{onBoundary = 1;}

        // Kokkos::Timer timer;
        if (ncellsBC > 0 && onBoundary == 1) {
            real Tembc;
			// interpolate if time-series boundary value is read
            if (htbctype == SUB_HT_BC_Dirichlet_T ) {Tembc = interpolateLinear(ts, gdom.etime);}
            // zero gradient if Q BC is specified
		  //20240617
            if (htbctype == SUB_HT_BC_Cauchy_CONST || htbctype == SUB_HT_BC_Cauchy_T || htbctype == SUB_HT_BC_FD)   {
                Kokkos::parallel_for("ht_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
                    int iGlob = bcells[ibc], iGhost = gcells[ibc], htiGhost = htgwcells[ibc];
					ht.T(iGhost,1) = ht.T(iGlob,1);
                });
            }
		//   20240617
			//一类边界 CONST or 三类边界CONST
          else if (htbctype == SUB_HT_BC_Dirichlet_CONST ){
				Kokkos::parallel_for("ht_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
                    int iGlob = bcells[ibc], iGhost = gcells[ibc], htiGhost = htgwcells[ibc]; //20240510
                         ht.T(htiGhost,1) = bcvals(ibc);
			// std::cout << "rt.c(iGhost,1) = " << rt.c(rtiGhost,1) << std::endl;	
                });
			}
            // H Time series
            else if (htbctype == SUB_HT_BC_Dirichlet_T)  {
                Kokkos::parallel_for("ht_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
                    int iGlob = bcells[ibc], iGhost = gcells[ibc], htiGhost = htgwcells[ibc]; //20240510
                    ht.T(htiGhost,1) = Tembc;
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
    inline void applyHTMatBC(HTState &ht, GwState &gw, GwDomain &gdom, Parallel &par) {
		// Check if on global boundaries
		bool onBoundary = 0;
		if (direction == 2 && par.px == 0)	{onBoundary = 1;}
		else if (direction == 1 && par.px == par.nproc_x-1)	{onBoundary = 1;}
		else if (direction == 4 && par.py == 0)	{onBoundary = 1;}
		else if (direction == 3 && par.py == par.nproc_y-1)	{onBoundary = 1;}
		else if (direction == 5 || direction == 6)	{onBoundary = 1;}

        if (ncellsBC > 0 && onBoundary == 1) {
            real Cauchybc;
            if (htbctype == SUB_HT_BC_Cauchy_T) {Cauchybc = interpolateLinear(ts, gdom.etime);}
	  	    switch (htbctype) {
                case SUB_HT_BC_Dirichlet_CONST:    
               //  case SUB_RT_BC_Cauchy_CONST:   
                case SUB_HT_BC_Dirichlet_T:    
               //  case SUB_RT_BC_Cauchy_T:
                    Kokkos::parallel_for("ht_bc_Dirichlet_const", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
                        int ii, jj, kk, idom, iGlobSW, iGlob = bcells[ibc], iGhost = gcells[ibc], htiGhost = htgwcells[ibc]; //20240510
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
						iGlobSW = jj*gdom.nxhc + ii;
						if (direction == 1)	{
							// rt.RTcoef(idom,1) = rt.RTcoef(idom,1) * 2.0;
							// rt.RTcoef(idom,7) -= rt.RTcoef(idom,1) * rt.c(iGlob+1,1);
							// rt.RTcoef(idom,1) = 0;

							ht.HTcoef(idom,7) = 1e50 * ht.HTcoef(idom,0) * ht.T(iGlob+1,1);
							ht.HTcoef(idom,0) = 1e50 * ht.HTcoef(idom,0);
						}
						else if (direction == 2)	{
							// rt.RTcoef(idom,2) = rt.RTcoef(idom,2) * 2.0;
							// rt.RTcoef(idom,7) -= rt.RTcoef(idom,2) * rt.c(iGlob-1,1);
							// rt.RTcoef(idom,2) = 0;

							ht.HTcoef(idom,7) = 1e50 * ht.HTcoef(idom,0) * ht.T(iGlob-1,1);
							ht.HTcoef(idom,0) = 1e50 * ht.HTcoef(idom,0);
						}
						else if (direction == 3)	{
							// rt.RTcoef(idom,3) = rt.RTcoef(idom,3) * 2.0;
							// rt.RTcoef(idom,7) -= rt.RTcoef(idom,3) * rt.c(iGlob-gdom.nxhc,1);
							// rt.RTcoef(idom,3) = 0;

							ht.HTcoef(idom,7) = 1e50 * ht.HTcoef(idom,0) * ht.T(iGlob-gdom.nxhc,1);
							ht.HTcoef(idom,0) = 1e50 * ht.HTcoef(idom,0);							
						}
						else if (direction == 4)	{
							// rt.RTcoef(idom,4) = rt.RTcoef(idom,4) * 2.0;
							// rt.RTcoef(idom,7) -= rt.RTcoef(idom,4) * rt.c(iGlob-gdom.nxhc,1);
							// rt.RTcoef(idom,4) = 0;

							ht.HTcoef(idom,7) = 1e50 * ht.HTcoef(idom,07) * ht.T(iGlob-gdom.nxhc,1);
							ht.HTcoef(idom,0) = 1e50 * ht.HTcoef(idom,0);
						}
						else if (direction == 5)	{
							// rt.RTcoef(idom,5) = rt.RTcoef(idom,5) * 2.0;
							// rt.RTcoef(idom,7) -= rt.RTcoef(idom,5) * rt.c(iGlob+gdom.nxhc*gdom.nyhc,1);
							// rt.RTcoef(idom,5) = 0;

							ht.HTcoef(idom,7) = 1e50 * ht.HTcoef(idom,0) * ht.T(iGlob+gdom.nxhc*gdom.nyhc,1);
							ht.HTcoef(idom,0) = 1e50 * ht.HTcoef(idom,0);
						}
						else if (direction == 6)	{
							// rt.RTcoef(idom,6) = rt.RTcoef(idom,6) * 2.0;
							// rt.RTcoef(idom,7) -= rt.RTcoef(idom,6) * rt.c(iGlob-gdom.nxhc*gdom.nyhc,1);
							// rt.RTcoef(idom,6) = 0;

							ht.HTcoef(idom,7) = 1e50 * ht.HTcoef(idom,0) * ht.T(iGlob-gdom.nxhc*gdom.nyhc,1);
							ht.HTcoef(idom,0) = 1e50 * ht.HTcoef(idom,0);
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
				case SUB_HT_BC_Cauchy_CONST:
					Kokkos::parallel_for("ht_bc_Cauchy_const", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
						int ii, jj, kk, idom, iGlobSW, iGlob = bcells[ibc];
						gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
						idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
						iGlobSW = jj*gdom.nxhc + ii;
						if (direction == 1)	{
							ht.HTcoef(idom,7) += gdom.dt * bcvals(ibc) / gdom.dx;
							ht.HTcoef(idom,1) = 0.0;
						}
						else if (direction == 2)	{
							ht.HTcoef(idom,7) -= gdom.dt * bcvals(ibc) / gdom.dx ;
							ht.HTcoef(idom,2) = 0.0;
						}
						else if (direction == 3)	{
							ht.HTcoef(idom,7) += gdom.dt * bcvals(ibc) / gdom.dx ;
							ht.HTcoef(idom,3) = 0.0;
						}
						else if (direction == 4)	{
							ht.HTcoef(idom,7) -= gdom.dt * bcvals(ibc) / gdom.dx;
							ht.HTcoef(idom,4) = 0.0;
						}
						else if (direction == 5)	{
							ht.HTcoef(idom,5) = 0.0;
						}
						else if (direction == 6)	{
							// rt.RTcoef(idom,7) -= gdom.dt * gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) ;
							ht.HTcoef(idom,7) -= gdom.dt * bcvals(ibc) / gdom.dz(iGlob);//20240617上边界不理解？
							// if (gdom.isRain)    {rt.RTcoef(idom,7) += gdom.dt * gdom.rainRate(iGlobSW) / gdom.dz(iGlob);}
							// if (gdom.isEvap)	{rt.RTcoef(idom,7) -= gdom.dt * gdom.evapRate(iGlobSW) / gdom.dz(iGlob);}
							ht.HTcoef(idom,6) = 0.0;
						}
					});
					break;
                case SUB_HT_BC_Cauchy_T:   
                    Kokkos::parallel_for("ht_bc_Cauchy_t", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
					// for (int ibc = 0; ibc < ncellsBC; ibc++)	{
                        int ii, jj, kk, idom, iGlobSW, iGlob = bcells[ibc];
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
						iGlobSW = jj*gdom.nxhc + ii;
						if (direction == 1)	{
							ht.HTcoef(idom,7) += gdom.dt * Cauchybc / gdom.dx;
							ht.HTcoef(idom,1) = 0.0;
						}
						else if (direction == 2)	{
							ht.HTcoef(idom,7) -= gdom.dt * Cauchybc / gdom.dx;
							ht.HTcoef(idom,2) = 0.0;
						}
						else if (direction == 3)	{
							ht.HTcoef(idom,7) += gdom.dt * Cauchybc / gdom.dx;
							ht.HTcoef(idom,3) = 0.0;
						}
						else if (direction == 4)	{
							ht.HTcoef(idom,7) -= gdom.dt * Cauchybc / gdom.dx;
							ht.HTcoef(idom,4) = 0.0;
						}
						else if (direction == 5)	{
							//20240617
							ht.HTcoef(idom,7) -= gdom.dt * Cauchybc / gdom.dz(iGlob);
							ht.HTcoef(idom,5) = 0;

						}
						else if (direction == 6)	{
							// rt.RTcoef(idom,7) -= gdom.dt * gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) / gdom.dz(iGlob);//20240217注释掉
							ht.HTcoef(idom,7) -= gdom.dt * Cauchybc / gdom.dz(iGlob);
							// if (gdom.isRain)    {rt.RTcoef(idom,7) += gdom.dt * gdom.rainRate(iGlobSW) / gdom.dz(iGlob);}
							// if (gdom.isEvap)	{rt.RTcoef(idom,7) -= gdom.dt * gdom.evapRate(iGlobSW) / gdom.dz(iGlob);}
							ht.HTcoef(idom,6) = 0.0;
						}
                    });
					// }
                    break;
				/*case SUB_RT_BC_FD:
					Kokkos::parallel_for("rt_bc_fd", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
					// for (int ibc = 0; ibc < ncellsBC; ibc++)	{
				    int ii, jj, kk, idom, iGlobSW, iGlob = bcells[ibc];
				    gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
				    idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
				    		iGlobSW = jj*gdom.nxhc + ii;
					double Up_Weighting_x,Up_Weighting_y,Up_Weighting_z;
							if (rt.aveV(iGlob, 0) > 0)
							{
								Up_Weighting_x = rt.Up_Weighting_vplus;
							}
							else
							{
								Up_Weighting_x = rt.Up_Weighting_vminus;
							}

					//y方向速度方向判断
							if (rt.aveV(iGlob, 1) > 0)
							{
								Up_Weighting_y = rt.Up_Weighting_vplus;
							}
							else
							{
								Up_Weighting_y = rt.Up_Weighting_vminus;
							}
					//z方向速度方向判断
							if (rt.aveV(iGlob, 2) > 0)
							{
								Up_Weighting_z = rt.Up_Weighting_vplus;
							}
							else
							{
								Up_Weighting_z = rt.Up_Weighting_vminus;
							}
						

						if (direction == 1) {


							// rt.RTcoef(idom, 0) = 1.0 * gw.wc(iGlob,0) 
							// + rt.c_difxx(iGlob)+Up_Weighting_x*rt.c_advxx(iGlob-1)-(1-Up_Weighting_x)*rt.c_advxx(iGlob)
							// + rt.c_difyy(iGlob)+rt.c_difyy(iGlob-gdom.nxhc)+Up_Weighting_y*rt.c_advyy(iGlob-gdom.nxhc)-(1-Up_Weighting_y)*rt.c_advyy(iGlob)
							// + rt.c_difzz(iGlob)+rt.c_difzz(iGlob-gdom.nxhc*gdom.nyhc)+Up_Weighting_z*rt.c_advzz(iGlob-gdom.nxhc*gdom.nyhc)-(1-Up_Weighting_z)*rt.c_advzz(iGlob);
							// rt.RTcoef(idom, 1) = - (rt.c_difxx(iGlob)  - (1 - Up_Weighting_x) * rt.c_advxx(iGlob)) ;//Ci+1
							// rt.RTcoef(idom, 2) = - (Up_Weighting_x * rt.c_advxx(iGlob - 1) ); // # Ci-1							
							// rt.RTcoef(idom, 7) = rt.c(iGlob, 1) * gw.wc(iGlob,0);


							// rt.RTcoef(idom, 0) = 1.0 * gw.wc(iGlob,0) 
							// +rt.c_difxx(iGlob-1)+Up_Weighting_x*rt.c_advxx(iGlob-1)-(1-Up_Weighting_x)*rt.c_advxx(iGlob)
							// + rt.c_difyy(iGlob)+rt.c_difyy(iGlob-gdom.nxhc)+Up_Weighting_y*rt.c_advyy(iGlob-gdom.nxhc)-(1-Up_Weighting_y)*rt.c_advyy(iGlob)
							// + rt.c_difzz(iGlob)+rt.c_difzz(iGlob-gdom.nxhc*gdom.nyhc)+Up_Weighting_z*rt.c_advzz(iGlob-gdom.nxhc*gdom.nyhc)-(1-Up_Weighting_z)*rt.c_advzz(iGlob);		
							// rt.RTcoef(idom, 1) = - (  - (1 - Up_Weighting_x) * rt.c_advxx(iGlob)) ;//Ci+1
							// rt.RTcoef(idom, 2) = - (Up_Weighting_x * rt.c_advxx(iGlob - 1) + rt.c_difxx(iGlob -1 )); // # Ci-1
							// rt.RTcoef(idom, 7) = rt.c(iGlob, 1) * gw.wc(iGlob,0);

							rt.RTcoef(idom, 0) = 1.0;		
							rt.RTcoef(idom, 1) = 0;//Ci+1
							rt.RTcoef(idom, 2) = -1; // # Ci-1
							rt.RTcoef(idom, 3) = 0;
							rt.RTcoef(idom, 4) = 0;
							rt.RTcoef(idom, 5) = 0;
							rt.RTcoef(idom, 6) = 0;
							rt.RTcoef(idom, 7) = 0;

						}
						else if (direction == 2) {
							rt.RTcoef(idom, 0) = 1.0;		
							rt.RTcoef(idom, 1) = -1;//Ci+1
							rt.RTcoef(idom, 2) = 0; // # Ci-1
							rt.RTcoef(idom, 3) = 0;//Cj+1
							rt.RTcoef(idom, 4) = 0;//Cj-1
							rt.RTcoef(idom, 5) = 0;// Ck+1 
							rt.RTcoef(idom, 6) = 0;// Ck-1 
							rt.RTcoef(idom, 7) = 0;							
						}
						else if (direction == 3) {
							rt.RTcoef(idom, 3) = 0.0;
						}
						else if (direction == 4) {
							// rt.RTcoef(idom, 0) = 1.0;		
							// rt.RTcoef(idom, 1) = 0;//Ci+1
							// rt.RTcoef(idom, 2) = 0; // # Ci-1
							// rt.RTcoef(idom, 3) = 0;//Cj+1
							// rt.RTcoef(idom, 4) = 0;//Cj-1
							// rt.RTcoef(idom, 5) = 0;// Ck+1 
							// rt.RTcoef(idom, 6) = -1;// Ck-1 
							// rt.RTcoef(idom, 7) = 0;
						}
						else if (direction == 5) {
							rt.RTcoef(idom, 0) = 1.0;		
							rt.RTcoef(idom, 1) = 0;//Ci+1
							rt.RTcoef(idom, 2) = 0; // # Ci-1
							rt.RTcoef(idom, 3) = 0;//Cj+1
							rt.RTcoef(idom, 4) = 0;//Cj-1
							rt.RTcoef(idom, 5) = 0;// Ck+1 
							rt.RTcoef(idom, 6) = -1;// Ck-1 
							rt.RTcoef(idom, 7) = 0;
						}
						else if (direction == 6) {
							rt.RTcoef(idom, 6) = 0.0;
						}
					});
					break;*/
            }
        }
    }
};

class HTSubsurfaceBoundaries{
// This class should not be invoked form a parallel region as it contains strings
public:
  	std::string BoundaryTypes[9] = {"NOFLOW","CONST_H","CONST_Q","CONST_WT","H_TIMESERIES","Q_TIMESERIES","WT_TIMESERIES","SWEXCHANGE","FREE_DRAINAGE"};
	std::vector<std::string> id;
	std::vector<HTBC> htgwbc;
};
#endif

#endif

