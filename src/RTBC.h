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
	real QMassTot, QMassInflow, QMassOutflow;

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
				std::cerr << RERROR << "No boundary cells found for subsurface reaction_transport boundary with id '" << rtid << "'" << std::endl;
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
            if (rtbctype == SUB_RT_BC_Dirichlet_T ) 
			{Conbc = interpolateLinear(ts, gdom.etime);}
            // zero gradient if Q BC is specified
		  //20240617
            if (rtbctype == SUB_RT_BC_Cauchy_CONST || rtbctype == SUB_RT_BC_Cauchy_T || rtbctype == SUB_RT_BC_FD)   {
                Kokkos::parallel_for("rt_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
                    int iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc];
					rt.c(rtiGhost,1) = rt.c(iGlob,1);
                });
            }
		//   20240617
			//一类边界 CONST 
          else if (rtbctype == SUB_RT_BC_Dirichlet_CONST ){
				Kokkos::parallel_for("rt_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
                    int iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc]; //20240510
                    rt.c(rtiGhost,1) = bcvals(ibc);
                });
			}
            // Dirichlet Time series
            else if (rtbctype == SUB_RT_BC_Dirichlet_T)  {
                Kokkos::parallel_for("rt_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
                    int iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc]; //20240510
                    rt.c(rtiGhost,1) = Conbc;
                });
            }

            // Surface-subsurface exchange for reaction-transport
            else if (rtbctype == SUB_RT_BC_SWE)  {

				
				// #if SERGHEI_SURFACE_TRANSPORT
				if (direction == 6)	{
					Kokkos::parallel_for("rt_bc_c", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
	                    int iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc];
	                    // printf("ncellsBC=%d, ibc=%d, iGlob=%d, iGhost=%d, rtiGhost=%d\n",ncellsBC, ibc, iGlob, iGhost, rtiGhost);
						int ii, jj, kk, iGlobSW;
						gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
						iGlobSW = jj*gdom.nxhc + ii;//with halo
						int iGlobSW2 = (jj-1)*gdom.nx + ii - 1;//without halo

						
						// !Get exchange flux
						// gw.qss(iGlobSW2)//地表地下水交换速率

						// printf("55555gw.qss(iGlobSW)=%e\n",gw.qss(iGlobSW));
						if (gw.qss(iGlobSW2) <= 0){//SW-GW流速向下（入渗）// 地表水→地下水
							//地下网格上边界浓度为地表网格浓度值
							rt.c(rtiGhost, 1) = rt.csw(iGlobSW);
							// rt.c(rtiGhost,1) = rt.csw(ibc);
							//! 如果存在降雨情景，当地表水水深小于hmin时，地表水浓度为0
							//! 此时假设降雨全部入渗地下水，即当存在降雨情景且地表水浓度为0时，入渗边界浓度直接取雨水浓度而不是读取地表水浓度
							if (rt.csw(iGlobSW) == 0 && gdom.isRain == 1 && gw.qss(iGlobSW2) != 0){
							rt.c(rtiGhost, 1) = 1e-15;//雨水浓度
							}
							// if (iGlobSW2 == 99){
							// printf("1111rtiGhost=%d, rt.c(rtiGhost, 1)= %e, rt.c(iGhost, 1)=%e, rt.csw(iGlobSW)=%e\n", rtiGhost, rt.c(rtiGhost, 1),rt.c(iGhost, 1), rt.csw(iGlobSW));
							// }
						}
						else {//地下水→地表水,地下网格上边界变为溶质流出三类边界，故边界浓度为地下网格上边界网格浓度
							rt.c(rtiGhost, 1) = rt.c(iGlob, 1);
						// 	if (iGlobSW2 == 99){
							// printf("!!!!2222iGlobSW2=%d, iGlob=%d, rtiGhost=%d, iGhost=%d, rt.c(rtiGhost, 1)= %e, rt.c(iGhost, 0)= %e, rt.c(iGhost, 1)=%e, rt.c(iGlob, 1)=%e\n", 
							// 	iGlobSW2, iGlob, rtiGhost, iGhost, rt.c(rtiGhost, 1), rt.c(iGhost, 0), rt.c(iGhost, 1), rt.c(iGlob, 1));
						// }
						}

						//计算地表地下溶质交换通量q*c
						rt.ConQss(iGlobSW2) = gw.qss(iGlobSW2) * rt.c(rtiGhost, 1);
						// if (ii == 10){
						// if (iGlobSW2 == 99){
						// // if (rt.c(iGhost, 1) < 0){
						// printf("rtiGhost=%d, iGlobSW2=%d, rt.ConQss(iGlobSW2)=%e, gw.qss(iGlobSW2)=%e,rt.c(rtiGhost,1)=%f,rt.csw(iGlobSW)=%f\n",rtiGhost, iGlobSW2, rt.ConQss(iGlobSW2) ,gw.qss(iGlobSW2),rt.c(rtiGhost,1),rt.csw(iGlobSW));
						// // printf("2222rtiGhost=%d, rt.c(rtiGhost, 1)= %e, rt.c(iGhost, 0)= %e, rt.c(iGhost, 1)=%e\n", rtiGhost, rt.c(rtiGhost, 1), rt.c(iGhost, 0), rt.c(iGhost, 1));

						// }
						// }
					
					//! 添加判断条件，如果地下水中溶质浓度小于0，则ConQss设为0

	                });
				}
				else {
					if (par.masterproc)	{std::cerr << RERROR "RTBC direction must be 6 for RT-SW-GW exchange boundary! " << "\n";}
				}
				// #endif
            }

		    /* --------------------------------------------------
		calculate_mass_balance: 计算边界质量通量
	-------------------------------------------------- */

			if (direction == 1) //右边界
			{
				Kokkos::parallel_reduce("reducexRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
					int iGlob = bcells[ibc];
					tmp += rt.aveVB(iGlob, 0) * rt.c(iGlob, 1) * gdom.dy * gdom.dz(iGlob); }, Kokkos::Sum<real>(QMassTot));
				if (QMassTot > 0)
				{
					QMassOutflow = QMassTot;
				}
				else
				{
					QMassInflow = -QMassTot;
				}
			}
			else if (direction == 2) //左边界
			{
				Kokkos::parallel_reduce("reducexRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
					int iGlob = bcells[ibc];
					tmp += rt.aveVB(iGlob, 0) * rt.c(iGlob, 1) * gdom.dy * gdom.dz(iGlob); }, Kokkos::Sum<real>(QMassTot));
				if (QMassTot > 0)
				{
					QMassInflow = QMassTot;
				}
				else
				{
					QMassOutflow = -QMassTot;
				}
			}
			// --- Y 方向边界 ---
	
			if (direction == 3) {
				Kokkos::parallel_reduce("reduceyRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
					int iGlob = bcells[ibc];
					tmp += rt.aveVB(iGlob, 1) * rt.c(iGlob, 1) * gdom.dx * gdom.dz(iGlob); }, Kokkos::Sum<real>(QMassTot));
				if (QMassTot > 0)
				{
					QMassOutflow = QMassTot;
				}
				else
				{
					QMassInflow = -QMassTot;
				}
			}
			else if (direction == 4) {	
				Kokkos::parallel_reduce("reduceyRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
					int iGlob = bcells[ibc];
					tmp += rt.aveVB(iGlob, 1) * rt.c(iGlob, 1) * gdom.dx * gdom.dz(iGlob); }, Kokkos::Sum<real>(QMassTot));
				if (QMassTot > 0)
				{
					QMassInflow = QMassTot;
				}
				else
				{
					QMassOutflow = -QMassTot;
				}
			}

			else if (direction == 5) //下边界
			{
		
				Kokkos::parallel_reduce("reducezRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
					int iGlob = bcells[ibc];
					tmp += rt.aveVB(iGlob, 2) * rt.c(iGlob, 1) * gdom.dx * gdom.dy; }, Kokkos::Sum<real>(QMassTot));
				if (QMassTot > 0)
				{
					QMassOutflow = QMassTot;
				}
				else
				{
					QMassInflow = -QMassTot;
				}
			}
			else if (direction == 6) //上边界
			{
				
				Kokkos::parallel_reduce("reducezRT", ncellsBC, KOKKOS_CLASS_LAMBDA(int ibc, real &tmp) {
					int iGlob = bcells[ibc];
					tmp += rt.aveVB(iGlob, 2) * rt.c(iGlob, 1) * gdom.dx * gdom.dy; }, Kokkos::Sum<real>(QMassTot));
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
            real Cauchybc;
            if (rtbctype == SUB_RT_BC_Cauchy_T) {Cauchybc = interpolateLinear(ts, gdom.etime);}
	  	    switch (rtbctype) {
                case SUB_RT_BC_Dirichlet_CONST:    
               //  case SUB_RT_BC_Cauchy_CONST:   
                case SUB_RT_BC_Dirichlet_T:    
               //  case SUB_RT_BC_Cauchy_T:
                    Kokkos::parallel_for("rt_bc_Dirichlet_const", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
                        int ii, jj, kk, idom,  iGlob = bcells[ibc], iGhost = gcells[ibc], rtiGhost = rtgwcells[ibc]; //20240510
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
						
						if (direction == 1)	{
							// rt.RTcoef(idom,1) = rt.RTcoef(idom,1) * 2.0;
							// rt.RTcoef(idom,7) -= rt.RTcoef(idom,1) * rt.c(iGlob+1,1);
							// rt.RTcoef(idom,1) = 0;

							rt.RTcoef(idom,7) = 1e50 * rt.RTcoef(idom,0) * rt.c(iGlob+1,1);
							rt.RTcoef(idom,0) = 1e50 * rt.RTcoef(idom,0);
						}
						else if (direction == 2)	{
							// rt.RTcoef(idom,7) = 1e50 * rt.RTcoef(idom,0) * rt.c(iGlob-1,1);
							// rt.RTcoef(idom,0) = 1e50 * rt.RTcoef(idom,0);	

							// if (kk<40){
							if (gw.h(iGlob, 1) >= 0) {// 含水层所在位置
								//! zzb 变水位边界条件时，随水位调节定浓度边界位置

							// 	// if (kk < 40) {
							if (gw.q(iGlob, 0) <= 0) {//! zzb 溶质流入含水层
								rt.RTcoef(idom,7) = 1e50 * rt.RTcoef(idom,0) * rt.c(iGlob-1,1);
								rt.RTcoef(idom,0) = 1e50 * rt.RTcoef(idom,0);						
								}
							else {//! zzb 含水层水位高于河流水位，溶质自由流出
								for (int i = 0; i < 30; i++) {
									rt.RTcoef(idom, i) = 0.0;
								}
								rt.RTcoef(idom, 0) = 1.0;		
								rt.RTcoef(idom, 1) = -1;//Ci+1
							}

							}
							
							// }
							else {//! zzb 不饱和带，渗流面溶质自由流出
							for (int i = 0; i < 30; i++) {
								rt.RTcoef(idom, i) = 0.0;
							}
							rt.RTcoef(idom, 0) = 1.0;		
							rt.RTcoef(idom, 1) = -1;//Ci+1
							}
							// }
						}
						else if (direction == 3)	{
							// rt.RTcoef(idom,3) = rt.RTcoef(idom,3) * 2.0;
							// rt.RTcoef(idom,7) -= rt.RTcoef(idom,3) * rt.c(iGlob-gdom.nxhc,1);
							// rt.RTcoef(idom,3) = 0;

							rt.RTcoef(idom,7) = 1e50 * rt.RTcoef(idom,0) * rt.c(iGlob+gdom.nxhc,1);
							rt.RTcoef(idom,0) = 1e50 * rt.RTcoef(idom,0);
						}
						else if (direction == 4)	{
							// rt.RTcoef(idom,4) = rt.RTcoef(idom,4) * 2.0;
							// rt.RTcoef(idom,7) -= rt.RTcoef(idom,4) * rt.c(iGlob-gdom.nxhc,1);
							// rt.RTcoef(idom,4) = 0;

							rt.RTcoef(idom,7) = 1e50 * rt.RTcoef(idom,0) * rt.c(iGlob-gdom.nxhc,1);
							rt.RTcoef(idom,0) = 1e50 * rt.RTcoef(idom,0);
						}
						else if (direction == 5)	{
							// rt.RTcoef(idom,5) = rt.RTcoef(idom,5) * 2.0;
							// rt.RTcoef(idom,7) -= rt.RTcoef(idom,5) * rt.c(iGlob+gdom.nxhc*gdom.nyhc,1);
							// rt.RTcoef(idom,5) = 0;

							rt.RTcoef(idom,7) = 1e50 * rt.RTcoef(idom,0) * rt.c(iGlob+gdom.nxhc*gdom.nyhc,1);
							rt.RTcoef(idom,0) = 1e50 * rt.RTcoef(idom,0);
						}
						else if (direction == 6)	{
							// rt.RTcoef(idom,6) = rt.RTcoef(idom,6) * 2.0;
							// rt.RTcoef(idom,7) -= rt.RTcoef(idom,6) * rt.c(iGlob-gdom.nxhc*gdom.nyhc,1);
							// rt.RTcoef(idom,6) = 0;

							rt.RTcoef(idom,7) = 1e50 * rt.RTcoef(idom,0) * rt.c(iGlob-gdom.nxhc*gdom.nyhc,1);
							rt.RTcoef(idom,0) = 1e50 * rt.RTcoef(idom,0);
							
							// if (gdom.isEvap)	{
							// 	rt.RTcoef(idom,7) -= gdom.dt * gdom.evapRate(iGlobSW) / gdom.dz(iGlob);
							// }
						}
                    });
                    break;
                case SUB_RT_BC_SWE:

				
				
                	#if SERGHEI_SURFACE_TRANSPORT
					
					if (direction == 6)	{

						
						
						Kokkos::parallel_for("rtgw_bc_swe_const", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
	                        int ii, jj, kk, idom, iGlobSW, iGlob = bcells[ibc], rtiGhost = rtgwcells[ibc];
	                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
	                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
							iGlobSW = (jj-1)*gdom.nx + ii - 1;//with halo
							int iGlobSW2 = (jj-1)*gdom.nx + ii - 1;//without halo
					

							// if (gw.qss(iGlobSW) <= 0) {//SW-GW流速向下（入渗）


							
							//sw-gw for rt 当做源汇项处理
							if (iGlobSW2 < 100) {

							rt.RTcoef(idom, 7) += gdom.dt  * (-1*rt.ConQss(iGlobSW2))  / gdom.dz(iGlob);
							}
							
							// rt.RTcoef(idom, 0) +=  gdom.dt * (-1*gw.qss(iGlobSW2))  / gdom.dz(iGlob); //主对角系数
							// rt.RTcoef(idom, 6) = 0.0;


// if (iGlobSW2 == 99){
// 							printf("333rt.ConQss(iGlobSW2)=%e,gdom.dt  * (-1*rt.ConQss(iGlobSW2))  / gdom.dz(iGlob)=%e\n",rt.ConQss(iGlobSW2), gdom.dt  * (-1*rt.ConQss(iGlobSW2))  / gdom.dz(iGlob));
// }							
							
							// }

							// else {
							// rt.RTcoef(idom,7) -= gdom.dt * rt.c(rtiGhost,1) * (-1 * gw.qss(iGlobSW))  / gdom.dz(iGlob);//20240617上边界不理解？
							
							// if (std::isnan(rt.RTcoef(idom,7))) {rt.RTcoef(idom,7) = 0.0;}
							// rt.RTcoef(idom,6) = 0.0;
							// printf("222rt.RTcoef(idom,7)=%e,rt.c(rtiGhost,1)=%e, gw.qss(iGlobSW)=%e\n",rt.RTcoef(idom,7), rt.c(rtiGhost,1), gw.qss(iGlobSW));
							// }



							
	                    });
					}
                    else {
						if (par.masterproc)	{std::cerr << RERROR "RTGWBC direction must be 6 for rt-sw-gw boundary! " << "\n";}
					}
                    #endif
                    break;
				case SUB_RT_BC_Cauchy_CONST:
					Kokkos::parallel_for("rt_bc_Cauchy_const", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
						int ii, jj, kk, idom,  iGlob = bcells[ibc], rtiGhost = rtgwcells[ibc];
						gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
						idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
						if (direction == 1)	{
							rt.RTcoef(idom,7) += gdom.dt * bcvals(ibc)  * 1000 / gdom.dx;
							rt.RTcoef(idom,1) = 0.0;
						}
						else if (direction == 2)	{

							rt.RTcoef(idom,7) -= gdom.dt * bcvals(ibc)  * 1000 / gdom.dx ;
							rt.RTcoef(idom,2) = 0.0;
							
						}
						else if (direction == 3)	{
							rt.RTcoef(idom,7) += gdom.dt * bcvals(ibc)  * 1000 / gdom.dx ;
							rt.RTcoef(idom,3) = 0.0;
						}
						else if (direction == 4)	{
							rt.RTcoef(idom,7) -= gdom.dt * bcvals(ibc)  * 1000 / gdom.dx;
							rt.RTcoef(idom,4) = 0.0;
						}
						else if (direction == 5)	{
							rt.RTcoef(idom,7) += gdom.dt * bcvals(ibc) * 1000 / gdom.dz(iGlob);
							rt.RTcoef(idom,5) = 0.0;
						}
						else if (direction == 6)	{
							
							//todo 需要结合流速方向判断溶质流入还是流出
							//第三类边界本质上是混合边界，需要同时修改主对角系数和相邻单元系数
							//当溶质边界为第三类边界时，input中bcval值代表浓度大小，流量取自gwbc.input文件中bcval值，两者相乘为溶质通量
							//rt.RTcoef(idom, 7)为右端项，量纲为mg/m^3					
							rt.RTcoef(idom, 7) += gw.wc(iGlob, 1) * gdom.dt  * (-1*gw.q(iGlob, 2)) * bcvals(ibc) / gdom.dz(iGlob);
							//rt.RTcoef(idom, 0)为主对角系数项，量纲为[-]
							//todo 忽略了弥散项贡献，theta*Dn/(dz)^2
							rt.RTcoef(idom, 0) += gw.wc(iGlob, 1) * gdom.dt * (-1*gw.q(iGlob, 2))  / gdom.dz(iGlob); //主对角系数
							rt.RTcoef(idom, 6) = 0.0;							
							// printf("gw.q(iGlob, 2)=%e,bcvals(ibc)=%e\n",gw.q(iGlob, 2),bcvals(ibc));

//流入情景
							// if (gw.q(iGlob, 2) <= 0){
							// 	rt.RTcoef(idom, 7) += gw.wc(iGlob, 1) * gdom.dt  * (-1*gw.q(iGlob, 2)) * bcvals(ibc) / gdom.dz(iGlob);
							// 	//rt.RTcoef(idom, 0)为主对角系数项，量纲为[-]
							// 	//todo 忽略了弥散项贡献，theta*Dn/(dz)^2
							// 	rt.RTcoef(idom, 0) += 0.0; //主对角系数
							// 	rt.RTcoef(idom, 6) = 0.0;							
	
							// }
							// //流出情景
							// else {
							// 	rt.RTcoef(idom, 7) += 0;
							// 	//rt.RTcoef(idom, 0)为主对角系数项，量纲为[-]
							// 	//todo 忽略了弥散项贡献，theta*Dn/(dz)^2
							// 	rt.RTcoef(idom, 0) += gw.wc(iGlob, 1) * gdom.dt * (-1*gw.q(iGlob, 2))  / gdom.dz(iGlob); //主对角系数
							// 	rt.RTcoef(idom, 6) = 0.0;							
	
							// }

						}
					});
					break;
                case SUB_RT_BC_Cauchy_T:   
                    Kokkos::parallel_for("rt_bc_Cauchy_t", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
					// for (int ibc = 0; ibc < ncellsBC; ibc++)	{
                        int ii, jj, kk, idom, iGlobSW, iGlob = bcells[ibc];
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
						iGlobSW = jj*gdom.nxhc + ii;
						if (direction == 1)	{
							rt.RTcoef(idom,7) += gdom.dt * Cauchybc  * 1000 / gdom.dx;
							rt.RTcoef(idom,1) = 0.0;
						}
						else if (direction == 2)	{
							rt.RTcoef(idom,7) -= gdom.dt * Cauchybc  * 1000 / gdom.dx;
							rt.RTcoef(idom,2) = 0.0;
						}
						else if (direction == 3)	{
							rt.RTcoef(idom,7) += gdom.dt * Cauchybc  * 1000 / gdom.dx;
							rt.RTcoef(idom,3) = 0.0;
						}
						else if (direction == 4)	{
							rt.RTcoef(idom,7) -= gdom.dt * Cauchybc  * 1000 / gdom.dx;
							rt.RTcoef(idom,4) = 0.0;
						}
						else if (direction == 5)	{
							//20240617
							rt.RTcoef(idom,7) -= gdom.dt * Cauchybc  * 1000 / gdom.dz(iGlob);
							rt.RTcoef(idom,5) = 0;

						}
						else if (direction == 6)	{
							// rt.RTcoef(idom,7) -= gdom.dt * gw.k(iGlob-gdom.nxhc*gdom.nyhc,2) / gdom.dz(iGlob);//20240217注释掉
							rt.RTcoef(idom,7) -= gdom.dt * Cauchybc  * 1000 / gdom.dz(iGlob);
							// if (gdom.isRain)    {rt.RTcoef(idom,7) += gdom.dt * gdom.rainRate(iGlobSW) / gdom.dz(iGlob);}
							// if (gdom.isEvap)	{rt.RTcoef(idom,7) -= gdom.dt * gdom.evapRate(iGlobSW) / gdom.dz(iGlob);}
							rt.RTcoef(idom,6) = 0.0;
						}
                    });
					// }
                    break;
				case SUB_RT_BC_FD:
					Kokkos::parallel_for("rt_bc_fd", ncellsBC, KOKKOS_CLASS_LAMBDA (int ibc){
					// for (int ibc = 0; ibc < ncellsBC; ibc++)	{
				    int ii, jj, kk, idom, iGlobSW, iGlob = bcells[ibc];
				    gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
				    idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
				    		iGlobSW = jj*gdom.nxhc + ii;
					double Up_Weighting_x,Up_Weighting_y,Up_Weighting_z;
							if (rt.aveVB(iGlob, 0) > 0)
							{
								Up_Weighting_x = rt.Up_Weighting_vplus;
							}
							else
							{
								Up_Weighting_x = rt.Up_Weighting_vminus;
							}

					//y方向速度方向判断
							if (rt.aveVB(iGlob, 1) > 0)
							{
								Up_Weighting_y = rt.Up_Weighting_vplus;
							}
							else
							{
								Up_Weighting_y = rt.Up_Weighting_vminus;
							}
					//z方向速度方向判断
							if (rt.aveVB(iGlob, 2) > 0)
							{
								Up_Weighting_z = rt.Up_Weighting_vplus;
							}
							else
							{
								Up_Weighting_z = rt.Up_Weighting_vminus;
							}
						

						if (direction == 1) {
							for (int i = 0; i < 30; i++) {
								rt.RTcoef(idom, i) = 0.0;
							}							
							rt.RTcoef(idom, 0) = 1.0;		
							// rt.RTcoef(idom, 1) = 0;//Ci+1
							rt.RTcoef(idom, 2) = -1; // # Ci-1
							// rt.RTcoef(idom, 3) = 0;
							// rt.RTcoef(idom, 4) = 0;
							// rt.RTcoef(idom, 5) = 0;
							// rt.RTcoef(idom, 6) = 0;
							// rt.RTcoef(idom, 7) = 0;

						}
						else if (direction == 2) {
							for (int i = 0; i < 30; i++) {
								rt.RTcoef(idom, i) = 0.0;
							}
							rt.RTcoef(idom, 0) = 1.0;		
							rt.RTcoef(idom, 1) = -1;//Ci+1
							// rt.RTcoef(idom, 2) = 0; // # Ci-1
							// rt.RTcoef(idom, 3) = 0;//Cj+1
							// rt.RTcoef(idom, 4) = 0;//Cj-1
							// rt.RTcoef(idom, 5) = 0;// Ck+1 
							// rt.RTcoef(idom, 6) = 0;// Ck-1 
							// rt.RTcoef(idom, 7) = 0;							
						}
						else if (direction == 3) {
							for (int i = 0; i < 30; i++) {
								rt.RTcoef(idom, i) = 0.0;
							}
							rt.RTcoef(idom, 0) = 1.0;		
							// rt.RTcoef(idom, 1) = -1;//Ci+1
							// rt.RTcoef(idom, 2) = 0; // # Ci-1
							// rt.RTcoef(idom, 3) = 0;//Cj+1
							rt.RTcoef(idom, 4) = -1;//Cj-1
							// rt.RTcoef(idom, 5) = 0;// Ck+1 
							// rt.RTcoef(idom, 6) = 0;// Ck-1 
							// rt.RTcoef(idom, 7) = 0;	
						}
						else if (direction == 4) {
							for (int i = 0; i < 30; i++) {
								rt.RTcoef(idom, i) = 0.0;
							}
							rt.RTcoef(idom, 0) = 1.0;		
							// rt.RTcoef(idom, 1) = -1;//Ci+1
							// rt.RTcoef(idom, 2) = 0; // # Ci-1
							rt.RTcoef(idom, 3) = -1;//Cj+1
							// rt.RTcoef(idom, 4) = 0;//Cj-1
							// rt.RTcoef(idom, 5) = 0;// Ck+1 
							// rt.RTcoef(idom, 6) = 0;// Ck-1 
							// rt.RTcoef(idom, 7) = 0;	
						}
						else if (direction == 5) {
							for (int i = 0; i < 30; i++) {
								rt.RTcoef(idom, i) = 0.0;
							}							
							rt.RTcoef(idom, 0) = 1.0;		
							// rt.RTcoef(idom, 1) = 0;//Ci+1
							// rt.RTcoef(idom, 2) = 0; // # Ci-1
							// rt.RTcoef(idom, 3) = 0;//Cj+1
							// rt.RTcoef(idom, 4) = 0;//Cj-1
							// rt.RTcoef(idom, 5) = 0;// Ck+1 
							rt.RTcoef(idom, 6) = -1;// Ck-1 
							// rt.RTcoef(idom, 7) = 0;
						}
						else if (direction == 6) {
							rt.RTcoef(idom, 6) = 0.0;
							for (int i = 0; i < 30; i++) {
								rt.RTcoef(idom, i) = 0.0;
							}
							rt.RTcoef(idom, 0) = 1.0;
							rt.RTcoef(idom, 5) = -1;

						}
					});
					break;
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

