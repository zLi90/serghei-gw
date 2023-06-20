/* -*- mode: c++ -*- */

#ifndef _BC_H_
#define _BC_H_

#include "define.h"
#include "Indexing.h"

// DCV 05.05.2021, left these here defined for generic use in exchange.h, not for hydraulics.
#define BC_PERIODIC 1
#define BC_REFLECTIVE 2
#define BC_TRANSMISSIVE 3

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

// subsurface bc types
#define SUB_BC_NOFLOW 1
#define SUB_BC_H_CONST 2
#define SUB_BC_Q_CONST 3
#define SUB_BC_H_FILE 4
#define SUB_BC_H_SWE 5


class ExtBC {

public:
	std::string id;
	int ncellsBC = 0; //number of bcells
	intArr bcells; //array of indexes of boundary cells
	real normalx, normaly; //direction set by user for inflow/outflow
	int location; //1->west, 2->north, 3->east, 4-> south 
	int bctype;
	int isInDomain;
  realArr bcvals;
	real outflowDischarge;
	real outflowAccumulated = 0;
	real inflowDischarge;
	real inflowAccumulated = 0;
	real adjustedVolume = 0;
	TimeSeries hydrograph;

  std::string BoundaryTypes[13] = {"NONE","PERIODIC","REFLECTIVE","TRANSMISSIVE","NONE","CRITICAL","CONSTANT DEPTH","CONSTANT INFLOW","CONSTANT WSELEVATION","FREE OUTFLOW","STAGE HYDROGRAPH INLET","STAGE HYDROGRAPH OUTLET", "HYDROGRAPH"};

	inline int find_bcells(State &state, ExtBC &extbc, Domain &dom, Parallel &par, int nPoly, realArr &xPoly, realArr &yPoly){
		std::vector<int> tmpbcells; //array of indexes of boundary cells

		for(int iGlob=0;iGlob<dom.nx*dom.ny;iGlob++){
			int i,j;
			unpackIndices(iGlob,dom.ny,dom.nx,j,i);
			int ii=(hc+j)*(dom.nx+2*hc)+hc+i;//index with the extended domain (including halo cells)

			if(!state.isnodata(ii)){
				if((j==0&&dom.iN) || (j==dom.ny-1&&dom.iS) || (i==0&&dom.iW) || (i==dom.nx-1&&dom.iE) ||
				state.isnodata(ii+1) || state.isnodata(ii-1) ||
				state.isnodata(ii-(dom.nx+2*hc)) || state.isnodata(ii+(dom.nx+2*hc))){
				//boundary domain || nodata neighbours
					real xCoord = dom.xll + ( par.i_beg + i + 0.5) * dom.dx;
					real yCoord = dom.yll + dom.ny_glob*dom.dx - ( par.j_beg + j + 0.5) * dom.dx;
					if(isInsidePoly(nPoly,xPoly, yPoly, xCoord, yCoord)){
						tmpbcells.push_back(ii);
						//it is important to add the outflow direction because there might be cells with double boundary walls
						//the value of -3000.0 is set as a very low number value. First, the elevation of the contiguous cell was imposed, but when dealing with no data and sawtooth pattern with extbc there might be some problems because a ghost cell can be the neighbour of more than one boundary cells with different elevations. Therefore, it is not clear which elevation is the best (should be the minimum of them to avoid wet/dry/solid wall problems. In this case, a very low number is chosen (-3000.0)
						if(state.isnodata(ii+1) && fabs(extbc.normalx)>0.0){
							state.z(ii+1)=-3000.0; //set the elevation of the neighbour (it was nodata) to allow the water to flow
						}
						if(state.isnodata(ii-1)&& fabs(extbc.normalx)>0.0){
							state.z(ii-1)=-3000.0; //set the elevation of the neighbour (it was nodata) to allow the water to flow
						}
						if(state.isnodata(ii-(dom.nx+2*hc)) && fabs(extbc.normaly)>0.0){
							state.z(ii-(dom.nx+2*hc))=-3000.0; //set the elevation of the neighbour (it was nodata) to allow the water to flow
						}
						if(state.isnodata(ii+(dom.nx+2*hc))&& fabs(extbc.normaly)>0.0){
							state.z(ii+(dom.nx+2*hc))=-3000.0; //set the elevation of the neighbour (it was nodata) to allow the water to flow
						}

					}

				}
			}
		}
		extbc.ncellsBC=int(tmpbcells.size());

		int ncells_all;
      MPI_Allreduce(&extbc.ncellsBC, &ncells_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		//we need the total boundary cells detected by all subdomain to launch an error otherwise
		if(ncells_all>0){
			extbc.bcells=intArr("bcells", extbc.ncellsBC);
			#ifdef __NVCC__
				cudaMemcpyAsync( extbc.bcells.data() , tmpbcells.data() , extbc.ncellsBC*sizeof(int) , cudaMemcpyHostToDevice );
				cudaDeviceSynchronize();
			#else
				std::memcpy(extbc.bcells.data(), tmpbcells.data(), extbc.ncellsBC*sizeof(int));
			#endif
		}
		else{
			if(par.masterproc){
				std::cerr << RERROR << "No boundary cells found for external boundary with id '" << id << "'" << std::endl;
			}
			return 0;
		}

		return 1;
	}



	//returns 1 if the coordinate is inside, 0 otherwise
	int isInsidePoly(int np, realArr &xPoly, realArr &yPoly, real &xCoord, real &yCoord) {

		int count;
		int i;
		double xInter;
		real x1,x2,y1,y2;

		count=0;
		for (i=0;i<np;i++){
			x1= xPoly(i);
			y1= yPoly(i);
			x2= xPoly((i+1)%np);
			y2= yPoly((i+1)%np);
			if ((yCoord > fmin(y1,y2)) && (yCoord <= fmax(y1,y2)) && (xCoord <= fmax(x1,x2))) {
			  if (y1 != y2) {
				 xInter = x1 + (yCoord-y1)*(x2-x1)/(y2-y1);
				 if (x1 == x2 || xCoord <= xInter){
						count++;
				 }
			  }
			}
		}

		if(count % 2 == 0){
			return 0;
		}else{
			return 1;
		}
	}



};


#endif
