/* -*- mode: c++ -*- */

#ifndef _TIMEINTEGRATOR_H_
#define _TIMEINTEGRATOR_H_

#include <stdlib.h>

#include "define.h"
#include "Parallel.h"
#include "Domain.h"
#include "State.h"
#include "BC.h"
#include "Edges.h"
#include "Indexing.h"
#include "FileIO.h"
#include "SWSourceSink.h"

class TimeIntegrator {

  Edges edge;
  Kokkos::Timer timer;

public :

  inline void stepForward(State &state, SourceSinkData &ss, std::vector<ExtBC> &extbc, Domain &dom, Exchange &exch, Parallel &par, FileIO &io) {
	computeDt(state,dom,io);

		edge.computeDeltaStateSW(state, dom, exch, par);

    ss.ComputeSWSourceSink(state,dom);

   	computeNewState(state, dom, ss);

		for (int k = 0; k < extbc.size(); k ++) { //should be done before the exchange (water depth might be modified).
		  extbc[k].apply(state,dom);
		}

		exch.exchangeMPIh(state,dom,exch,par); //only neccesary to exchange the h (for wet-dry) but for the moment we exchange everything

		wetDryCorrection( state, dom);


		exch.exchangeMPIhuhv(state,dom,exch,par);//neccesary to exchange again because of the wet/dry correction

		for (int k = 0; k < extbc.size(); k ++) { //after getting the final values, the discharge is integrated at every BC. The reason for not doing this before is because the previous kernels could eventually modify the boundary cell values.
			extbc[k].integrate(state,dom);
		}

	}

inline void computeNewState(State &state , const Domain &dom, const SourceSinkData &ss) {
    #if SERGHEI_DEBUG_WORKFLOW
    std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
    #endif
    timer.reset();
    Kokkos::parallel_for("computeNewState",dom.nCell , KOKKOS_LAMBDA (int iGlob) {
    int ii = dom.getIndex(iGlob);

		real z=state.z(ii);
		real hold=state.h(ii);
		real huold=state.hu(ii);
		real hvold=state.hv(ii);
		bool nodata=state.isnodata(ii);

		real hf,huf,hvf;
      int ncells =dom.nCellMem;

		hf = hold - dom.dt * (state.dsw0(ii)+state.dsw1(ii))/dom.dx();

		if(dom.isRain) hf += ss.rainRate(ii)*dom.dt;
    if(ss.inf.model) {

	 	ss.inf.rate(ii)=min(ss.inf.rate(ii),hf/dom.dt); //correct infiltration arte according to the available water volume
		ss.inf.rate(ii)=max(ss.inf.rate(ii),0.0); //avoid negative (in the order of machine accuracy) infiltration rates

      hf -= ss.inf.rate(ii)*dom.dt;
    }

		if(hf<TOL_MACHINE_ACCURACY || nodata){
			//reduction or remove. Should be in the order of machine accuracy
			hf=0.0;

		}
		if(hf<state.hmin){
			huf=0.0;
			hvf=0.0;
		}else{
			real mx= huold - (state.dsw0(ii+ncells)+state.dsw1(ii+ncells))*dom.dt/dom.dx();
			real my= hvold - (state.dsw0(ii+2*ncells)+state.dsw1(ii+2*ncells))*dom.dt/dom.dx();
			#if POINTWISE_FRICTION
				real nsq= state.roughness(ii)*state.roughness(ii);
				real modM=sqrt(mx*mx/hold/hold+my*my/hold/hold);
				if(nsq>TOL12 && hold>=state.hmin && modM>TOL12){
					real tt=dom.dt*GRAV*nsq*modM/(hold*cbrt(hold));
					real ff=sqrt(1.0+4.0*tt);
					huf=-0.5*(mx-mx*ff)/tt;
					hvf=-0.5*(my-my*ff)/tt;
				}else{
					huf=mx;
					hvf=my;
				}
			#else
				huf=mx;
				hvf=my;
			#endif


		}

		if(fabs(huf)<TOL_ZERO_MOMENTUM){
			huf=0.0;
		}
		if(fabs(hvf)<TOL_ZERO_MOMENTUM){
			hvf=0.0;
		}

		state.h(ii)=hf;
		state.hu(ii)=huf;
		state.hv(ii)=hvf;

    //reset the contributions
		state.dsw0(ii)=0.0;
		state.dsw0(ii+ncells)=0.0;
		state.dsw0(ii+2*ncells)=0.0;
		state.dsw1(ii)=0.0;
		state.dsw1(ii+ncells)=0.0;
		state.dsw1(ii+2*ncells)=0.0;

		#if SERGHEI_MAXFLOOD
			if(hf > state.hMax(ii)){
				state.hMax(ii) = hf;
				state.time_hMax(ii) = dom.etime;
			}
			real mom = mysqrt(huf*huf+hvf*hvf);
			if(mom > state.momentumMax(ii)) state.momentumMax(ii) = mom;
		#endif

    });

		dom.timers.swe += timer.seconds();
  }

  inline void computeDt(State &state, Domain &dom, FileIO &io) {
  timer.reset();

  dom.dt = 1.e7;

	Kokkos::parallel_reduce("reduceDt",dom.nCell , KOKKOS_LAMBDA (int iGlob, real &dt) {
    int ii = dom.getIndex(iGlob);
	 	real h=state.h(ii);
	 	real hu=state.hu(ii);
	 	real hv=state.hv(ii);
		dt=min(dt,1.e6);
		if(h>TOL12){
			dt=min(dt,dom.dx()/(fabs(hu/h)+sqrt(GRAV*h)));
			dt=min(dt,dom.dx()/(fabs(hv/h)+sqrt(GRAV*h)));
		}
	} , Kokkos::Min<real>(dom.dt) );

	  Kokkos::fence();

	  real dtloc = dom.dt;
	  int ierr = MPI_Allreduce(&dtloc, &dom.dt, 1, SERGHEI_MPI_REAL , MPI_MIN, MPI_COMM_WORLD);

	  dom.dt*=dom.cfl;

#if SERGHEI_DEBUG_DT
	std::cout << "time = " << dom.etime << "\tdt_cfl = " << dom.dt << "\tdom.cfl = " << dom.cfl << std::endl;
#endif
	 if(dom.dt>1e5){ //it means that evertyhing is dry.
     if(dom.isRain){
       // if there is rain, we impose a time step equivalent for h=1
       // this is to make sure we capture the start of the rain
       // TODO: improve this using the known rainfall signal
       dom.dt=dom.dx()/(1+sqrt(GRAV));		// dx() purposely used here to fail upon compilation when moving to adaptive mesh
#if SERGHEI_DEBUG_DT
	std::cout << "time = " << dom.etime << "\tdt_rain = " << dom.dt << std::endl;
#endif
     }
    }
    // correction to match output times
		//std::cout << GGD << dom.etime << "\t" << dom.etime+dom.dt << "\t" << io.numOut << "\t" << io.outFreq*io.numOut << "\t" << io.numOut*io.outFreq + dom.startTime << std::endl;
    if (dom.etime + dom.dt > dom.startTime +  io.numOut*io.outFreq) dom.dt = io.numOut*io.outFreq + dom.startTime - dom.etime; 
    if (dom.etime + dom.dt > dom.endTime) { dom.dt = dom.endTime - dom.etime; }


#if SERGHEI_DEBUG_DT
	std::cout << "time = " << dom.etime << "\tdt_cor = " << dom.dt << std::endl;
#endif

  dom.timers.dt=timer.seconds();
  }


	inline void wetDryCorrection(State &state, Domain &dom) {
	#if SERGHEI_DEBUG_WORKFLOW
    	std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
	#endif
    timer.reset();

    	Kokkos::parallel_for("wetDryCorrection", dom.nCellMem , KOKKOS_LAMBDA (int iGlob) {
			int i, j;
      	unpackIndicesUniformGrid(iGlob,dom.ny+2*hc,dom.nx+2*hc,j,i);
			if(i>hc-1 && i<dom.nx+hc && j>hc-1 && j<dom.ny+hc){ //computational domain without halo cells
				real hij=state.h(iGlob);
				real zij=state.z(iGlob);
				int isB=state.isBound(iGlob);
				if(hij >= state.hmin){
					if(((hij+zij< state.z(iGlob+1)) && state.h(iGlob+1)<TOL_WETDRY) || ((hij+zij<state.z(iGlob-1)) && state.h(iGlob-1)<TOL_WETDRY) || (isB==0 && state.isnodata(iGlob+1)) || (isB==0 && state.isnodata(iGlob-1) )){
						state.hu(iGlob)=0.0;
					}
					if(((hij+zij<state.z(iGlob+dom.nx+2*hc)) && state.h(iGlob+dom.nx+2*hc)<TOL_WETDRY)  || ((hij+zij<state.z(iGlob-(dom.nx+2*hc))) && state.h(iGlob-(dom.nx+2*hc))<TOL_WETDRY) || (isB==0 && state.isnodata(iGlob+dom.nx+2*hc)) || (isB==0 && state.isnodata(iGlob-(dom.nx+2*hc)))){
						state.hv(iGlob)=0.0;
					}
				}
			}

		});

    dom.timers.swe += timer.seconds();
  }


}; // end of TimeIntegrator class

#endif
