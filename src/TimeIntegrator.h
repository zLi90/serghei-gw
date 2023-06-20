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

public :

  inline void stepForward(State &state, SourceSinkData &ss, std::vector<ExtBC> &extbc, Domain &dom, Exchange &exch, Parallel &par, FileIO &io, SergheiTimers &timers) {
	  timers.swe.reset();
        #if !SERGHEI_SUBSURFACE_MODEL
	    computeDt(state,dom,io);
        real dtloc = dom.dt;
        // int ierr = MPI_Allreduce(&dtloc, &dom.dt, 1, MPI_DOUBLE , MPI_MIN, MPI_COMM_WORLD);
        #endif
        //dom.dt = 0.1;

		edge.computeDeltaStateSW(state, dom, exch, par);

        // #if !SERGHEI_SUBSURFACE_MODEL
		// adjustDt(dom,io);
        // #endif
		timers.Tsweflux += timers.swe.seconds();

      ss.ComputeSWSourceSink(state,dom);

	   timers.swe.reset();
   	computeNewState(state, dom, ss);
		timers.Tswe += timers.swe.seconds();

		for (int k = 0; k < extbc.size(); k ++) { //should be done before the exchange (water depth might be modified).
		  applyExtBC(state,extbc[k],dom);
		}

		exch.exchangeMPIh(state,dom,exch,par); //only neccesary to exchange the h (for wet-dry) but for the moment we exchange everything
	   timers.swe.reset();

		wetDryCorrection( state, dom);

		timers.Tswe += timers.swe.seconds();

		exch.exchangeMPIhuhv(state,dom,exch,par);//neccesary to exchange again because of the wet/dry correction
	   timers.swe.reset();

		for (int k = 0; k < extbc.size(); k ++) { //after getting the final values, the discharge is integrated at every BC. The reason for not doing this before is because the previous kernels could eventually modify the boundary cell values.
			integrateExtBC(state,extbc[k],dom);
		}

		timers.Tswe += timers.swe.seconds();


        // Kokkos::parallel_for( dom.ncells , KOKKOS_LAMBDA (int iGlob) {
        //   int i, j, ncells;
    	// 	int id1,id2;
        //     if (iGlob == 124)  {
        //         unpackIndices(iGlob,dom.ny+2*hc,dom.nx+2*hc,j,i);
        //         //id2 = iGlob - dom.nx - 2*hc;
        //         id2 = iGlob;
        //         printf(" (%d,%d) : H = %f, V = %f, Q = %f\n",i,j,state.h(id2),state.hv(id2),20.0*60.0*state.hv(id2));
        //     }
        // });

	}

  inline void computeNewState(State &state , const Domain &dom, const SourceSinkData &ss) {
    #if DEBUG_WORKFLOW
    std::cout << GGD << __func__ << std::endl;
    #endif
    Kokkos::parallel_for( dom.nCellDomain , KOKKOS_LAMBDA (int iGlob) {
    int ii = dom.getIndex(iGlob);

		real z=state.z(ii);
		real hold=state.h(ii);
		real huold=state.hu(ii);
		real hvold=state.hv(ii);
		bool nodata=state.isnodata(ii);

		real hf,huf,hvf;
      int ncells =dom.ncells;

		hf = hold - dom.dt * (state.dsw0(ii)+state.dsw1(ii))/dom.dx;

        // for now, disable surface rainfall when subsurface is activated, ZhiLi20220414
        // #if !SERGHEI_SUBSURFACE_MODEL
		if(dom.isRain) {hf += ss.rainRate(ii)*dom.dt;}
        // #endif

    if(ss.inf.model) {

	 	ss.inf.rate(ii)=fmin(ss.inf.rate(ii),hf/dom.dt); //correct infiltration arte according to the available water volume
		ss.inf.rate(ii)=fmax(ss.inf.rate(ii),0.0); //avoid negative (in the order of machine accuracy) infiltration rates

      hf -= ss.inf.rate(ii)*dom.dt;
    }
    #if SERGHEI_SUBSURFACE_MODEL
    // printf(" SW -%d- : hf=%f, qss=%f, qt=%f\n",ii,1e5*hf,1e5*state.qss(ii),1e5*state.qss(ii) * dom.dt);
    hf += state.qss(ii) * dom.dt;
    #endif

		if(hf<TOL12 || nodata){
			//reduction or remove. Should be in the order of machine accuracy
			hf=0.0;

		}
		if(hf<state.hmin){
			huf=0.0;
			hvf=0.0;
		}else{
			real mx= huold - (state.dsw0(ii+ncells)+state.dsw1(ii+ncells))*dom.dt/dom.dx;
			real my= hvold - (state.dsw0(ii+2*ncells)+state.dsw1(ii+2*ncells))*dom.dt/dom.dx;
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

		if(fabs(huf)<TOL12){
			huf=0.0;
		}
		if(fabs(hvf)<TOL12){
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


    });

  }


  inline void computeGwExchange(State &state , const Domain &dom) {
    #if DEBUG_WORKFLOW
    std::cout << GGD << __func__ << std::endl;
    #endif
    Kokkos::parallel_for( dom.nCellDomain , KOKKOS_LAMBDA (int iGlob) {
        int ii = dom.getIndex(iGlob);
        state.h(ii) += state.qss(ii) * dom.dt;
    	if(state.h(ii)<TOL12){state.h(ii)=0.0;}
    });
  }



  inline void computeDt(State &state, Domain &dom, FileIO &io) {

  dom.dt = 1.e7;

	Kokkos::parallel_reduce( dom.nCellDomain , KOKKOS_LAMBDA (int iGlob, real &dt) {
    int ii = dom.getIndex(iGlob);
	 	real h=state.h(ii);
	 	real hu=state.hu(ii);
	 	real hv=state.hv(ii);
		dt=fmin(dt,1e3);
		if(h>TOL12){
			dt=fmin(dt,dom.dx/(fabs(hu/h)+sqrt(GRAV*h)));
			dt=fmin(dt,dom.dx/(fabs(hv/h)+sqrt(GRAV*h)));
		}
	} , Kokkos::Min<real>(dom.dt) );

	  Kokkos::fence();

	  dom.dt*=dom.cfl;

#if DEBUG_DT
	std::cout << "time = " << dom.etime << "\tdt_cfl = " << dom.dt << std::endl;
#endif
	 if(dom.dt>1e2){ //it means that evertyhing is dry.

     if(dom.isRain){
       // if there is rain, we impose a time step equivalent for h=1
       // this is to make sure we capture the start of the rain
       // TODO: improve this using the known rainfall signal
       dom.dt=dom.dx/(1+sqrt(GRAV));
        #if DEBUG_DT
        	std::cout << "time = " << dom.etime << "\tdt_rain = " << dom.dt << std::endl;
        #endif
     }
    }
    // Zhi Li
    //if (dom.dt > 0.1)   {dom.dt = 0.1;}

    // correction to match output times
    if (dom.etime + dom.dt > io.numOut*io.outFreq) { dom.dt = io.numOut*io.outFreq - dom.etime; }
    if (dom.etime + dom.dt > dom.simLength) { dom.dt = dom.simLength - dom.etime; }

    // real dtloc = dom.dt;
    // int ierr = MPI_Allreduce(&dtloc, &dom.dt, 1, MPI_DOUBLE , MPI_MIN, MPI_COMM_WORLD);


#if DEBUG_DT
	std::cout << "time = " << dom.etime << "\tdt_cor = " << dom.dt << std::endl;
#endif


  }


	inline void wetDryCorrection(State &state, Domain &dom) {

    	Kokkos::parallel_for( dom.ncells , KOKKOS_LAMBDA (int iGlob) {
			int i, j;
      	unpackIndices(iGlob,dom.ny+2*haloc,dom.nx+2*haloc,j,i);
			if(i>haloc-1 && i<dom.nx+haloc && j>haloc-1 && j<dom.ny+haloc){ //computational domain without halo cells
				real hij=state.h(iGlob);
				real zij=state.z(iGlob);
				real isB=state.isBound(iGlob);
				if(hij >= state.hmin){
					if(((hij+zij< state.z(iGlob+1)) && state.h(iGlob+1)<TOL12) || ((hij+zij<state.z(iGlob-1)) && state.h(iGlob-1)<TOL12) || (isB==0 && state.isnodata(iGlob+1)) || (isB==0 && state.isnodata(iGlob-1) )){
						state.hu(iGlob)=0.0;
					}
					if(((hij+zij<state.z(iGlob+dom.nx+2*haloc)) && state.h(iGlob+dom.nx+2*haloc)<TOL12)  || ((hij+zij<state.z(iGlob-(dom.nx+2*haloc))) && state.h(iGlob-(dom.nx+2*haloc))<TOL12) || (isB==0 && state.isnodata(iGlob+dom.nx+2*haloc)) || (isB==0 && state.isnodata(iGlob-(dom.nx+2*haloc)))){
						state.hv(iGlob)=0.0;
					}
				}
			}

		});


  }

  inline void applyExtBC(State &state, ExtBC &extbc, Domain &dom) {

    real extraMass=0.0;

    if(extbc.ncellsBC > 0){
	     switch (extbc.bctype){
         default:
	        std::cerr << RERROR "Boundary type: " << extbc.bctype << " not recognised." << std::endl;
	        std::cerr << RERROR "No boundary condition applied." << std::endl;
	        exit(EXIT_FAILURE);
          break;

         case SWE_BC_CRITICAL: // critical flow boundary condition
	         Kokkos::parallel_for(extbc.ncellsBC , KOKKOS_LAMBDA (int iGlob){
		           int ii=extbc.bcells[iGlob];
		           real h=state.h(ii);
		           if( h>=state.hmin) {
                 real hu=state.hu(ii);
		             real hv=state.hv(ii);
		             real vel=1.0*sqrt(GRAV*h); //Froude 1.0 (critical)
		             hu=vel*h*extbc.normalx;
		             hv=vel*h*extbc.normaly;
		             state.hu(ii)=hu;
		             state.hv(ii)=hv;
		          }
				   });
        break;

        case SWE_BC_WSE_CONST: // constant free surface elevation
          Kokkos::parallel_reduce(extbc.ncellsBC, KOKKOS_LAMBDA (int iGlob, real &sumM){
            int ii = extbc.bcells[iGlob];
            real h = state.h(ii);
				real hu= state.hu(ii);
				real hv= state.hv(ii);
            real z = state.z(ii);
            state.h(ii) = fmax(extbc.bcvals(0) - z, 0.0); // enforce water depth positivity
            sumM += (state.h(ii)-h)*dom.dx*dom.dx;
		        //orientation wrt to the outflow normal direction
				real modQ=sqrt(hu*hu+hv*hv);
		        state.hu(ii)=extbc.normalx*modQ;
		        state.hv(ii)=extbc.normaly*modQ;
          }, Kokkos::Sum<real>(extraMass));
	       break;

	       case SWE_BC_H_CONST: // constant depth boundary condition
	        Kokkos::parallel_reduce(extbc.ncellsBC, KOKKOS_LAMBDA (int iGlob, real &sumM){
	           int ii = extbc.bcells[iGlob];
			       real h = state.h(ii);
					 real hu= state.hu(ii);
					 real hv= state.hv(ii);

             state.h(ii) = extbc.bcvals(0);
			       sumM += (state.h(ii)-h)*dom.dx*dom.dx;
			       //orientation wrt to the outflow normal direction
					real modQ=sqrt(hu*hu+hv*hv);
		        	state.hu(ii)=extbc.normalx*modQ;
		        	state.hv(ii)=extbc.normaly*modQ;

           }, Kokkos::Sum<real>(extraMass) );
        break;

        case SWE_BC_Q_CONST: // constant inflow discharge boundary condition
          Kokkos::parallel_reduce(extbc.ncellsBC, KOKKOS_LAMBDA (int iGlob, real &sumM){
          int ii = extbc.bcells[iGlob];
          real h = state.h(ii);
          real qbc_x = extbc.bcvals(1);
 			    real qbc_y = extbc.bcvals(2);

			    //check if water depth is subcritical (Froude
			    //number less than 0.99). Otherwise impose
			    //boundary condition water depth
			    real inletFr=0.99;
			    real hcr=cbrt((qbc_x*qbc_x+qbc_y*qbc_y)/(GRAV*inletFr*inletFr));

			    //TODO: add this mass in case h is less than hcr
			    //state.h(ii) = fmax(h, hcr);
			    if (hcr > h) state.h(ii) = extbc.bcvals(0);

          sumM += (state.h(ii)-h)*dom.dx*dom.dx;

			    state.hu(ii) = qbc_x;
			    state.hv(ii) = qbc_y;
        }, Kokkos::Sum<real>(extraMass));
        break;

        case SWE_BC_FREE_OUTFLOW: // zero gradient or free boundary
          /*Kokkos::parallel_for(extbc.ncellsBC , KOKKOS_LAMBDA (int iGlob){
            int ii=extbc.bcells[iGlob];
				    //orientation wrt to the outflow direction
					 real hu= state.hu(ii);
				 	 real hv= state.hv(ii);
					 real modQ=sqrt(hu*hu+hv*hv);
		          state.hu(ii)=extbc.normalx*modQ;
		          state.hv(ii)=extbc.normaly*modQ;

          });*/
        break;

        case SWE_BC_HZ_T_INLET: // stage hydrograph inlet
          {
          real hzBC = interpolateLinear(extbc.hydrograph,dom.etime);
          #if DEBUG_BOUNDARY
            std::cout << "time = " << dom.etime << "\th+z = " << hzBC << std::endl;
          #endif
	         Kokkos::parallel_reduce(extbc.ncellsBC, KOKKOS_LAMBDA (int iGlob, real &sumM){
             int ii = extbc.bcells[iGlob];
             real h = state.h(ii);
				 real hu= state.hu(ii);
				 real hv= state.hv(ii);
                 if (hzBC-state.z(ii) > 0.0)    {state.h(ii) = hzBC-state.z(ii);}
                 else {state.h(ii) = 0.0;}
             sumM += (state.h(ii)-h)*dom.dx*dom.dx;

			    //orientation wrt to the outflow normal direction
				 real modQ=sqrt(hu*hu+hv*hv);
		       state.hu(ii)=extbc.normalx*modQ;
		       state.hv(ii)=extbc.normaly*modQ;

				}, Kokkos::Sum<real>(extraMass) );
         }
        break;

        case SWE_BC_HZ_T_OUTLET: // stage hydrograph outlet
          {
          real hzBC = interpolateLinear(extbc.hydrograph,dom.etime);
          #if DEBUG_BOUNDARY
            std::cout << "time = " << dom.etime << "\th+z = " << hzBC << std::endl;
          #endif
	         Kokkos::parallel_reduce(extbc.ncellsBC, KOKKOS_LAMBDA (int iGlob, real &sumM){
             int ii = extbc.bcells[iGlob];
             real h = state.h(ii);
             if (hzBC-state.z(ii) > 0.0)    {state.h(ii) = hzBC-state.z(ii);}
             else {state.h(ii) = 0.0;}
             sumM += (state.h(ii)-h)*dom.dx*dom.dx;

			    //no orientation wrt to the outflow normal direction to allow tidal wave coming into the domain

				}, Kokkos::Sum<real>(extraMass) );
         }
        break;




	  } // end switch
  } // endif ncellsBC


	 extbc.adjustedVolume=extraMass; //extraMass per bc

  }


  inline void integrateExtBC(State &state, ExtBC &extbc, Domain &dom) {

    real outflowDischarge=0.0;
    real inflowDischarge=0.0;
    real totalDischarge=0.0;


    if(extbc.ncellsBC > 0){
			//discharge integration
			Kokkos::parallel_reduce(extbc.ncellsBC, KOKKOS_LAMBDA (int iGlob, real &sumD){
				int ii = extbc.bcells[iGlob];
                //printf(" BC cell = %d\n",ii);
		      if( state.h(ii)>=state.hmin) {
					//the integration is done over all boundary walls according to the outflow direction
					sumD += (state.hu(ii)*sgn(extbc.normalx) + state.hv(ii)*sgn(extbc.normaly)) * dom.dx;
				}


			}, Kokkos::Sum<real>(totalDischarge));

	     switch (extbc.bctype){
			case SWE_BC_CRITICAL:
			case SWE_BC_H_CONST:
			case SWE_BC_WSE_CONST:
			case SWE_BC_FREE_OUTFLOW:
      	case SWE_BC_HZ_T_OUTLET:

				outflowDischarge=totalDischarge;
			break;

			case SWE_BC_Q_CONST:
			case SWE_BC_HZ_T_INLET:
			case SWE_BC_Q_T:
				inflowDischarge=totalDischarge;

			break;

		  default:
			 std::cerr << RERROR " Boundary type: " << extbc.bctype << " not recognized." << std::endl;
			 std::cerr << RERROR " No boundary condition applied." << std::endl;
			 exit(EXIT_FAILURE);
		  }
	}


    extbc.inflowDischarge = inflowDischarge; //inflowdischarge local per bc
    extbc.outflowDischarge = outflowDischarge; //outflowdischarge local per bc
    extbc.inflowAccumulated += inflowDischarge*dom.dt; //inflowAccumulated local per bc
	 extbc.outflowAccumulated += outflowDischarge * dom.dt; //outflowaccumulated local per bc


	}
};

#endif
