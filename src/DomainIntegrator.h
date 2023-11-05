#ifndef _DOMAIN_INTEGRATOR_H_
#define _DOMAIN_INTEGRATOR_H_

#include "State.h"
#include "SWSourceSink.h"
#include "Domain.h"
#include "Indexing.h"

class surfaceIntegrator {

  Kokkos::Timer timer;

  public:

  // integrated variables
  real surfaceVolume ;     // surface water volume in domain [L^3] (local)
  real rainFlux ;   // total rain flux [L^3 / T] (local)
  real rainAccum = 0.;   // accumulated rainfall in simulation [L^3] (local)
  real infFlux ;    // total infiltration flux [L^3/T] (local)
  real infAccum = 0.;    // accumulated infiltration in simulation [L^3] (local)

  real surfaceVolumeG ;     // surface water volume in domain [L^3] (global)
  real rainFluxG ;   // total rain flux [L^3 / T] (global)
  real rainAccumG = 0.;   // accumulated rainfall in simulation [L^3] (global)
  real infFluxG ;    // total infiltration flux [L^3/T] (global)
  real infAccumG = 0.;    // accumulated infiltration in simulation [L^3] (global)

  // pointers
  SourceSinkData *ss;
  State *state;
  Domain *dom;


  void initialize(State &state_, Domain &dom_, SourceSinkData &ss_){
    state = &state_;
    dom = &dom_;
    ss = &ss_;
  }

  void integrate(State const &state, Domain const &dom, SourceSinkData &ss){
    timer.reset();
    real A = dom.dx * dom.dx;  // WARNING UCM  - Uniform Cartesian Mesh

	 int ierr=0;

    surfaceVolume = 0;
	  Kokkos::parallel_reduce( dom.nCellDomain , KOKKOS_LAMBDA (int iGlob, real &valUpdate) {
        int ii = dom.getIndex(iGlob);
		  bool nodata=state.isnodata(ii);
		  if(!nodata)
        valUpdate +=  state.h(ii) * A;
    } , Kokkos::Sum<real>(surfaceVolume) );
		Kokkos::fence();

	if(dom.etime<TOL12){ //change by initial time when hotstart is implemented

	 	rainFlux=0.0;
	 	rainAccum=0.0;
	 	infFlux=0.0;
	 	infAccum=0.0;

	}else{

		 rainFlux = 0.0;
		 if(dom.isRain){
			Kokkos::parallel_reduce( dom.nCellDomain , KOKKOS_LAMBDA (int iGlob, real &valUpdate) {
			  int ii = dom.getIndex(iGlob);
			  bool nodata=state.isnodata(ii);
			  if(!nodata)
			  valUpdate += ss.rainRate(ii) * A;
			} , Kokkos::Sum<real>(rainFlux) );
			  Kokkos::fence();
		    rainAccum += rainFlux * dom.dt;
		 }

		 infFlux = 0.0;
		 if(ss.inf.model){
			Kokkos::parallel_reduce( dom.nCellDomain , KOKKOS_LAMBDA (int iGlob, real &valUpdate) {
			  real inffluxlocal;
			  int ii = dom.getIndex(iGlob);
			  bool nodata=state.isnodata(ii);
			  if(!nodata)
			  inffluxlocal = ss.inf.rate(ii)*A;
			  valUpdate += inffluxlocal;
			  ss.inf.infVol(ii) += inffluxlocal;
			} , Kokkos::Sum<real>(infFlux) );
			  Kokkos::fence();
		 }
		 infAccum += infFlux * dom.dt;

	}

	surfaceVolumeG=0.0;
	rainFluxG=0.0;
	rainAccumG=0.0;
	infFluxG=0.0;
	infAccumG=0.0;
	ierr=MPI_Allreduce(&surfaceVolume, &surfaceVolumeG, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
	ierr=MPI_Allreduce(&rainFlux, &rainFluxG, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
	ierr=MPI_Allreduce(&rainAccum, &rainAccumG, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
	ierr=MPI_Allreduce(&infFlux, &infFluxG, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
	ierr=MPI_Allreduce(&infAccum, &infAccumG, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

  dom.timers.integrate += timer.seconds();
  }
};

class boundaryIntegrator
{

Kokkos::Timer timer;

public:

  int ncellsBC, ncellsBCG;

  real adjustedVolume ;  // boundary water volume in domain [L^3] adjusted (e.g. impose water depth) (local)
  real outflowDischarge; // boundary outflow discharge in domain [L^3/T] (local)
  real outflowAccumulated = 0.0; // boundary accumulated outflow volume in domain [L^3] (local)
  real inflowDischarge; // boundary inflow discharge in domain [L^3/T] (local)
  real inflowAccumulated = 0.0; // boundary accumulated inflow volume in domain [L^3] (local)

  real adjustedVolumeG ;  // boundary water volume in domain [L^3] adjusted (e.g. impose water depth) (global)
  real outflowDischargeG; // boundary outflow discharge in domain [L^3/T] (global)
  real outflowAccumulatedG = 0.0; // boundary accumulated outflow volume in domain [L^3] (global)
  real inflowDischargeG; // boundary inflow discharge in domain [L^3/T] (global)
  real inflowAccumulatedG = 0.0; // boundary accumulated inflow volume in domain [L^3] (global)


  std::vector<ExtBC>* extbc;

  void initialize (std::vector<ExtBC> &extbc_)
  {
    extbc = &extbc_;
  }

  void integrate (std::vector<ExtBC> &extbc, Domain const &dom, int mode){
    timer.reset();

	 int ierr=0;

	 //mode is a flag to integrate extra mass or boundary flows

	 if(mode==0){

	 	adjustedVolume = 0.0;
		for (int i = 0; i < extbc.size(); i ++) {
			adjustedVolume += extbc[i].adjustedVolume;
		}
	 	adjustedVolumeG = 0.0;
		ierr=MPI_Allreduce(&adjustedVolume, &adjustedVolumeG, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);

	 }else{

		 ncellsBC = 0;
		 int _ncellsBC;
		 inflowDischarge = 0.0;
		 outflowDischarge = 0.0;
     inflowAccumulated= 0.0;
     outflowAccumulated=0.0;


		 // integrate over all the open external boundaries
		 // no MPI reduction is necessary, as they flows and volumes are already computed per open boundary in ExtBC::integrate
		 for (int i = 0; i < extbc.size(); i ++) {
			MPI_Allreduce(&(extbc[i].ncellsBC), &_ncellsBC, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			ncellsBC += extbc[i].ncellsBC;
			inflowDischarge += extbc[i].inflowDischarge;
			inflowAccumulated+= extbc[i].inflowAccumulated;
			outflowDischarge += extbc[i].outflowDischarge;
			outflowAccumulated+= extbc[i].outflowAccumulated;
			#if SERGHEI_DEBUG_BOUNDARY
				std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "\tBC " << i << "\tinflowDischarge = " << extbc[i].inflowDischarge << std::endl;
				std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "\tinflowDischarge = " << inflowDischarge << std::endl;
				std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "\tBC " << i << "\tinflowDischarge = " << extbc[i].outflowDischarge << std::endl;
				std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "\toutflowDischarge = " << outflowDischarge << std::endl;
			#endif
		 }

		inflowDischargeG=0.0;
		outflowDischargeG=0.0;
		inflowAccumulatedG=0.0;
		outflowAccumulatedG=0.0;
        ncellsBCG = 0;
		ierr=MPI_Allreduce(&inflowDischarge, &inflowDischargeG, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
		ierr=MPI_Allreduce(&outflowDischarge, &outflowDischargeG, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
		ierr=MPI_Allreduce(&inflowAccumulated, &inflowAccumulatedG, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
		ierr=MPI_Allreduce(&outflowAccumulated, &outflowAccumulatedG, 1, SERGHEI_MPI_REAL , MPI_SUM, MPI_COMM_WORLD);
        ierr=MPI_Allreduce(&ncellsBC, &ncellsBCG, 1, MPI_INT , MPI_SUM, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}

  dom.timers.integrate += timer.seconds();
  }

};

#endif
