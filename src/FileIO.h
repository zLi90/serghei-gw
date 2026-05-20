#ifndef _FILEIO_H_
#define _FILEIO_H_

#include <chrono>
#include <ctime>
#include <unistd.h>
#ifndef __NVCC__ 
#include <cpuid.h>
#endif
#include "const.h"
#include "define.h"
#include "State.h"
#include "SourceSink.h"
#include "BC.h"
#include "pnetcdf.h"
#include "mpi.h"
#include "Indexing.h"
#include "DomainIntegrator.h"
#include "tools.h"
#include "ParticleTracking.h"

#include "GwState.h"
#include "GwDomain.h"
#include "GwIntegrator.h"
#if SERGHEI_DRAINAGE_MODEL
#include "Drainage/DrainageDywave.h"
#include "Drainage/funcs.h"
#endif

#ifndef SERGHEI_NC_MODE
#define SERGHEI_NC_MODE NC_CLOBBER
#endif

#ifndef SERGHEI_WRITE_HZ
#define SERGHEI_WRITE_HZ 1
#endif
#ifndef SERGHEI_WRITE_SUBDOMS
#define SERGHEI_WRITE_SUBDOMS 0
#endif

#ifndef SERGHEI_NC_REAL
  #if SERGHEI_REAL == SERGHEI_DOUBLE
    #define SERGHEI_NC_REAL NC_DOUBLE
  #elif SERGHEI_REAL == SERGHEI_FLOAT
    #define SERGHEI_NC_REAL NC_FLOAT
  #endif
#endif

#if SERGHEI_NC_REAL == NC_DOUBLE
	#define ncmpi_put_vara_real_all ncmpi_put_vara_double_all
  #define ncmpi_put_att_real ncmpi_put_att_double
#endif
#if SERGHEI_NC_REAL == NC_FLOAT
	#define ncmpi_put_vara_real_all ncmpi_put_vara_float_all
  #define ncmpi_put_att_real ncmpi_put_att_float
#endif

#define SERGHEI_NC_FILL_VALUE_KEY "_FillValue"
#ifndef SERGHEI_NC_ENABLE_NAN
#define SERGHEI_NC_ENABLE_NAN 1
#endif
#ifndef SERGHEI_INPUT_NETCDF
#define SERGHEI_INPUT_NETCDF 0

#endif

#ifdef __NVCC__
	#if SERGHEI_NC_REAL==NC_FLOAT && SERGHEI_REAL==SERGHEI_DOUBLE
  	typedef Kokkos::View<float*     ,Kokkos::LayoutRight,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaUVMSpace>> ncArr;
	#else
  	typedef Kokkos::View<real*     ,Kokkos::LayoutRight,Kokkos::Device<Kokkos::Cuda,Kokkos::CudaUVMSpace>> ncArr;
	#endif
#else
	#if SERGHEI_NC_REAL==NC_FLOAT && SERGHEI_REAL==SERGHEI_DOUBLE
  	typedef Kokkos::View<float*     ,Kokkos::LayoutRight> ncArr;
	#else
  	typedef Kokkos::View<real*     ,Kokkos::LayoutRight> ncArr;
	#endif
#endif
#if SERGHEI_NC_REAL==NC_FLOAT
	typedef float ncreal;
#elif SERGHEI_NC_REAL==NC_DOUBLE
	typedef double ncreal;
#endif

/*
#define EMPTY()
#define DEFER(x) x EMPTY()
#define PASTER(x,y) x ## y
#define EVALUATOR(x,y) PASTER(x,y)
#define CHOOSE_BACKEND(name,id) EVALUATOR(name,id)
*/



class ncStream{
	public:
	std::string fname;
	int id;
	int ndims;
	int nvars;
	int ngatts;
	int unlimited;
	int dimids[3];
	int ndata;

};

class FileIO {
public: 
	ncStream ncin; 
  bool writeRain = 0;


protected:

  int ncid;
  int tDim, xDim, yDim, zDim;
//   int tVar, xVar, yVar, hVar, hzVar, uVar, vVar, zVar, z3Var, hdVar, wcVar, qVar;
  int tVar, xVar, yVar, hVar, hzVar, uVar, vVar, zVar, z3Var, hdVar, wcVar, qVar, DqVar, VqVar;
  int infVar,infVolVar,infGaFVar;
  int rainVar;
  std::ofstream domainOutputFile;
  std::ofstream SubsurfaceOutputFile;
  std::ofstream DrainageOutputFile;
  std::ofstream logFile;
	int nOut; // expected number of spatial outputs
	#if SERGHEI_MAXFLOOD > 0
	int hMaxVar, momMaxVar, timehMaxVar;
	#endif
	#if SERGHEI_WRITE_SUBDOMS
		int subdomVar;
	#endif
	#if SERGHEI_LPT
		std::string lpt_filename;
		int tpDim,pDim;
		int ncidP,pVar,tpVar,xpVar,ypVar,zpVar;
		std::ofstream particlesFile;
	#endif


private:
  Kokkos::Timer timer;
  #if SERGHEI_LPT
  Kokkos::Timer timer_particles;
  #endif

/*
	//TODO probably should be a C++ template
	inline int ncmpi_put_vara_real_all( int ncid , int varid , const MPI_Offset start[], const MPI_Offset count[], const real *buf){
		#if SERGHEI_NC_REAL == NC_DOUBLE
    	return( ncmpi_put_vara_double_all( ncid , varid, start , count , buf ));
		#endif
		#if SERGHEI_NC_REAL == NC_FLOAT
			#if SERGHEI_REAL == SERGHEI_DOUBLE

				return( ncmpi_put_vara_float_all( ncid , varid, start , count , (float*) buf));
			#else
				return( ncmpi_put_vara_float_all( ncid , varid, start , count , buf));
			#endif
		#endif
	}
*/

public:
	bool allowIni=1;
  real outFreq;
  int outFormat;
  int numOut=0; // spatial output counter

  real obsFreq;
  int numObs = 0;
  int nScreen = 100; //if this value is not specified, the information is displayed every 1000 iterations

// Writes spatial fields for initial state
	void outputIni(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par, std::string dir){
		nOut = floor(dom.simLength / outFreq);
		numOut=0;
		if(outFormat==OUT_NETCDF){
			outputInitNETCDF(state,dom,ss,par,dir);
		}
		if(outFormat==OUT_VTK){
			outputVTK(state,dom,ss,par,dir);
		}
		if(outFormat==OUT_BIN){
			initBIN(state,dom,par,dir);
			outputBIN(state,dom,par,dir);
		}
		numOut++;
	}

  #if SERGHEI_SUBSURFACE_MODEL
  void outputIniSub(const GwState &gw, GwDomain const &gdom, Parallel const &par, std::string dir){
		nOut = floor(gdom.simLength / outFreq);
		outputInitNETCDFSub(gw, gdom, par, dir);
	}

  void outputSubsurface(const GwState &gw, GwDomain const &gdom, Parallel const &par, std::string dir){
    timer.reset();
    #if SERGHEI_SWE_MODEL
    numOut--;
    #endif
    outputNETCDFSubsurface(gw, gdom, par, dir);
    numOut++;
    gdom.timers.out += timer.seconds();
	}
  #endif

	// Writes spatial fields
	#if SERGHEI_LPT==0
	void output(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par,std::string dir){
    #else
    void output(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par,std::string dir, ParticleTracker const &particles){
    #endif
    	timer.reset();

		if(outFormat==OUT_NETCDF){
			outputNETCDF(state,dom,ss,par,dir);
			#if SERGHEI_LPT
				timer_particles.reset();
				#if SERGHEI_PARTICLE_NO_OUTPUT
			    #else
            		outputNETCDF_particles(particles,dom,par,dir); //Particle Tracking
				#endif
				dom.timers.output_particles += timer_particles.seconds();
            #endif
		}
		if(outFormat==OUT_VTK){
			outputVTK(state,dom,ss,par,dir);
			
			#if SERGHEI_LPT
              timer_particles.reset();
			  #if SERGHEI_PARTICLE_NO_OUTPUT
			  #else
			  outputParticle(dom,particles,dir); //Particle Tracking
			  #endif
			  dom.timers.output_particles += timer_particles.seconds();
            #endif
		}
		if(outFormat==OUT_BIN){
			outputBIN(state,dom,par,dir);
		}
		#if SERGHEI_LPT
			timer_particles.reset();
			#if SERGHEI_PARTICLE_NO_OUTPUT
			#else
				#if SERGHEI_LPT
				writeParticle_TemporalFile(dom, particles, dir);
				#endif
			#endif
			dom.timers.output_particles += timer_particles.seconds();
		#endif

		numOut++;
    	dom.timers.out += timer.seconds();
	}

	void writeNetCDFfield(const Domain &dom, int ncid, int ncvar, MPI_Offset *st, MPI_Offset *ct, const boolArr &mask, const realArr &myview, ncArr &data){
		Kokkos::parallel_for("writeNCDFfield", dom.nCell , KOKKOS_LAMBDA(int iGlob) {
			int ii = dom.getIndex(iGlob);
      data(iGlob) = myview(ii);
			#if SERGHEI_NC_ENABLE_NAN
			if(mask(ii)) data(iGlob) = SERGHEI_NAN;
			#endif
  	});
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_real_all( ncid , ncvar  , st , ct , data.data() ) , __LINE__);
	};

	// NetCDF initialiser and writer
  void outputInitNETCDF(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par, std::string dir) {
    int dimids[3];
    MPI_Offset st[3], ct[3];
    doubleArr xCoord = doubleArr("xCoord",dom.nx);
    doubleArr yCoord = doubleArr("yCoord",dom.ny);
    ncArr data   = ncArr("data",dom.nCell);
		static char timeUnits[] = "seconds" ;
		std::string filename;
		std::string longname;
		std::string units;

	filename=dir+"output.nc";
    // Create the file
    ncwrap( ncmpi_create( MPI_COMM_WORLD , filename.c_str() , SERGHEI_NC_MODE , MPI_INFO_NULL , &ncid ) , __LINE__,par.myrank );

    // Create the dimensions
    ncwrap( ncmpi_def_dim( ncid , "t" , (MPI_Offset) NC_UNLIMITED , &tDim ) , __LINE__,par.myrank );
    ncwrap( ncmpi_def_dim( ncid , "x" , (MPI_Offset) dom.nx_glob  , &xDim ) , __LINE__,par.myrank );
    ncwrap( ncmpi_def_dim( ncid , "y" , (MPI_Offset) dom.ny_glob  , &yDim ) , __LINE__,par.myrank );
    // Create the variables
    dimids[0] = tDim;
    ncwrap( ncmpi_def_var( ncid , "t"      , NC_DOUBLE , 1 , dimids , &tVar ) , __LINE__,par.myrank );
	 	ncwrap( ncmpi_put_att_text (ncid, tVar, "units",strlen(timeUnits), timeUnits), __LINE__,par.myrank );
    dimids[0] = xDim;
    ncwrap( ncmpi_def_var( ncid , "x"      , NC_DOUBLE , 1 , dimids , &xVar ) , __LINE__,par.myrank );
    dimids[0] = yDim;
    ncwrap( ncmpi_def_var( ncid , "y"      , NC_DOUBLE , 1 , dimids , &yVar ) , __LINE__,par.myrank );
		// time dependend variables
    dimids[0] = tDim; dimids[1] = yDim; dimids[2] = xDim;
    ncwrap( ncmpi_def_var( ncid , "h" , SERGHEI_NC_REAL , 3 , dimids , &hVar  ) , __LINE__,par.myrank );
#if SERGHEI_WRITE_HZ
    ncwrap( ncmpi_def_var( ncid , "h+z" , SERGHEI_NC_REAL , 3 , dimids , &hzVar  ) , __LINE__,par.myrank );
#endif
    ncwrap( ncmpi_def_var( ncid , "u"      , SERGHEI_NC_REAL , 3 , dimids , &uVar  ) , __LINE__,par.myrank );
    ncwrap( ncmpi_def_var( ncid , "v"      , SERGHEI_NC_REAL , 3 , dimids , &vVar  ) , __LINE__,par.myrank );
    if(ss.inf.model){
      ncwrap( ncmpi_def_var( ncid , "inf" , SERGHEI_NC_REAL , 3 , dimids , &infVar  ) , __LINE__,par.myrank );
      ncwrap( ncmpi_def_var( ncid , "infVol" , SERGHEI_NC_REAL , 3 , dimids , &infVolVar  ) , __LINE__,par.myrank );
      if(ss.inf.model == INF_GREENAMPT){
        ncwrap( ncmpi_def_var( ncid , "infGaF" , SERGHEI_NC_REAL , 3 , dimids , &infGaFVar  ) , __LINE__,par.myrank );
      }
    }
    if(writeRain) ncwrap( ncmpi_def_var( ncid , "rain" , SERGHEI_NC_REAL , 3 , dimids , &rainVar  ) , __LINE__,par.myrank );

    #if SERGHEI_SUBSURFACE_MODEL
    ncwrap( ncmpi_def_var( ncid , "qss" , SERGHEI_NC_REAL , 3 , dimids , &qVar  ) , __LINE__ );
    #endif
	#if SERGHEI_DRAINAGE_MODEL
    ncwrap( ncmpi_def_var( ncid , "Dqss" , SERGHEI_NC_REAL , 3 , dimids , &DqVar  ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "Vqss" , SERGHEI_NC_REAL , 3 , dimids , &VqVar  ) , __LINE__ );
    #endif

		#if SERGHEI_MAXFLOOD > 0
		int nc_ndims = 2+SERGHEI_MAXFLOOD-1;
		if(nc_ndims==2){
			dimids[0] = yDim;
			dimids[1] = xDim;
		}
		ncwrap( ncmpi_def_var( ncid , "hMax"    , SERGHEI_NC_REAL , nc_ndims , dimids , &hMaxVar ) , __LINE__ , par.myrank);
		ncwrap( ncmpi_def_var( ncid , "momMax"    , SERGHEI_NC_REAL , nc_ndims , dimids , &momMaxVar ) , __LINE__ , par.myrank);
		ncwrap( ncmpi_def_var( ncid , "timehMax"    , SERGHEI_NC_REAL , nc_ndims , dimids , &timehMaxVar ) , __LINE__ , par.myrank);

		longname.assign("Maximum water depth");
		ncwrap( ncmpi_put_att_text(ncid, hMaxVar, "long_name", longname.length(),longname.c_str()),__LINE__);
		units.assign("m");
		ncwrap( ncmpi_put_att_text(ncid, hMaxVar, "units", units.length(),units.c_str()),__LINE__);

		longname.assign("Maximum momentum");
		ncwrap( ncmpi_put_att_text(ncid, momMaxVar, "long_name", longname.length(),longname.c_str()),__LINE__);
		units.assign("m^2/s");
		ncwrap( ncmpi_put_att_text(ncid, momMaxVar, "units", units.length(),units.c_str()),__LINE__);

		longname.assign("Time to maximum depth");
		ncwrap( ncmpi_put_att_text(ncid, timehMaxVar, "long_name", longname.length(),longname.c_str()),__LINE__);
		units.assign("s");
		ncwrap( ncmpi_put_att_text(ncid, timehMaxVar, "units", units.length(),units.c_str()),__LINE__);
		#endif

		#if SERGHEI_WRITE_SUBDOMS
			dimids[0] = yDim; dimids[1] = xDim;
			ncwrap( ncmpi_def_var( ncid , "subdom"    , NC_INT , 2 , dimids , &subdomVar ) , __LINE__ );
		#endif

		dimids[0] = yDim; dimids[1] = xDim;
    ncwrap( ncmpi_def_var( ncid , "z"    , SERGHEI_NC_REAL , 2 , dimids , &zVar ) , __LINE__,par.myrank );

		// define global attributes
		static char title[] = "SERGHEI simulation";
		ncwrap( ncmpi_put_att(ncid, NC_GLOBAL, "title", NC_CHAR,strlen(title)+1,title), __LINE__,par.myrank);

		std::string source = std::string("SERGHEI ");
        // std::string source = std::string("SERGHEI ") + SERGHEI_GIT_VERSION;
		ncwrap( ncmpi_put_att(ncid, NC_GLOBAL, "source", NC_CHAR,source.length()+1,source.c_str()), __LINE__,par.myrank);

		auto nowtime = std::chrono::system_clock::now();
		std::time_t now_time = std::chrono::system_clock::to_time_t(nowtime);
		source = std::string("Simulation started on ") + std::ctime(&now_time);
		ncwrap( ncmpi_put_att(ncid, NC_GLOBAL, "history", NC_CHAR,source.length()-1,source.c_str()), __LINE__,par.myrank);

		// define variable attributes
		// seems unnecessary if NAN is used
		/*
		#if SERGHEI_NC_ENABLE_NAN
		ncreal nan_value = SERGHEI_NAN;
  	ncwrap(ncmpi_put_att_real(ncid, zVar, SERGHEI_NC_FILL_VALUE_KEY, SERGHEI_NC_REAL, 1, &nan_value),__LINE__);
  	ncwrap(ncmpi_put_att_real(ncid, hVar, SERGHEI_NC_FILL_VALUE_KEY, SERGHEI_NC_REAL, 1, &nan_value),__LINE__);
  	ncwrap(ncmpi_put_att_real(ncid, uVar, SERGHEI_NC_FILL_VALUE_KEY, SERGHEI_NC_REAL, 1, &nan_value),__LINE__);
  	ncwrap(ncmpi_put_att_real(ncid, vVar, SERGHEI_NC_FILL_VALUE_KEY, SERGHEI_NC_REAL, 1, &nan_value),__LINE__);
		#endif
		*/

		longname.assign("projection_x_coordinate");
		ncwrap( ncmpi_put_att_text(ncid, xVar, "standard_name", longname.length(),longname.c_str()),__LINE__,par.myrank);
		longname.assign("x coordinate of projection");
		ncwrap( ncmpi_put_att_text(ncid, xVar, "long_name", longname.length(),longname.c_str()),__LINE__,par.myrank);
		units.assign("m");
		ncwrap( ncmpi_put_att_text(ncid, xVar, "units", units.length(),units.c_str()),__LINE__,par.myrank);

		longname.assign("projection_y_coordinate");
		ncwrap( ncmpi_put_att_text(ncid, yVar, "standard_name", longname.length(),longname.c_str()),__LINE__,par.myrank);
		longname.assign("y coordinate of projection");
		ncwrap( ncmpi_put_att_text(ncid, yVar, "long_name", longname.length(),longname.c_str()),__LINE__,par.myrank);
		units.assign("m");
		ncwrap( ncmpi_put_att_text(ncid, yVar, "units", units.length(),units.c_str()),__LINE__,par.myrank);

		longname.assign("Bed elevation");
		ncwrap( ncmpi_put_att_text(ncid, zVar, "long_name", longname.length(),longname.c_str()),__LINE__,par.myrank);
		units.assign("m");
		ncwrap( ncmpi_put_att_text(ncid, zVar, "units", units.length(),units.c_str()),__LINE__,par.myrank);

		longname.assign("Water depth");
		ncwrap( ncmpi_put_att_text(ncid, hVar, "long_name", longname.length(),longname.c_str()),__LINE__,par.myrank);
		units.assign("m");
		ncwrap( ncmpi_put_att_text(ncid, hVar, "units", units.length(),units.c_str()),__LINE__,par.myrank);

		#if SERGHEI_WRITE_HZ
		longname.assign("Water elevation");
		ncwrap( ncmpi_put_att_text(ncid, hzVar, "long_name", longname.length(),longname.c_str()),__LINE__,par.myrank);
		units.assign("m");
		ncwrap( ncmpi_put_att_text(ncid, hzVar, "units", units.length(),units.c_str()),__LINE__,par.myrank);
		#endif

		longname.assign("Velocity x-component");
		ncwrap( ncmpi_put_att_text(ncid, uVar, "long_name", longname.length(),longname.c_str()),__LINE__,par.myrank);
		units.assign("m/s");
		ncwrap( ncmpi_put_att_text(ncid, uVar, "units", units.length(),units.c_str()),__LINE__,par.myrank);

		longname.assign("Velocity y-component");
		ncwrap( ncmpi_put_att_text(ncid, vVar, "long_name", longname.length(),longname.c_str()),__LINE__,par.myrank);
		units.assign("m/s");
		ncwrap( ncmpi_put_att_text(ncid, vVar, "units", units.length(),units.c_str()),__LINE__,par.myrank);

    if(ss.inf.model){
			longname.assign("Infiltration rate");
			ncwrap( ncmpi_put_att_text(ncid, infVar, "long_name", longname.length(),longname.c_str()),__LINE__,par.myrank);
			units.assign("m/s");
			ncwrap( ncmpi_put_att_text(ncid, infVar, "units", units.length(),units.c_str()),__LINE__,par.myrank);

			longname.assign("Accumulated infiltration volume");
			ncwrap( ncmpi_put_att_text(ncid, infVolVar, "long_name", longname.length(),longname.c_str()),__LINE__,par.myrank);
			units.assign("m^3");
			ncwrap( ncmpi_put_att_text(ncid, infVolVar, "units", units.length(),units.c_str()),__LINE__,par.myrank);
      if(ss.inf.model == INF_GREENAMPT){
        longname.assign("Green-Ampt cumulative infiltrated depth per cell");
        ncwrap( ncmpi_put_att_text(ncid, infGaFVar, "long_name", longname.length(),longname.c_str()),__LINE__,par.myrank);
        units.assign("m");
        ncwrap( ncmpi_put_att_text(ncid, infGaFVar, "units", units.length(),units.c_str()),__LINE__,par.myrank);
      }
		}
    if(writeRain){
			longname.assign("Rainfall rate");
			ncwrap( ncmpi_put_att_text(ncid, rainVar, "long_name", longname.length(),longname.c_str()),__LINE__,par.myrank);
			units.assign("m/s");
			ncwrap( ncmpi_put_att_text(ncid, rainVar, "units", units.length(),units.c_str()),__LINE__,par.myrank);
    }

		#if SERGHEI_WRITE_SUBDOMS
		longname.assign("Subdomain index");
		ncwrap( ncmpi_put_att_text(ncid, subdomVar , "long_name", longname.length(),longname.c_str()),__LINE__,par.myrank);
		units.assign("NA");
		ncwrap( ncmpi_put_att_text(ncid, subdomVar, "units", units.length(),units.c_str()),__LINE__,par.myrank);
		static int range[] = {0,par.nranks-1};
		ncwrap( ncmpi_put_att_int(ncid, subdomVar, "valid_range", NC_INT, 2, range),__LINE__,par.myrank);
		#endif


    // End "define" mode
    ncwrap( ncmpi_enddef( ncid ) , __LINE__,par.myrank );

		#if SERGHEI_DEBUG_OUTPUT
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << " NetCDF define mode completed" << std::endl; ;
		#endif

    // Compute x, y coordinates
    Kokkos::parallel_for("compute_grid_coord_x", dom.nx , KOKKOS_LAMBDA(int i) {
      xCoord(i) = dom.xll + ( par.i_beg + i + 0.5) * dom.dxConst;
    });
    Kokkos::parallel_for("compute_grid_coord_x", dom.ny , KOKKOS_LAMBDA(int j) {
      yCoord(j) = dom.yll + dom.ny_glob*dom.dxConst - ( par.j_beg + j + 0.5) * dom.dxConst;
    });
    Kokkos::fence();

    // Write out x, y coordinates
    st[0] = par.i_beg;
    ct[0] = dom.nx;
    ncwrap( ncmpi_put_vara_double_all( ncid , xVar , st , ct , xCoord.data() ) , __LINE__,par.myrank );
    st[0] = par.j_beg;
    ct[0] = dom.ny;
    ncwrap( ncmpi_put_vara_double_all( ncid , yVar , st , ct , yCoord.data() ) , __LINE__,par.myrank );

    st[0] = par.j_beg; st[1] = par.i_beg;
    ct[0] = dom.ny   ; ct[1] = dom.nx   ;

		#if SERGHEI_NC_ENABLE_MISSING_VALUE
		real missing_value = dom.MISSING_VALUE;
  	ncwrap(nc_put_att_float(ncid, zVar, SERGHEI_NC_MISSING_VALUE, SERGHEI_NC_REAL, 1, &missing_value));
		#endif
		#if SERGHEI_DEBUG_OUTPUT
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << " NetCDF coordinates written" << std::endl; ;
		#endif

		// write elevation
		writeNetCDFfield(dom,ncid,zVar,st,ct,state.isnodata,state.z,data);
		#if SERGHEI_DEBUG_OUTPUT
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << " NetCDF z written" << std::endl; ;
		#endif

		#if SERGHEI_WRITE_SUBDOMS
			dimids[0] = yDim; dimids[1] = xDim;
			realArr subdom  = realArr("subdom",dom.nCellMem);
			Kokkos::parallel_for("subdom",dom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
				subdom(iGlob) = dom.id;
			});
			Kokkos::fence();
			writeNetCDFfield(dom,ncid,subdomVar,st,ct,state.isnodata,subdom,data);
			#if SERGHEI_DEBUG_OUTPUT
				std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << " NetCDF subdom written" << std::endl; ;
			#endif
		#endif


    writeStateNETCDF(state, dom, ss, par);

    ncwrap( ncmpi_close(ncid) , __LINE__,par.myrank );


  }

/////////////////////////////////////////////////////////
	// NetCDF writer - requires the initialisation
  void outputNETCDF(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par, std::string dir) {

		std::string filename;
		filename=dir+"output.nc";




    // Create the file
    ncwrap( ncmpi_open( MPI_COMM_WORLD , filename.c_str() , NC_WRITE , MPI_INFO_NULL , &ncid ) , __LINE__,par.myrank );
    ncwrap( ncmpi_inq_varid( ncid , "h" , &hVar  ) , __LINE__,par.myrank );
#if SERGHEI_WRITE_HZ
    ncwrap( ncmpi_inq_varid( ncid , "h+z" , &hzVar  ) , __LINE__,par.myrank );
#endif
    ncwrap( ncmpi_inq_varid( ncid , "u"      , &uVar  ) , __LINE__,par.myrank );
    ncwrap( ncmpi_inq_varid( ncid , "v"      , &vVar  ) , __LINE__,par.myrank );
    if(ss.inf.model){
      ncwrap( ncmpi_inq_varid( ncid , "inf"      , &infVar  ) , __LINE__,par.myrank );
      ncwrap( ncmpi_inq_varid( ncid , "infVol"      , &infVolVar  ) , __LINE__,par.myrank );
      if(ss.inf.model == INF_GREENAMPT){
        ncwrap( ncmpi_inq_varid( ncid , "infGaF"      , &infGaFVar  ) , __LINE__,par.myrank );
      }
    }
    if(writeRain) ncwrap( ncmpi_inq_varid( ncid , "rain"      , &rainVar  ) , __LINE__,par.myrank );
    #if SERGHEI_SUBSURFACE_MODEL
    ncwrap( ncmpi_inq_varid( ncid , "qss" , &qVar  ) , __LINE__ );
    #endif
	#if SERGHEI_DRAINAGE_MODEL
    ncwrap( ncmpi_inq_varid( ncid , "Dqss" , &DqVar  ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "Vqss" , &VqVar  ) , __LINE__ );
    #endif

		#if SERGHEI_MAXFLOOD > 0
    	ncwrap( ncmpi_inq_varid( ncid , "hMax" , &hMaxVar  ) , __LINE__ );
    	ncwrap( ncmpi_inq_varid( ncid , "momMax" , &momMaxVar ) , __LINE__ );
    	ncwrap( ncmpi_inq_varid( ncid , "timehMax" , &timehMaxVar  ) , __LINE__ );
		#endif

	  writeStateNETCDF(state, dom, ss, par);

    ncwrap( ncmpi_close(ncid) , __LINE__,par.myrank );

  }

	// Builds the NetCDF dataset for the state
  void writeStateNETCDF(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par) {
		#if SERGHEI_DEBUG_OUTPUT
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl; ;
		#endif
  	ncArr data = ncArr("data",dom.nCell);

    MPI_Offset st[3], ct[3];
	 	double timeIter[numOut+1];

	 	//write t. As the first one is written in the first iteration we should add +1
	 	for (int i=0; i<numOut+1; i++) { timeIter[i] = i*1.0*outFreq + dom.startTime;}

    st[0] = 0;
	 	ct[0] = numOut+1;

    ncwrap( ncmpi_put_vara_double_all( ncid , tVar ,  st , ct , timeIter ) , __LINE__,par.myrank );

    st[0] = numOut; st[1] = par.j_beg; st[2] = par.i_beg;
    ct[0] = 1     ; ct[1] = dom.ny   ; ct[2] = dom.nx   ;

    // Write out depth
		writeNetCDFfield(dom,ncid,hVar,st,ct,state.isnodata,state.h,data);

#if SERGHEI_WRITE_HZ
    Kokkos::parallel_for("ncwrap_h+z", dom.nCell , KOKKOS_LAMBDA(int iGlob) {
	 	 int i,j;
		dom.unpackIndices(iGlob,j,i);
		int ii=(hc+j)*(dom.nx+2*hc)+hc+i;//index with the extended domain (including halo cells)
      data(iGlob) = state.h(ii)+state.z(ii);
			#if SERGHEI_NC_ENABLE_NAN
			if(state.isnodata(ii)) data(iGlob) = SERGHEI_NAN;
			#endif
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_real_all( ncid , hzVar , st , ct , data.data() ) , __LINE__,par.myrank );
#endif

    #if SERGHEI_SUBSURFACE_MODEL
    Kokkos::parallel_for("ncwrap_qss", dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
	 	 int i,j;
		dom.unpackIndices(iGlob,j,i);
		int ii=(hc+j)*(dom.nx+2*hc)+hc+i;//index with the extended domain (including halo cells)
      // data(iGlob) = state.qss(ii);
      data(iGlob) = state.qss(iGlob);
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_real_all( ncid , qVar , st , ct , data.data() ) , __LINE__ );
    #endif

    // Write out x-velocity

    Kokkos::parallel_for("ncwrap_u", dom.nCell , KOKKOS_LAMBDA(int iGlob) {
			int ii = dom.getIndex(iGlob);
		if(state.h(ii)>TOL12){
     	data(iGlob) = state.hu(ii)/state.h(ii);
		}else{
			data(iGlob)=0.0;
		}
		#if SERGHEI_NC_ENABLE_NAN
		if(state.isnodata(ii)){
			data(iGlob) = SERGHEI_NAN;
		}
		#endif
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_real_all( ncid , uVar , st , ct , data.data() ) , __LINE__,par.myrank );

    // Write out y-velocity
	 Kokkos::parallel_for("ncwrap_v", dom.nCell , KOKKOS_LAMBDA(int iGlob) {
		int ii = dom.getIndex(iGlob);
		if(state.h(ii)>TOL12){
  	   	data(iGlob) = state.hv(ii)/state.h(ii);
		}else{
			data(iGlob)=0.0;
		}
		#if SERGHEI_NC_ENABLE_NAN
		if(state.isnodata(ii)){
			data(iGlob) = SERGHEI_NAN;
		}
		#endif

    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_real_all( ncid , vVar , st , ct , data.data() ) , __LINE__,par.myrank );

    if(ss.inf.model){
      // infiltration rates
			writeNetCDFfield(dom,ncid,infVar,st,ct,state.isnodata,ss.inf.rate,data);
      // accumulated infiltration
			writeNetCDFfield(dom,ncid,infVolVar,st,ct,state.isnodata,ss.inf.infVol,data);
      if(ss.inf.model == INF_GREENAMPT){
        writeNetCDFfield(dom,ncid,infGaFVar,st,ct,state.isnodata,ss.inf.gaF,data);
      }
    }

    if(writeRain) writeNetCDFfield(dom,ncid,rainVar,st,ct,state.isnodata,ss.rainRate,data);

			#if SERGHEI_MAXFLOOD > 0
				#if SERGHEI_MAXFLOOD == 1		// write only at the end, otherwise writes at every time step
				if(numOut == nOut){
    			st[0] = par.j_beg; st[1] = par.i_beg;
    			ct[0] = dom.ny   ; ct[1] = dom.nx   ;
				#endif
				writeNetCDFfield(dom,ncid,hMaxVar,st,ct,state.isnodata,state.hMax,data);
				writeNetCDFfield(dom,ncid,momMaxVar,st,ct,state.isnodata,state.momentumMax,data);
				writeNetCDFfield(dom,ncid,timehMaxVar,st,ct,state.isnodata,state.time_hMax,data);
			#if SERGHEI_MAXFLOOD == 1
			}
			#endif
		#endif
}

/////////////////////////////////////////

  //Error reporting routine for the PNetCDF I/O
  void ncwrap( int ierr , int line, int rank) {
    if (ierr != NC_NOERR) {
      std::cerr<< RERROR "NetCDF reports error from rank " << rank << " at " << __FILE__ << ":" << line << std::endl << RERROR << ncmpi_strerror(ierr) << std::endl;
      exit(-1);
    }
  }
  void ncwrap( int ierr , int line) {
    if (ierr != NC_NOERR) {
      std::cerr<< RERROR "NetCDF reports error at " << __FILE__ << ":" << line << std::endl << RERROR << ncmpi_strerror(ierr) << std::endl;
      exit(-1);
    }
  }

#if SERGHEI_LPT == 1

  	void outputParticle(Domain const &dom, ParticleTracker const &particles, std::string dir){

      std::string filename;
      filename = dir+"particles"+std::to_string(numOut)+".vtk";

      std::ofstream particleFile(filename);

      //particleFile.open(filename);
      if (particleFile.is_open()){
        particleFile << "# vtk DataFile Version 2.0\n";
        particleFile << "Particles\n";
        particleFile << "ASCII\n";
        particleFile << "DATASET UNSTRUCTURED_GRID\n";

        int num_points;
        int noactive = 0;
        num_points = particles.N_par;
        for (int i = 0; i < num_points; i++){
          if(particles.particles[i].active==0){
              noactive+=1;
          }
        }
        particleFile << "POINTS " << num_points-noactive << " float"<<std::endl;
        for (int i = 0; i < num_points; i++) {
          if(particles.particles[i].active==1){
            particleFile << std::setprecision(9) << particles.particles[i].x(0) << " " <<  particles.particles[i].x(1) << " " << particles.particles[i].z <<std::endl;
          }
        }

        //Write information about particles as cells
        int num_cells = num_points-noactive;
        particleFile << "CELLS " << num_cells << " " << 2*num_cells <<std::endl;
        for (int i = 0; i < num_cells; i++) {
          particleFile << "1 " << i <<std::endl;
        }

        //Write information about the type of cell (in this case, points)
        particleFile << "CELL_TYPES " << num_cells <<std::endl;
        for (int i = 0; i < num_cells; i++){
          particleFile << "1"<<std::endl;
        }
        particleFile.close();
     }
  }

    void outputIniParticle(Domain const &dom, ParticleTracker const &particles, std::string dir){

    std::string filename;
    int numOutIni;
    numOutIni=0;
    filename = dir+"particles"+std::to_string(numOutIni)+".vtk";

    std::ofstream particleFile(filename);

    //particleFile.open(filename);
    if (particleFile.is_open()){
      particleFile << "# vtk DataFile Version 2.0\n";
      particleFile << "Particles\n";
      particleFile << "ASCII\n";
      particleFile << "DATASET UNSTRUCTURED_GRID\n";

      int num_points;
      num_points = particles.N_par;
      int noactive = 0;
      num_points = particles.N_par;
      for (int i = 0; i < num_points; i++){
        if(particles.particles[i].active==0){
            noactive+=1;
        }
      }
      particleFile << "POINTS " << num_points-noactive << " float"<<std::endl;
      for (int i = 0; i < num_points; i++){
		if(particles.particles[i].active==1){
            particleFile << std::setprecision(9) << particles.particles[i].x(0) << " " << particles.particles[i].x(1) << " " << particles.particles[i].z <<std::endl;
		}
	  }


      //Write information about particles as cells
      int num_cells = num_points;
      particleFile << "CELLS " << num_cells-noactive << " " << 2*(num_cells-noactive) <<std::endl;
      for (int i = 0; i < num_cells-noactive; i++) {
        particleFile << "1 " << i <<std::endl;
      }

      //Write information about the type of cell (in this case, points)
      particleFile << "CELL_TYPES " << num_cells-noactive <<std::endl;
      for (int i = 0; i < num_cells-noactive; i++){
        particleFile << "1"<<std::endl;
      }
      particleFile.close();
     }
  }

	void writeParticleFile(Domain const &dom, ParticleTracker const &particles, Parallel const &par, std::string dir){
		std::string filename = dir + "particles.out";
		particlesFile.open(filename);
		//if (particlesFile.is_open()){
		particlesFile << "ID Particle \t CoveredDistance \t TotalTime \t TravelTime \t StopTime" << std::endl;
		int num_points;
		num_points = particles.N_par;
		for(int i = 0; i < num_points; i++){
			particlesFile << std::setprecision(6) << i << " " << particles.particles[i].CoveredDistance << " " << particles.particles[i].TotalLife << " " << particles.particles[i].TravelTime << " " << particles.particles[i].StopTime << std::endl;
		}
		
		particlesFile.close();
		//}
    }

	void writeParticle_TemporalFile(Domain const &dom, ParticleTracker const &particles, std::string dir){
		std::string filename = dir + "particles_temporal.out";
		particlesFile.open(filename, std::ios::app); 
		//particlesFile.open(filename);
		if (particlesFile.is_open()){
			//particlesFile << "ID Particle \t CoveredDistance \t TotalTime \t TravelTime \t StopTime" << std::endl;
			int num_points;
			real total_distance;
			num_points = particles.N_par;
			total_distance=0.0;
			for(int i = 0; i < num_points; i++){
				total_distance+=particles.particles[i].CoveredDistance;
			}
			particlesFile << dom.etime << " " << total_distance << " " << std::endl;
			particlesFile.close();
		}
    }

	void outputInitParticlesNETCDF(ParticleTracker const &particles, Domain const &dom, Parallel const &par, std::string dir) {

		int dimids[2];
		MPI_Offset st[3], ct[3];
		ncArr data   = ncArr("data",particles.N_par);
		static char timeUnits[] = "seconds" ;
		std::string filename;
		std::string longname;
		std::string units;

		lpt_filename = dir + "lpt.nc";

		// Create the file
		ncwrap( ncmpi_create( MPI_COMM_WORLD , lpt_filename.c_str() , SERGHEI_NC_MODE , MPI_INFO_NULL , &ncidP ) , __LINE__ );

		// Create the dimensions
		ncwrap( ncmpi_def_dim( ncidP , "t" , (MPI_Offset) NC_UNLIMITED , &tpDim ) , __LINE__ );
		ncwrap( ncmpi_def_dim( ncidP , "p" , (MPI_Offset) particles.N_par  , &pDim ) , __LINE__ );
		// Create the variables for dimensions
		dimids[0] = tDim;
		ncwrap( ncmpi_def_var( ncidP , "t"      , NC_DOUBLE , 1 , dimids , &tpVar ) , __LINE__ );
		ncwrap( ncmpi_put_att_text (ncidP, tVar, "units",strlen(timeUnits), timeUnits), __LINE__ );
		dimids[0] = pDim;
		ncwrap( ncmpi_def_var( ncidP , "p"      , NC_INT , 1 , dimids , &pVar ) , __LINE__ );
		ncwrap( ncmpi_put_att_text (ncidP, pVar, "units",strlen(timeUnits), timeUnits), __LINE__ );
		// time dependent variables
		dimids[0] = tDim; dimids[1] = pDim;
		ncwrap( ncmpi_def_var( ncidP , "x"      , SERGHEI_NC_REAL , 2 , dimids , &xpVar ) , __LINE__ );
		ncwrap( ncmpi_def_var( ncidP , "y"      , SERGHEI_NC_REAL , 2 , dimids , &ypVar ) , __LINE__ );
		//ncwrap( ncmpi_def_var( ncidP , "z" , SERGHEI_NC_REAL , 1 , dimids , &zpVar  ) , __LINE__ );


		// define global attributes
		static char title[] = "SERGHEI-LPT simulation";
		ncwrap( ncmpi_put_att(ncidP, NC_GLOBAL, "title", NC_CHAR,strlen(title)+1,title), __LINE__);

		longname.assign("x-coordinate");
		ncwrap( ncmpi_put_att_text(ncidP, xpVar, "long_name", longname.length(),longname.c_str()),__LINE__);
		units.assign("m");
		ncwrap( ncmpi_put_att_text(ncidP, xpVar, "units", units.length(),units.c_str()),__LINE__);

		longname.assign("y-coordinate");
		ncwrap( ncmpi_put_att_text(ncidP, ypVar, "long_name", longname.length(),longname.c_str()),__LINE__);
		units.assign("m");
		ncwrap( ncmpi_put_att_text(ncidP, ypVar, "units", units.length(),units.c_str()),__LINE__);

		//longname.assign("z-coordinate");
		//ncwrap( ncmpi_put_att_text(ncidP, zpVar, "long_name", longname.length(),longname.c_str()),__LINE__);
		//units.assign("m");
		//ncwrap( ncmpi_put_att_text(ncidP, zpVar, "units", units.length(),units.c_str()),__LINE__);


		// End "define" mode
		ncwrap( ncmpi_enddef( ncidP ) , __LINE__ );

		#if SERGHEI_DEBUG_OUTPUT
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << " NetCDF define mode completed" << std::endl; ;
		#endif

		st[0] = 0; 	// DCV REVIEW: this needs revision
		ct[0] = particles.N_par;			// DCV REVIEW: this needs revision

			intArr pID = intArr("particleID",particles.N_par);

			Kokkos::parallel_for("generate_particle_id", particles.N_par , KOKKOS_LAMBDA(int ii) {
		pID(ii) = ii;
		});
		ncwrap( ncmpi_put_vara_int_all( ncidP , pVar , st , ct , pID.data() ) , __LINE__ );

		writeStateNETCDF_particles(particles, dom, par);

		ncwrap( ncmpi_close(ncidP) , __LINE__ );

	}

	void writeStateNETCDF_particles(ParticleTracker const &particles, Domain const &dom, Parallel const &par) {
		#if SERGHEI_DEBUG_OUTPUT
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl; ;
		#endif


		MPI_Offset st[2], ct[2];
		double timeIter[numOut+1];

		//write t. As the first one is written in the first iteration we should add +1
		for (int i=0; i<numOut+1; i++) { timeIter[i] = i*1.0*outFreq + dom.startTime;}
		st[0] = 0;
		ct[0] = numOut+1;

		ncwrap( ncmpi_put_vara_double_all( ncidP , tpVar ,  st , ct , timeIter ) , __LINE__ );

		st[0] = numOut; st[1] = 0;
		ct[0] = 1     ; ct[1] = particles.N_par;

		// DCV REVIEW: why not do this when particles.particles[i].active is computed?
		// DCV REVIEW: implemented parallel_reduce
		// Calculate particles no active
		int inactive=0;
		//	for(int i=0; i<particles.N_par;i++){
		Kokkos::parallel_reduce("particle_count_active", particles.N_par, KOKKOS_LAMBDA(int ii, int inactive_l){
			if(particles.particles(ii).active!=1) inactive_l++;
		}, Kokkos::Sum<int>(inactive) );
		// DCV REVIEW: do we need MPI reduction?

		int Np;
		Np = particles.N_par;
		//Np = particles.N_par - inactive;
		ncArr data = ncArr("data",Np);	// DCV REVIEW size was wrong. Should be number of particles. Do we really need the data buffer?




		// Write out x-coordinate
		Kokkos::parallel_for("ncwrap_x_par",Np, KOKKOS_LAMBDA(int ii) {
			//if(particles.particles(ip).active)
			data(ii) = particles.particles(ii).x(0);
		});
		Kokkos::fence();

		ncwrap( ncmpi_put_vara_real_all( ncidP , xpVar , st , ct , data.data() ) , __LINE__ );

		// Write out y-coordinate
		Kokkos::parallel_for("ncwrap_y_par", Np , KOKKOS_LAMBDA(int ii) {
			//if(particles.particles[iGlob].active==1)
			data(ii) = particles.particles(ii).x(1);
		});
		Kokkos::fence();

		ncwrap( ncmpi_put_vara_real_all( ncidP , ypVar , st , ct , data.data() ) , __LINE__ );

		// DCV REVIEW: not sure we need this
		// Write out z-coordinate
		
		//Kokkos::parallel_for("ncwrap_z_par", Np , KOKKOS_LAMBDA(int ii) {
			//if(particles.particles[iGlob].active==1)
			//data(ii) = particles.particles(ii).z;
		//});
		//Kokkos::fence();
		//ncwrap( ncmpi_put_vara_real_all( ncidP , zpVar , st , ct , data.data() ) , __LINE__ );
		
	}
		
	void outputNETCDF_particles(ParticleTracker const &particles, Domain const &dom, Parallel const &par, std::string dir) {
		lpt_filename=dir+"lpt.nc";
		// Create the file
		ncwrap( ncmpi_open( MPI_COMM_WORLD , lpt_filename.c_str() , NC_WRITE , MPI_INFO_NULL , &ncidP ) , __LINE__ );
		ncwrap( ncmpi_inq_varid( ncidP , "x" , &xpVar  ) , __LINE__ );
		ncwrap( ncmpi_inq_varid( ncidP , "y" , &ypVar  ) , __LINE__ );
		//ncwrap( ncmpi_inq_varid( ncidP , "z" , &zpVar  ) , __LINE__ );
		writeStateNETCDF_particles(particles, dom, par);
		ncwrap( ncmpi_close(ncidP) , __LINE__ );

	}
#endif
	// Writes VTK file
	void outputVTK(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par, std::string dir){

		std::string filename;
		//filename = dir+"result_."+std::to_string(par.myrank)+".vtk";
		filename = dir+"result"+std::to_string(numOut)+".vtk";
	 	real *xCoord;
  	real *yCoord;
		real *data_cpu;

		int ncells=dom.nCell;

		int nVars=4; //z, h, (u,v)
		#if SERGHEI_WRITE_HZ
		nVars++;
		#endif
		#if SERGHEI_DEBUG_BOUNDARY
		nVars++;
		#endif
		#if SERGHEI_LPT
		nVars++;
		#endif
		#if SERGHEI_VERTICAL_VELOCITY
		nVars++;
		#endif
    if(ss.inf.model){
      nVars = nVars+2;  // inf, infVol
      if(ss.inf.model == INF_GREENAMPT) nVars = nVars+1;  // infGaF (cumulative depth F)
    }

		realArr data  = realArr("data",nVars*ncells);
		#ifdef __NVCC__
			cudaMallocHost( &data_cpu , nVars*ncells*sizeof(real) );
		#else
			data_cpu = data.data();
		#endif

		xCoord=(real*) malloc((dom.nx+1)*sizeof(real));
		yCoord=(real*) malloc((dom.ny+1)*sizeof(real));

		#if SERGHEI_DEBUG_OUTPUT
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << " VTK data allocated" << std::endl; ;
		#endif

    	// Compute x, y coordinates
    	for (int i=0; i<dom.nx+1; i++) {
   		xCoord[i] = dom.xll + ( par.i_beg + i) * dom.dxConst;
    	};

	 	for (int j=0; j<dom.ny+1; j++) {
      	yCoord[j] = dom.yll + dom.ny_glob*dom.dxConst - ( par.j_beg + j) * dom.dxConst;
    	};

		#if SERGHEI_DEBUG_OUTPUT
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << " VTK coordinate data generated" << std::endl; ;
		#endif
		int nnodes=(dom.nx+1)*(dom.ny+1);
		int of_sw = 3;
		int offset_z = 0;
		int offset_h = ncells;
		int offset_hu = 2*ncells;
		int offset_hv = 3*ncells;
		#if SERGHEI_LPT && SERGHEI_VERTICAL_VELOCITY==0
		int offset_particles = 4 * ncells;
		#elif SERGHEI_LPT && SERGHEI_VERTICAL_VELOCITY
		int offset_particles = 4 * ncells; 
		int offset_vertical_velocity = 5 * ncells;
		#elif SERGHEI_LPT==0 && SERGHEI_VERTICAL_VELOCITY==1
		int offset_vertical_velocity = 4 * ncells;
		#endif
		int of_bc = of_sw;
		#if SERGHEI_DEBUG_BOUNDARY
		of_bc++;
		int offset_bc = of_bc*ncells;
		#endif
		int of_inf = of_bc + 1;
		int offset_infrate = of_inf * ncells;
		int offset_infVol = (of_inf+1) * ncells;
		int offset_gaF = (of_inf+2) * ncells;
		#if SERGHEI_DEBUG_OUTPUT
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << " VTK offsets generated" << std::endl;
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << " offset_z" << std::endl;
		#endif

    	Kokkos::parallel_for("vtkwrap_all", ncells , KOKKOS_LAMBDA(int iGlob) {
			int ii = dom.getIndex(iGlob);

			data(iGlob+offset_z) = state.z(ii);
			data(iGlob+offset_h) = state.h(ii);
			data(iGlob+offset_hu) = state.hu(ii);
			data(iGlob+offset_hv) = state.hv(ii);
			if(ss.inf.model){
				data(iGlob+offset_infrate) = ss.inf.rate(ii);
				data(iGlob+offset_infVol) = ss.inf.infVol(ii);
				if(ss.inf.model == INF_GREENAMPT){
					data(iGlob+offset_gaF) = ss.inf.gaF(ii);
				}
			}
			#if SERGHEI_LPT && SERGHEI_VERTICAL_VELOCITY
			if(dom.etime==0.0){
			  data(iGlob + offset_particles) = state.tparticles(ii);
			}else{
			  data(iGlob + offset_particles) = state.tparticles(ii)/(dom.etime);
			}
			data(iGlob+offset_vertical_velocity) = state.w(ii);
			#elif SERGHEI_LPT==0 && SERGHEI_VERTICAL_VELOCITY==1
			data(iGlob+offset_vertical_velocity) = state.w(ii);
			#endif
			#if SERGHEI_DEBUG_BOUNDARY
			data(iGlob + offset_bc) = state.isBound(ii);
			#endif
      // WARNING if you implement a new variable, you have to handle the offsets in a general case (yes, you!), to handle the possibility of different variable combinations
    	});
    	Kokkos::fence();

		#if SERGHEI_DEBUG_OUTPUT
			std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << " VTK data wrappers set" << std::endl; ;
		#endif

		#ifdef __NVCC__
			cudaMemcpyAsync( data_cpu , data.data() , nVars*ncells*sizeof(real) , cudaMemcpyDeviceToHost );
			cudaDeviceSynchronize();
		#endif

		int i,j,iGlob;

		std::ofstream fOutStream(filename);
		if (fOutStream.is_open()){
			//fOutStream << "# vtk DataFile Version 3.0.\nOutput file " << filename <<"\nASCII\nDATASET UNSTRUCTURED_GRID\n";
			fOutStream << "# vtk DataFile Version 3.0.\nOutputfile\nASCII\nDATASET UNSTRUCTURED_GRID\n";

			fOutStream << "POINTS " << nnodes  << SERGHEI_VTK_REAL << "\n";
			for(j=0;j<=dom.ny;j++){
				for(i=0;i<=dom.nx;i++){
					fOutStream << std::setprecision(9) << xCoord[i] << " " << yCoord[j]  << " 0.0\n";
				}
			}
			fOutStream << "CELLS " << ncells << " " << 5*ncells <<"\n";
			for(iGlob=0;iGlob<ncells;iGlob++){
			   int ii,jj;
				dom.unpackIndices(iGlob,jj,ii);
				fOutStream << "4 " << jj*(dom.nx+1)+ii << " " << (jj+1)*(dom.nx+1)+ii << " " << (jj+1)*(dom.nx+1)+ii+1 << " " << jj*(dom.nx+1)+ii+1 << "\n";
			}

			fOutStream << "CELL_TYPES " << ncells << "\n";
			for(i=0;i<ncells;i++){
				fOutStream << "9\n";
			}

			fOutStream << "CELL_DATA " << ncells << "\n";
			fOutStream << "SCALARS z " << SERGHEI_VTK_REAL << "\n";
			fOutStream << "LOOKUP_TABLE default\n";
			for(iGlob=0;iGlob<ncells;iGlob++){
				fOutStream << std::setprecision(9) << data_cpu[iGlob] << "\n";
			}

			fOutStream << "SCALARS h " << SERGHEI_VTK_REAL <<"\n";
			fOutStream << "LOOKUP_TABLE default\n";
			for(iGlob=0;iGlob<ncells;iGlob++){
				fOutStream << std::setprecision(9) << data_cpu[iGlob+offset_h]  << "\n";
			}

#if SERGHEI_WRITE_HZ
			fOutStream << "SCALARS h+z "<< SERGHEI_VTK_REAL << "\n";
			fOutStream << "LOOKUP_TABLE default\n";
			for(iGlob=0;iGlob<ncells;iGlob++){
				fOutStream << std::setprecision(9) << data_cpu[iGlob] + data_cpu[iGlob+offset_h] << "\n";
			}
#endif

			fOutStream << "VECTORS velocity " << SERGHEI_VTK_REAL << "\n";
			for(iGlob=0;iGlob<ncells;iGlob++){
				if(data_cpu[iGlob+ncells]>TOL12){
					fOutStream << std::setprecision(9) << data_cpu[iGlob+offset_hu]/data_cpu[iGlob+offset_h] << " " << data_cpu[iGlob+offset_hv]/data_cpu[iGlob+offset_h] << " 0.0\n";
				}else{
					fOutStream << "0.0 0.0 0.0\n" ;
				}
			}

      if(ss.inf.model){
        fOutStream << "SCALARS inf " << SERGHEI_VTK_REAL << "\n";
        fOutStream << "LOOKUP_TABLE default\n";
        for(iGlob=0;iGlob<ncells;iGlob++){
          fOutStream << std::setprecision(9) << data_cpu[iGlob+offset_infrate]  << "\n";
        }

        fOutStream << "SCALARS infVol " << SERGHEI_VTK_REAL << "\n";
        fOutStream << "LOOKUP_TABLE default\n";
        for(iGlob=0;iGlob<ncells;iGlob++){
          fOutStream << std::setprecision(9) << data_cpu[iGlob+offset_infVol]  << "\n";
        }
        if(ss.inf.model == INF_GREENAMPT){
          fOutStream << "SCALARS infGaF " << SERGHEI_VTK_REAL << "\n";
          fOutStream << "LOOKUP_TABLE default\n";
          for(iGlob=0;iGlob<ncells;iGlob++){
            fOutStream << std::setprecision(9) << data_cpu[iGlob+offset_gaF]  << "\n";
          }
        }
      }
	  #if SERGHEI_LPT
		fOutStream << "SCALARS t_particles " << SERGHEI_VTK_REAL << "\n";
		fOutStream << "LOOKUP_TABLE default\n";
		for(iGlob=0;iGlob<ncells;iGlob++){
			fOutStream << std::setprecision(9) << data_cpu[iGlob+offset_particles]  << "\n";
		}
	  #endif
	  #if SERGHEI_VERTICAL_VELOCITY
		fOutStream << "SCALARS w_velocity " << SERGHEI_VTK_REAL << "\n";
		fOutStream << "LOOKUP_TABLE default\n";
		for(iGlob=0;iGlob<ncells;iGlob++){
			fOutStream << std::setprecision(9) << data_cpu[iGlob+offset_vertical_velocity]  << "\n";
		}
	  #endif

			#if SERGHEI_DEBUG_BOUNDARY
        fOutStream << "SCALARS bc int\n";
        fOutStream << "LOOKUP_TABLE default\n";
        for(iGlob=0;iGlob<ncells;iGlob++){
          fOutStream << data_cpu[iGlob+offset_bc]  << "\n";
        }

			#endif
      fOutStream.close();
		}

	}






	/*void initBIN(Domain const &dom, Parallel const &par){

		cur_proc_data_size = dom.nx*dom.ny;

		if (par.myrank == 0)
			recvcounts = new int[par.nranks];
			MPI_Gather(&cur_proc_data_size, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

		if (par.myrank == 0){
			displs = new int[par.nranks];
			displs[0] = 0;
			total_data_size += recvcounts[0];

			for (int i = 1; i < par.nranks; i++){
				total_data_size += recvcounts[i];
				displs[i] = displs[i-1] + recvcounts[i-1];
			}
			total_data_arr = new real[total_data_size];
		}

	}*/

	// Writes a binary raster file - initialisation
	void initBIN(const State &state, Domain const &dom, Parallel const &par, std::string dir){

		std::string filename;
    	realArr data   = realArr("data",dom.nCell);


		//header
		filename = dir+std::to_string(par.myrank)+"result.hdr";
		std::ofstream fOutStream1(filename);
		if (fOutStream1.is_open()){
			//fOutStream << "# vtk DataFile Version 3.0.\nOutput file " << filename <<"\nASCII\nDATASET UNSTRUCTURED_GRID\n";
			fOutStream1 << "ncols " << dom.nx << std::endl;
			fOutStream1 << "nrows " << dom.ny << std::endl;
			fOutStream1 << "xllcorner " << dom.xll + par.i_beg * dom.dxConst << std::endl;
			fOutStream1 << "yllcorner " << dom.yll + dom.ny_glob*dom.dxConst - par.j_beg * dom.dxConst << std::endl;
			fOutStream1 << "cellsize " << dom.dxConst << std::endl;
			fOutStream1 << "nodata_value -9999" << std::endl;
			fOutStream1 << "byteorder msbfirst" << std::endl;
			fOutStream1.close();
		}

		Kokkos::fence();



		//z
    	Kokkos::parallel_for("binwrap_init_z", dom.nCell , KOKKOS_LAMBDA(int iGlob) {
				int ii = dom.getIndex(iGlob);
      	data(iGlob) = state.z(ii);
    	});

		Kokkos::fence();




		filename = dir+std::to_string(par.myrank)+"elevation.bin";
		std::ofstream fOutStream2(filename.c_str(), std::ios::binary);
		fOutStream2.write((char*)&data[0], dom.nCell * sizeof(real));
		fOutStream2.close();

		Kokkos::fence();


		//h
    	Kokkos::parallel_for("binwrap_init_h", dom.nCell , KOKKOS_LAMBDA(int iGlob) {
				int ii = dom.getIndex(iGlob);
      	data(iGlob) = state.h(ii);
    	});
		Kokkos::fence();



		filename = dir+std::to_string(par.myrank)+"result"+std::to_string(numOut)+".bin";
		std::ofstream fOutStream3(filename.c_str(), std::ios::binary);
		fOutStream3.write((char*)&data[0], dom.nCell * sizeof(real));
		fOutStream3.close();

		Kokkos::fence();

		if (par.nranks > 1){
			MPI_Barrier(MPI_COMM_WORLD);
		}

	}

	// Writes a binary raster file
	void outputBIN(const State &state, Domain const &dom, Parallel const &par, std::string dir){

		std::string filename;
    	realArr data   = realArr("data",dom.nCell);

    Kokkos::parallel_for("binwrap_h", dom.nCell , KOKKOS_LAMBDA(int iGlob) {
			int ii = dom.getIndex(iGlob);
      data(iGlob) = state.h(ii);
    });

    Kokkos::fence();

		//MPI_Gatherv(&data[0],cur_proc_data_size, SERGHEI_MPI_REAL, total_data_arr, recvcounts, displs, SERGHEI_MPI_REAL, 0, MPI_COMM_WORLD);

		/*if (rank_ == 0){
			MPI_Gatherv(state.h.get_address_at(0, 0), cur_proc_data_size, SERGHEI_MPI_REAL, total_data_arr, recvcounts, displs, SERGHEI_MPI_REAL, 0, MPI_COMM_WORLD);
		}else{
			MPI_Gatherv(state.h.get_address_at(1, 0), cur_proc_data_size, SERGHEI_MPI_REAL, total_data_arr, recvcounts, displs, SERGHEI_MPI_REAL, 0, MPI_COMM_WORLD);
		}*/

		filename = dir+std::to_string(par.myrank)+"result"+std::to_string(numOut)+".bin";
		std::ofstream fOutStream(filename.c_str(), std::ios::binary);
		fOutStream.write((char*)&data[0], dom.nCell * sizeof(real));
		fOutStream.close();

		if (par.nranks > 1){
			MPI_Barrier(MPI_COMM_WORLD);
		}

	}

  // Write header and initial state for time series files
  int writeTimeSeriesIni (const State &state, Domain const &dom, Parallel const &par,
		      SourceSinkData &ss, surfaceIntegrator &sint, boundaryIntegrator &bint,
		      std::vector<ExtBC> &extbc, std::string dir){

    numObs = 0;
    std::string filename = dir + "domainTimeSeries.out";
    domainOutputFile.open (filename);

    if (domainOutputFile.is_open ()){

	// Write the header
	domainOutputFile << "Time ";
	domainOutputFile << "SurfaceVolume ";



#if SERGHEI_DEBUG_BOUNDARY

	std::cerr << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "Hello from rank: " << par.myrank << "/" << par.nranks << "\n";
	for (int i = 0; i < extbc.size(); i++) {
	  std::cerr << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "boundary cells: " << extbc[i].ncellsBC << "\n";
	  std::cerr << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "inflow discharge: " << extbc[i].inflowDischarge << "\n";
	  std::cerr << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "inflow accumulated: " << extbc[i].inflowAccumulated << "\n";
	  std::cerr << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "outflow discharge: " << extbc[i].outflowDischarge << "\n";
	  std::cerr << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "outflow accumulated: " << extbc[i].outflowAccumulated << "\n";
	}

#endif

	bint.integrate(extbc,dom,1);

#if SERGHEI_DEBUG_BOUNDARY
	if (par.masterproc){
	    std::cerr << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "boundary cells (reduced): " << bint.ncellsBC << "\n";
	    std::cerr << GGD  << GRAY << __PRETTY_FUNCTION__ << RESET << "inflow discharge (integrated) " << bint.inflowDischargeG << "\n";
	    std::cerr << GGD  << GRAY << __PRETTY_FUNCTION__ << RESET << "inflow accumulated (integrated) " << bint.inflowAccumulatedG << "\n";
	    std::cerr << GGD  << GRAY << __PRETTY_FUNCTION__ << RESET << "outflow discharge (integrated) " << bint.outflowDischargeG << "\n";
	    std::cerr << GGD  << GRAY << __PRETTY_FUNCTION__ << RESET << "outflow accumulated (integrated) " << bint.outflowAccumulatedG << std::endl;
	  }
#endif



	if(bint.ncellsBC){
		// global inflow/outflow
	  domainOutputFile << "BoundaryInflow ";
	  domainOutputFile << "BoundaryInflowAccum ";
	  domainOutputFile << "BoundaryOutflow ";
	  domainOutputFile << "BoundaryOutflowAccum ";
		// boundary flow per BC
		for (int i = 0; i < extbc.size(); i ++) {
			domainOutputFile << "BoundaryFlow_" << i << " " ;
			domainOutputFile << "BoundaryAccum_" << i << " " ;
		}
	}

	if (dom.isRain)
	  {
	    domainOutputFile << "RainFlux ";
	    domainOutputFile << "RainAccum ";
	  }

	if (ss.inf.model)
	  {
	    domainOutputFile << "InfFlux ";
	    domainOutputFile << "InfAccum ";
	  }

	domainOutputFile << std::endl;

      }

    else
      {
	std::cerr << RERROR "Could not create domainTimeSeries.out file" << std::endl ;
	return 0;
      }

    // write the data
    writeTimeSeries (state, dom, par, sint, bint, *bint.extbc);

    numObs++;

    return 1;

  }
  #if SERGHEI_DRAINAGE_MODEL
   int writeDrainageTimeSeriesIni (Parallel const &par,surfaceIntegrator &sint, DrainageDywave &Drain, std::string dir, ReadFileFuncs &RFF, Link &Tlink, Node &Tnode){
    numObs = 0;
    std::string filename = dir + "DrainageTimeSeries.out";
    DrainageOutputFile.open (filename);
    if (DrainageOutputFile.is_open ()){
    DrainageOutputFile << "Time ";
    // DrainageOutputFile << "DrainageVolume[L^3] ";
    DrainageOutputFile << "two2oneAccum[L^3] ";
    DrainageOutputFile << "one2twoAccum[L^3] ";
    // DrainageOutputFile << "Vrouted[L^3] ";
    DrainageOutputFile << "OutfallInflow[L^3] ";
    DrainageOutputFile << "OutfallDischarge[L^3/s] ";
    for (int ii = 0; ii<Nobjects[NODE];ii++)
            {
               
				if(Tnode.typee(ii) == JUNCTION)
				{
				  DrainageOutputFile << "JUNCTION_" << ii << "2D-1D"<<"[L^3/s]" << "" ;
				  DrainageOutputFile << "JUNCTION_" << ii << "1D-2D"<<"[L^3/s]" << ""  ;
				}
				else if (Tnode.typee(ii) == INLET)
				{
				  DrainageOutputFile << "INLET_" << ii <<"2D-1D"<< "[L^3/s]" << "" ;
				  DrainageOutputFile << "INLET_" << ii <<"1D-2D"<< "[L^3/s]" << ""  ;
				} 
				else if (Tnode.typee(ii) == OUTFALL)
                {
                  DrainageOutputFile << "OUTFALL_" << ii << "[L^3/s]" << "" ;
                } 
    }
	DrainageOutputFile << std::endl;}

    else
      {
	  std::cerr << RERROR "Could not create DrainageTimeSeries.out file" << std::endl ;
	  return 0;
      }

    // write the data
     if(par.masterproc)  {writeDrainageOutputTimeSeries (par, sint, Drain, RFF, Tlink, Tnode);}

    numObs++;

    return 1;

  }



  void writeDrainageOutputTimeSeries (Parallel const &par, surfaceIntegrator &sint, DrainageDywave &Drain, ReadFileFuncs &RFF, Link &Tlink, Node &Tnode){

    // Write the data
    std::cout.precision(OUTPUT_PRECISION);
    DrainageOutputFile << std::scientific << RFF.dt << " ";
    // DrainageOutputFile << std::scientific << sint.drainageVolumeG << " ";
    DrainageOutputFile << std::scientific << Drain.drainage_two2one << " ";
    DrainageOutputFile << std::scientific << Drain.drainage_one2two << " ";
    // DrainageOutputFile << std::scientific << Drain.Vrouted << " ";
    DrainageOutputFile << std::scientific << Drain.outfallInflow << " ";
    DrainageOutputFile << std::scientific << Drain.outfallDischarge << " ";
    for (int ii = 0; ii<Nobjects[NODE];ii++)
            {
				if(Tnode.typee(ii) == JUNCTION)
				{
				  DrainageOutputFile << std::scientific << Tnode.SDinflow(ii) << " ";
				  DrainageOutputFile << std::scientific << Tnode.SDoutflow(ii) << " ";
				}
				else if (Tnode.typee(ii) == INLET)
				{
				  DrainageOutputFile << std::scientific << Tnode.SDinflow(ii) << " ";
				  DrainageOutputFile << std::scientific << Tnode.SDoutflow(ii) << " ";
				}
                else if (Tnode.typee(ii) == OUTFALL)
                {
                  DrainageOutputFile << std::scientific << Tnode.inflow(ii) << " ";
                } 
    }
    DrainageOutputFile << std::endl;

  }

#endif

  // Writes time series files
  void writeTimeSeries (const State &state, Domain const &dom, Parallel const &par, surfaceIntegrator &sint, boundaryIntegrator &bint, std::vector<ExtBC> &extbc){
    timer.reset();
		for (int i = 0; i < extbc.size(); i ++) extbc[i].reduce(par);
    if(par.masterproc) writeDomainTimeSeries(state,dom,par,sint,bint,extbc);
    numObs++;
    dom.timers.out += timer.seconds();
  }

  void writeDomainTimeSeries (State const &state, Domain const &dom,
			 Parallel const &par, surfaceIntegrator &sint,
			 boundaryIntegrator &bint, const std::vector<ExtBC> &extbc){

    // Write the data
    std::cout.precision(OUTPUT_PRECISION);
    domainOutputFile << std::scientific << dom.etime << " ";
    domainOutputFile << std::scientific << sint.surfaceVolumeG << " ";


    if(bint.ncellsBC){
	  	domainOutputFile << std::scientific << bint.inflowDischargeG << " ";
	   	domainOutputFile << std::scientific << bint.inflowAccumulatedG << " ";
      domainOutputFile << std::scientific << bint.outflowDischargeG << " ";
	   	domainOutputFile << std::scientific << bint.outflowAccumulatedG << " ";
			for (int i = 0; i < extbc.size(); i ++) domainOutputFile << extbc[i].netQ << " ";
			for (int i = 0; i < extbc.size(); i ++) domainOutputFile << extbc[i].netVol << " ";
    }

    if (dom.isRain)
      {
	domainOutputFile << std::scientific << sint.rainFluxG << " ";
	domainOutputFile << std::scientific << sint.rainAccumG << " ";
      }

    if (sint.ss->inf.model)
      {
	domainOutputFile << std::scientific << sint.infFluxG << " ";
	domainOutputFile << std::scientific << sint.infAccumG << " ";
      }

    domainOutputFile << std::endl;

  }

  void closeOutputStreams(){
  	domainOutputFile.close();
	}



  void writeTimerRank(const Parallel &par, const real &time, const std::string &name){
          real *logdata;
          logdata = (real*) malloc(par.nranks * sizeof(real));

          MPI_Gather(&time,1,MPI_DOUBLE,logdata,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
          if(par.masterproc) logFile << name <<"\t\t"  ; for(int ii=0; ii<par.nranks; ii++) logFile << " :\t"<< logdata[ii] ; logFile << std::endl;
          real avg=0;
          for(int ii=0; ii<par.nranks; ii++) avg += logdata[ii];
          avg = avg/par.nranks;
          if(par.masterproc) logFile << name << "_lb\t\t"  ; for(int ii=0; ii<par.nranks; ii++) logFile << " :\t"<< logdata[ii]/avg; logFile << std::endl;

  }


void writeLogFile(Domain const &dom, Parallel const &par, std::string dir){
      real sim, ratio;

      if(par.masterproc){
      std::string filename = dir + "log.out";
      logFile.open(filename);
      if (domainOutputFile.is_open()){
          logFile << "DomainArea [m2]: " << dom.areaGlobal << std::endl;
          logFile << "nCell(SW) : " << dom.nCellGlobal << std::endl;
          logFile << "nCellValid(SW) : " << dom.nCellValidGlobal << std::endl;
          // logFile << "nCell(GW) : " << dom.nCellGlobalGW << std::endl;
          // logFile << "nCellValid(GW) : " << dom.nCellValidGlobalGW << std::endl;

			    ratio = dom.timers.total/dom.timers.total;
			    logFile << "runTime : " << dom.timers.total << " : " << ratio << std::endl;
			    sim = dom.timers.total-dom.timers.init;
			    ratio = sim/dom.timers.total;
			    logFile << "simTime : " << sim << " : " << ratio << std::endl;
			    ratio = dom.timers.init/dom.timers.total;
			    logFile << "initTime : " << dom.timers.init << " : " <<ratio << std::endl;
			    ratio = dom.timers.out/dom.timers.total;
			    logFile << "outputTime : " << dom.timers.out << " : " << ratio <<std::endl;
			    ratio = dom.timers.sweflux/dom.timers.total;
			    logFile << "sweFluxTime : " << dom.timers.sweflux << " : " << ratio << std::endl;
			    ratio = dom.timers.raininf/dom.timers.total;
			    logFile << "rainInfTime : " << dom.timers.raininf << " : " << ratio << std::endl;
			    ratio = dom.timers.swe/dom.timers.total;
			    logFile << "sweNotFluxTime : " << dom.timers.swe << " : " << ratio << std::endl;
				  ratio = dom.timers.exchange / dom.timers.total;
				  logFile << "exchangeTime : " << dom.timers.exchange << " : " << ratio << std::endl;
		      dom.timers.halo -= dom.timers.exchange;
      	  ratio = dom.timers.halo / dom.timers.total;
      	  logFile << "haloTime : " << dom.timers.halo << " : " << ratio << std::endl;
				  ratio = dom.timers.integrate / dom.timers.total;
				  logFile << "integrateTime : " << dom.timers.integrate << " : " << ratio << std::endl;
      	  ratio = dom.timers.integrateMPI / dom.timers.total;
      	  logFile << "integrateMPITime : " << dom.timers.integrateMPI << " : " << ratio << std::endl;
				  ratio = dom.timers.dt / dom.timers.total;
				  logFile << "dtComputeTime : " << dom.timers.dt << " : " << ratio << std::endl;
				  logFile << "sweBCtime : " << dom.timers.sweBC << " : " << dom.timers.sweBC/dom.timers.total << std::endl;

          #if SERGHEI_SUBSURFACE_MODEL
          logFile << " -------- GW Times --------" << std::endl;
          ratio = dom.timers.gw/dom.timers.total;
          logFile << "gwTime : " << dom.timers.gw << " : " << ratio << std::endl;
          ratio = dom.timers.gwlinsys/dom.timers.total;
          logFile << "gwLinSysTime : " << dom.timers.gwlinsys << " : " << ratio << std::endl;
          ratio = dom.timers.gwlinsol/dom.timers.total;
          logFile << "gwLinSolTime : " << dom.timers.gwlinsol << " : " << ratio << std::endl;
          ratio = dom.timers.gwUpdateK/dom.timers.total;
          logFile << "gwUpdateKTime : " << dom.timers.gwUpdateK << " : " << ratio << std::endl;
          ratio = dom.timers.gwUpdateQ/dom.timers.total;
          logFile << "gwUpdateQTime : " << dom.timers.gwUpdateQ << " : " << ratio << std::endl;
          ratio = dom.timers.gwUpdateWC/dom.timers.total;
          logFile << "gwUpdateWCTime : " << dom.timers.gwUpdateWC << " : " << ratio << std::endl;
          ratio = dom.timers.gwBC/dom.timers.total;
          logFile << "gwBCTime : " << dom.timers.gwBC << " : " << ratio << std::endl;
          ratio = dom.timers.gwIntegrate/dom.timers.total;
          logFile << "gwIntegrateTime : " << dom.timers.gwIntegrate << " : " << ratio << std::endl;
          ratio = dom.timers.gwMPI/dom.timers.total;
          logFile << "gwExchangeTime : " << dom.timers.gwMPI << " : " << ratio << std::endl;
          #endif
				#if SERGHEI_LPT
					logFile << "InitializationParticles : " << dom.timers.initialization_particles << " : " << dom.timers.initialization_particles/dom.timers.total << std::endl;
					logFile << "UpdateParticles : " << dom.timers.update_particles << " : " << dom.timers.update_particles/dom.timers.total << std::endl;
					logFile << "OutputParticles : " << dom.timers.output_particles << " : " << dom.timers.output_particles/dom.timers.total << std::endl;
					logFile << "TotalParticles : " << dom.timers.initialization_particles+dom.timers.update_particles+dom.timers.output_particles << " : " << (dom.timers.initialization_particles+dom.timers.update_particles+dom.timers.output_particles)/dom.timers.total << std::endl;
				#endif
          }
      }

          logFile << "TIMERS PER RANK" << std::endl;
          // timers per rank
          writeTimerRank(par,dom.timers.total,"runTime");
    sim = dom.timers.total-dom.timers.init;
          writeTimerRank(par,sim,"simTime");
          writeTimerRank(par,dom.timers.init,"initTime");
          writeTimerRank(par,dom.timers.out,"outputTime");
          writeTimerRank(par,dom.timers.sweflux,"sweFluxTime");
          writeTimerRank(par,dom.timers.raininf,"rainInfTime");
          writeTimerRank(par,dom.timers.swe,"sweNotFluxTime");
			    writeTimerRank(par,dom.timers.exchange,"exchangeTime");
			    writeTimerRank(par,dom.timers.halo,"haloTime");
			    writeTimerRank(par,dom.timers.integrate,"integrateTime");
			    writeTimerRank(par,dom.timers.integrateMPI,"integrateMPITime");
			    writeTimerRank(par,dom.timers.dt,"dtComputeTime");
			    writeTimerRank(par,dom.timers.dt,"sweBCtime");
          #if SERGHEI_SUBSURFACE_MODEL
          writeTimerRank(par,dom.timers.gw,"gwTime");
          writeTimerRank(par,dom.timers.gwlinsys,"gwLinSysTime");
          writeTimerRank(par,dom.timers.gwlinsol,"gwLinSolTime");
		  writeTimerRank(par,dom.cg_iter,"gwLinSolIter");
           writeTimerRank(par,dom.timers.gwUpdateK,"gwUpdateKTime");
           writeTimerRank(par,dom.timers.gwUpdateQ,"gwUpdateQTime");
           writeTimerRank(par,dom.timers.gwUpdateWC,"gwUpdateWCTime");
           writeTimerRank(par,dom.timers.gwBC,"gwBCTime");
           writeTimerRank(par,dom.timers.gwIntegrate,"gwIntegrateTime");
           writeTimerRank(par,dom.timers.gwMPI,"gwExchangeTime");
          #endif
		  #if SERGHEI_DRAINAGE_MODEL
          writeTimerRank(par,dom.timers.drainage,"drainageTime");
          writeTimerRank(par,dom.timers.drainageExchange,"drainageExchangeTime");
          writeTimerRank(par,dom.timers.drainageLinkFlows,"drainageLinkFlowsTime");
          writeTimerRank(par,dom.timers.drainageNodeDepths,"drainageNodeDepthsTime");
          #endif
          if(par.masterproc){
              auto nowtime = std::chrono::system_clock::now();
              std::time_t now_time = std::chrono::system_clock::to_time_t(nowtime);

              logFile << std::endl << "DateTime : " << std::ctime(&now_time) << std::endl;

              char hostbuffer[256];
              int hostname = gethostname(hostbuffer, sizeof(hostbuffer));
              logFile << "Machine : " << hostbuffer << std::endl;
              #ifdef __NVCC__
              cudaDeviceProp deviceProp;
              cudaError_t result = cudaGetDeviceProperties(&deviceProp, 0);
              logFile << "GPU : " << deviceProp.name << std::endl;
              #else
              std::string CPUBrandString;
              CPUBrandString.resize(49);
              uint *CPUInfo = reinterpret_cast<uint*>(CPUBrandString.data());
              for (uint i=0; i<3; i++) __cpuid(0x80000002+i, CPUInfo[i*4+0], CPUInfo[i*4+1], CPUInfo[i*4+2], CPUInfo[i*4+3]);
              CPUBrandString.assign(CPUBrandString.data()); // correct null terminator
              logFile << "CPU : " << CPUBrandString << std::endl;
              #endif
              logFile << "nTasks : " << par.nranks << std::endl;

				// write compilation setup
        logFile << "\nSERGHEI_GIT_VERSION : " << SERGHEI_GIT_VERSION << std::endl;
				logFile << "\n-------------------------\nMODEL COMPONENT SETUP" <<std::endl;
				logFile << "SERGHEI_TOOLS : " << SERGHEI_TOOLS << std::endl;
				logFile << "SERGHEI_SUBSURFACE_MODEL : " << SERGHEI_SUBSURFACE_MODEL << std::endl;
				logFile << "SERGHEI_PARTICLE_TRACKING : " << SERGHEI_LPT << std::endl;
				logFile << "SERGHEI_VEGETATION_MODEL : " << SERGHEI_VEGETATION_MODEL << std::endl;
				logFile << "SERGHEI_FRICTION_MODEL : " << SERGHEI_FRICTION_MODEL << std::endl;
				logFile << "SERGHEI_POINTWISE_FRICTION : " << SERGHEI_POINTWISE_FRICTION << std::endl;
				logFile << "SERGHEI_MAXFLOOD : " << SERGHEI_MAXFLOOD << std::endl;
        logFile << "SERGHEI_REAL : " << SERGHEI_REAL << std::endl;
        logFile << "SERGHEI_NC_MODE : " << SERGHEI_NC_MODE << std::endl;
        logFile << "SERGHEI_NC_REAL : " ;
	      if(SERGHEI_NC_REAL == NC_FLOAT) logFile << "NC_FLOAT";
	      if(SERGHEI_NC_REAL == NC_DOUBLE) logFile << "NC_DOUBLE";
 	      logFile << std::endl;
        logFile << "SERGHEI_WRITE_HZ : " << SERGHEI_WRITE_HZ << std::endl;
				logFile << "SERGHEI_WRITE_SUBDOMS : " << SERGHEI_WRITE_SUBDOMS << std::endl;
		logFile << "SERGHEI_INPUT_NETCDF : " << SERGHEI_INPUT_NETCDF << std::endl;
		    logFile << "SERGHEI_NC_ENABLE_NAN : " << SERGHEI_NC_ENABLE_NAN << std::endl;

      logFile << "\n-------------------------\nDEBUG FLAGS" <<std::endl;
      logFile << "SERGHEI_DEBUG_PARALLEL_DECOMPOSITION : " << SERGHEI_DEBUG_PARALLEL_DECOMPOSITION << std::endl;
      logFile << "SERGHEI_DEBUG_WORKFLOW : " << SERGHEI_DEBUG_WORKFLOW << std::endl;
      logFile << "SERGHEI_DEBUG_KOKKOS_SETUP : " << SERGHEI_DEBUG_KOKKOS_SETUP << std::endl;
      logFile << "SERGHEI_DEBUG_BOUNDARY : " << SERGHEI_DEBUG_BOUNDARY << std::endl;
      logFile << "SERGHEI_DEBUG_DT : " << SERGHEI_DEBUG_DT << std::endl;
      logFile << "SERGHEI_DEBUG_TOOLS : " << SERGHEI_DEBUG_TOOLS << std::endl;
      logFile << "SERGHEI_DEBUG_MASS_CONS : " << SERGHEI_DEBUG_MASS_CONS << std::endl;
      logFile << "SERGHEI_DEBUG_INFILTRATION : " << SERGHEI_DEBUG_INFILTRATION << std::endl;
      logFile << "SERGHEI_DEBUG_MPI : " << SERGHEI_DEBUG_MPI << std::endl;
			logFile << "SERGHEI_DEBUG_OUTPUT : " << SERGHEI_DEBUG_OUTPUT << std::endl;
          }

  }


  // Writes VTK file for subsurface
    #if SERGHEI_SUBSURFACE_MODEL

    void outputInitNETCDFSub(const GwState &gw, GwDomain const &gdom, Parallel const &par, std::string dir) {
        int dimids[4], nxhalo, nyhalo;
        MPI_Offset st[3], ct[3];
        realArr xCoord = realArr("xCoord",gdom.nx);
        realArr yCoord = realArr("yCoord",gdom.ny);
        realArr zCoord = realArr("zCoord",gdom.nz);
        realArr data   = realArr("data",gdom.nCell);
        static char title[] = "seconds" ;
        std::string filename;

        filename=dir+"output_subsurface.nc";

        // Create the file
        ncwrap( ncmpi_create( MPI_COMM_WORLD , filename.c_str() , NC_CLOBBER , MPI_INFO_NULL , &ncid ) , __LINE__ );

        // Create the dimensions
        ncwrap( ncmpi_def_dim( ncid , "t" , (MPI_Offset) NC_UNLIMITED , &tDim ) , __LINE__ );
        ncwrap( ncmpi_def_dim( ncid , "x" , (MPI_Offset) gdom.nx_glob  , &xDim ) , __LINE__ );
        ncwrap( ncmpi_def_dim( ncid , "y" , (MPI_Offset) gdom.ny_glob  , &yDim ) , __LINE__ );
        ncwrap( ncmpi_def_dim( ncid , "z" , (MPI_Offset) gdom.nz_glob  , &zDim ) , __LINE__ );
        // Create the variables
        dimids[0] = tDim;
        ncwrap( ncmpi_def_var( ncid , "t"      , NC_DOUBLE , 1 , dimids , &tVar ) , __LINE__ );
        ncwrap( ncmpi_put_att_text (ncid, tVar, "units",strlen(title), title), __LINE__ );
        dimids[0] = xDim;
        ncwrap( ncmpi_def_var( ncid , "x"      , NC_DOUBLE , 1 , dimids , &xVar ) , __LINE__ );
        dimids[0] = yDim;
        ncwrap( ncmpi_def_var( ncid , "y"      , NC_DOUBLE , 1 , dimids , &yVar ) , __LINE__ );
        dimids[0] = zDim;
        ncwrap( ncmpi_def_var( ncid , "z"      , NC_DOUBLE , 1 , dimids , &zVar ) , __LINE__ );

        dimids[0] = zDim; dimids[1] = yDim; dimids[2] = xDim;
        ncwrap( ncmpi_def_var( ncid , "z3d" , NC_DOUBLE , 3 , dimids , &z3Var  ) , __LINE__ );
        dimids[0] = tDim; dimids[1] = zDim; dimids[2] = yDim; dimids[3] = xDim;
        ncwrap( ncmpi_def_var( ncid , "hd" , NC_DOUBLE , 4 , dimids , &hdVar  ) , __LINE__ );
        ncwrap( ncmpi_def_var( ncid , "wc" , NC_DOUBLE , 4 , dimids , &wcVar  ) , __LINE__ );

        // End "define" mode
        ncwrap( ncmpi_enddef( ncid ) , __LINE__ );

        // Compute x, y, z coordinates
        Kokkos::parallel_for( gdom.nx , KOKKOS_LAMBDA(int i) {
            xCoord(i) = gdom.xll + ( par.i_beg + i + 0.5) * gdom.dx;
        });
        Kokkos::parallel_for( gdom.ny , KOKKOS_LAMBDA(int j) {
            yCoord(j) = gdom.yll + gdom.ny_glob*gdom.dx - ( par.j_beg + j + 0.5) * gdom.dx;
        });
        Kokkos::parallel_for( gdom.nz , KOKKOS_LAMBDA(int k) {
            zCoord(k) = -(k+0.5)*gdom.dz(k);
        });
        Kokkos::fence();

        // Write out x, y coordinates
        st[0] = par.i_beg;
        ct[0] = gdom.nx;
        ncwrap( ncmpi_put_vara_double_all( ncid , xVar , st , ct , xCoord.data() ) , __LINE__ );
        st[0] = par.j_beg;
        ct[0] = gdom.ny;
        ncwrap( ncmpi_put_vara_double_all( ncid , yVar , st , ct , yCoord.data() ) , __LINE__ );
        st[0] = 0;
        ct[0] = gdom.nz;
        ncwrap( ncmpi_put_vara_double_all( ncid , zVar , st , ct , zCoord.data() ) , __LINE__ );

        // Write z for the 3D domain
        nxhalo = gdom.nx + 2*hc;
        nyhalo = gdom.ny + 2*hc;
        st[0] = 0;          st[1] = par.j_beg;  st[2] = par.i_beg;
        ct[0] = gdom.nz;    ct[1] = gdom.ny;    ct[2] = gdom.nx;
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*nxhalo*nyhalo + (hc+jj)*nxhalo + ii + hc;
            data(idom) = gdom.z(iGlob);
        });
        Kokkos::fence();
        ncwrap( ncmpi_put_vara_double_all( ncid , z3Var  , st , ct , data.data() ) , __LINE__ );

        writeGwNETCDF(gw, gdom, par);
        ncwrap( ncmpi_close(ncid) , __LINE__ );
    }

    void outputNETCDFSubsurface(const GwState &gw, GwDomain const &gdom, Parallel const &par, std::string dir) {
        std::string filename;
        filename=dir+"output_subsurface.nc";
        // Create the file
        ncwrap( ncmpi_open( MPI_COMM_WORLD , filename.c_str() , NC_WRITE , MPI_INFO_NULL , &ncid ) , __LINE__ );
        ncwrap( ncmpi_inq_varid( ncid , "hd" , &hdVar  ) , __LINE__ );
        ncwrap( ncmpi_inq_varid( ncid , "wc" , &wcVar  ) , __LINE__ );

        writeGwNETCDF(gw, gdom, par);

        ncwrap( ncmpi_close(ncid) , __LINE__ );
    }

    void writeGwNETCDF(const GwState &gw, GwDomain const &gdom, Parallel const &par) {
        realArr data = realArr("data",gdom.nCell);
        MPI_Offset st[4], ct[4];
        double timeIter[numOut+1];
        int nxhalo = gdom.nx + 2*hc, nyhalo = gdom.ny + 2*hc;
        //write t. As the first one is written in the first iteration we should add +1
        for (int i=0; i<numOut+1; i++) { timeIter[i] = i*1.0;}
        st[0] = 0;
     	ct[0] = numOut+1;
        ncwrap( ncmpi_put_vara_double_all( ncid , tVar ,  st , ct , timeIter ) , __LINE__ );

        // st[0] = numOut; st[1] = par.i_beg; st[2] = par.j_beg;   st[3] = 0;
        // ct[0] = 1     ; ct[1] = gdom.nx  ; ct[2] = gdom.ny  ;   ct[3] = gdom.nz;

        st[0] = numOut; st[1] = 0;       st[2] = par.j_beg;  st[3] = par.i_beg;
        ct[0] = 1     ; ct[1] = gdom.nz; ct[2] = gdom.ny  ;  ct[3] = gdom.nx  ;

        // Kokkos::parallel_for(gdom.nCell, kernel_output<dspace>(gw.h, data, gdom));
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*nxhalo*nyhalo + (hc+jj)*nxhalo + ii + hc;
            data(idom) = gw.h(iGlob,1);
        });
        Kokkos::fence();
        ncwrap( ncmpi_put_vara_double_all( ncid , hdVar , st , ct , data.data() ) , __LINE__ );

        // Kokkos::parallel_for(gdom.nCell, kernel_output<dspace>(gw.wc, data, gdom));
        Kokkos::parallel_for( gdom.nCell , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            gdom.unpackIndicesGw(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*nxhalo*nyhalo + (hc+jj)*nxhalo + ii + hc;
            data(idom) = gw.wc(iGlob,1);
        });
        Kokkos::fence();
        ncwrap( ncmpi_put_vara_double_all( ncid , wcVar , st , ct , data.data() ) , __LINE__ );


    }


  // VTK writer for the subsurface
  void outputVTKsubsurface(const GwState &gw, GwDomain const &gdom, Parallel const &par, std::string dir){

      std::string filename;
      filename = dir+"result_subsurface"+std::to_string(numOut)+".vtk";
      real *xCoord;
      real *yCoord;
      real *zCoord;
      real *data_cpu;

      int ncells=gdom.ny*gdom.nx*gdom.nz;

      int nVars=2; //2 variables to write: h, wc

      realArr data  = realArr("data",nVars*ncells);
      #ifdef __NVCC__
          cudaMallocHost( &data_cpu , nVars*ncells*sizeof(real) );
      #else
          data_cpu = data.data();
      #endif

      xCoord=(real*) malloc((gdom.nx+1)*sizeof(real));
      yCoord=(real*) malloc((gdom.ny+1)*sizeof(real));
      zCoord=(real*) malloc((gdom.nz+1)*sizeof(real));

      // Compute x, y coordinates
      for (int i=0; i<gdom.nx+1; i++) {
      xCoord[i] = gdom.xll + ( par.i_beg + i) * gdom.dx;
      };

      for (int j=0; j<gdom.ny+1; j++) {
      yCoord[j] = gdom.yll + gdom.ny_glob*gdom.dx - ( par.j_beg + j) * gdom.dy;
      };

      int nnodes=(gdom.nx+1)*(gdom.ny+1)*(gdom.nz+1);

      // dualDbl::t_host h_h = gstate.h.h_view;
      // dualDbl::t_host h_wc = gstate.wc.h_view;
      Kokkos::parallel_for( ncells , KOKKOS_LAMBDA(int idom) {
          int i,j,k;
          int nxhalo = gdom.nx+2*hc;
          int nyhalo = gdom.ny+2*hc;
          gdom.unpackIndicesGw(idom,gdom.nz,gdom.ny,gdom.nx,k,j,i);
          int iGlob = (hc+k)*nxhalo*nyhalo + (hc+j)*nxhalo + i + hc;
          data(idom) = gw.h(iGlob,1);
          data(idom+ncells) = gw.wc(iGlob,1);
          // get z coordinate
          if (i == 0 & j == 0) {zCoord[k] = gdom.z(iGlob);}
      });
      Kokkos::fence();

      #ifdef __NVCC__
          cudaMemcpyAsync( data_cpu , data.data() , nVars*ncells*sizeof(real) , cudaMemcpyDeviceToHost );
          cudaDeviceSynchronize();
      #endif

      int i,j,k,iGlob;

      std::ofstream fOutStream(filename);
      if (fOutStream.is_open()){
          //fOutStream << "# vtk DataFile Version 3.0.\nOutput file " << filename <<"\nASCII\nDATASET UNSTRUCTURED_GRID\n";
          fOutStream << "# vtk DataFile Version 3.0.\nOutputfile\nASCII\nDATASET UNSTRUCTURED_GRID\n";

          fOutStream << "POINTS " << nnodes  << " double\n";
          for (k=0; k<=gdom.nz; k++){
              for(j=0;j<=gdom.ny;j++){
                  for(i=0;i<=gdom.nx;i++){
                      fOutStream << std::setprecision(9) << xCoord[i] << " " << yCoord[j] <<" "<< zCoord[k]<< " 0.0\n";
                  }
              }
          }
          fOutStream << "CELLS " << ncells << " " << 2*ncells <<"\n";
          for(iGlob=0;iGlob<ncells;iGlob++){
             int ii,jj,kk;
              gdom.unpackIndicesGw(iGlob,gdom.nz,gdom.ny,gdom.nx,kk,jj,ii);
              fOutStream << "4 " << jj*(gdom.nx+1)+ii << " " << (jj+1)*(gdom.nx+1)+ii << " " << (jj+1)*(gdom.nx+1)+ii+1 << " " << jj*(gdom.nx+1)+ii+1 << "\n";
          }

          fOutStream << "CELL_TYPES " << ncells << "\n";
          for(i=0;i<ncells;i++){
              fOutStream << "9\n";
          }

          fOutStream << "CELL_DATA " << ncells << "\n";
          fOutStream << "SCALARS h double\n";
          fOutStream << "LOOKUP_TABLE default\n";
          for(iGlob=0;iGlob<ncells;iGlob++){
              fOutStream << std::setprecision(9) << data_cpu[iGlob] << "\n";
          }

          fOutStream << "SCALARS wc double\n";
          fOutStream << "LOOKUP_TABLE default\n";
          for(iGlob=0;iGlob<ncells;iGlob++){
              fOutStream << std::setprecision(9) << data_cpu[iGlob+ncells]  << "\n";
          }

    fOutStream.close();
      }


  }



  // Write header and initial state for time series files
  int writeSubTimeSeriesIni (GwDomain const &gdom, GwIntegrator const &gint, Parallel const &par, std::string dir){

        numObs = 0;
        std::string filename = dir + "SubsurfaceTimeSeries.out";
        SubsurfaceOutputFile.open (filename);

        if (SubsurfaceOutputFile.is_open ()){
            // Write the header
            SubsurfaceOutputFile << "Time ";
            SubsurfaceOutputFile << "SubSurfaceVolume [m3] ";
            SubsurfaceOutputFile << "ExchangeRate [m3/s] ";
            SubsurfaceOutputFile << "BoundaryInflow ";
            SubsurfaceOutputFile << "BoundaryOutflow ";
            SubsurfaceOutputFile << "Source/SinkInflow [m3/s] ";
            SubsurfaceOutputFile << "Source/SinkOutflow [m3/s] ";
            SubsurfaceOutputFile << std::endl;
        }
        else
        {std::cerr << RERROR "Could not create domainTimeSeries.out file" << std::endl; return 0;}

        if(par.masterproc)  {writeSubsurfaceTimeSeries(gdom, gint);}
        numObs++;
        return 1;
  }



  void writeSubsurfaceTimeSeries (GwDomain const &gdom, GwIntegrator const &gint){
    // Write the data
    std::cout.precision(OUTPUT_PRECISION);
    SubsurfaceOutputFile << std::scientific << gdom.etime << " ";
    SubsurfaceOutputFile << std::scientific << gint.Vtot_glob << " ";
    SubsurfaceOutputFile << std::scientific << gint.Vexch_glob << " ";
    SubsurfaceOutputFile << std::scientific << gint.QinBC_glob << " ";
    SubsurfaceOutputFile << std::scientific << gint.QoutBC_glob << " ";
    SubsurfaceOutputFile << std::scientific << gint.QinSS_glob << " ";
    SubsurfaceOutputFile << std::scientific << gint.QoutSS_glob << " ";
    SubsurfaceOutputFile << std::endl;

  }

  #endif

// parallel netcdf input functionality

// generic error handling function
static void handle_error(Parallel &par, int status, int lineno) {

  if (par.masterproc) std::cerr << RERROR << "Error at line: " << lineno << ": " << ncmpi_strerror(status) << std::endl;
  MPI_Abort(MPI_COMM_WORLD, 1);

}

#if SERGHEI_INPUT_NETCDF
int readNetCDFheader(const Parallel &par, ncStream &nc, Domain &dom){
	MPI_Offset dimsize;

  // open the netcdf file
  ncwrap(ncmpi_open(MPI_COMM_WORLD, nc.fname.c_str(), NC_NOWRITE, MPI_INFO_NULL, &nc.id),__LINE__,par.myrank);

  if (par.masterproc) std::cout << GOK << "Read header from " << nc.fname << "(rank " << par.myrank << ")" << std::endl;

  ncwrap(ncmpi_inq(nc.id, &nc.ndims, &nc.nvars, &nc.ngatts, &nc.unlimited),__LINE__,par.myrank);

  #if SERGHEI_DEBUG_INPUT_NETCDF
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "\tNetCDF ndims: " << nc.ndims << std::endl;
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "\tNetCDF nvars: " << nc.nvars << std::endl;
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "\tNetCDF ngatts: " << nc.ngatts << std::endl;
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "\tNetCDF unlimited: " << nc.unlimited << std::endl;
	#endif

	// get the size of the domain
	ncwrap(ncmpi_inq_dimid(nc.id, "x", &nc.dimids[2]),__LINE__,par.myrank);
	ncwrap(ncmpi_inq_dimid(nc.id, "y", &nc.dimids[1]),__LINE__,par.myrank);
	ncwrap(ncmpi_inq_dimid(nc.id, "t", &nc.dimids[0]),__LINE__,par.myrank);

	// read time dimension
	ncwrap(ncmpi_inq_dimlen(nc.id, nc.dimids[0], &dimsize),__LINE__,par.myrank);
	if(dimsize != 1){
		std::cerr << RERROR << "NetCDF input file " << nc.fname << " has a time dimension of " << dimsize << " and it should be 1." << std::endl;
		std::cerr << RED << "Try using ncks to extract the time slice that you are interested, e.g.," << std::endl << "\tncks -d t,10,10,1 input.nc input.nc" << RESET << std::endl;
		return 0;
	}
	int varid;
	ncwrap(ncmpi_inq_varid(nc.id,"t",&varid),__LINE__,par.myrank);
	ncwrap(ncmpi_get_var_double_all(nc.id,varid,&dom.startTime),__LINE__,par.myrank);
	if(par.masterproc) std::cout << BDASH << "Current time read from NetCDF input: " << dom.startTime << std::endl; 

	// read x dimension
	ncwrap(ncmpi_inq_dimlen(nc.id, nc.dimids[2], &dimsize),__LINE__,par.myrank);
	dom.nx_glob = (int) dimsize;

	// read y dimension
	ncwrap(ncmpi_inq_dimlen(nc.id, nc.dimids[1], &dimsize),__LINE__,par.myrank);
	dom.ny_glob = (int) dimsize;
	if(par.masterproc) std::cout << GOK << "Domain dimensions: " << dom.nx_glob << " " << dom.ny_glob << std::endl;

	nc.ndata = dom.nx_glob * dom.ny_glob;

	return 1;
}

int readNetCDFcoordinates(const Parallel &par, ncStream &nc, Domain &dom){
	int varid;
  char varname[NC_MAX_NAME+1];
  nc_type vtype;
  int var_ndims, var_natts;
  int dimids[3];


	doubleArr data = realArr("databuffer",dom.nx_glob);

	ncwrap(ncmpi_inq_varid(nc.id,"x",&varid),__LINE__,par.myrank);
  ncwrap(ncmpi_inq_var(nc.id, varid, varname, &vtype, &var_ndims, dimids, &var_natts),__LINE__,par.myrank);

  //ncwrap(ncmpi_get_var_double_all(ncin,varid,data),__LINE__,par.myrank);

  ncwrap(ncmpi_get_var_double_all(nc.id,varid,data.data()),__LINE__,par.myrank);

	// find the westmost corner, xll
	Kokkos::parallel_reduce( dom.nx_glob , KOKKOS_LAMBDA (int ii, double &x) {
       x = min(x,data(ii));
    } , Kokkos::Min<double>(dom.xll) );
	Kokkos::fence();

	dom.dxConst = data(1)-data(0);

	// the y-coordinates
	Kokkos::resize(data,dom.ny_glob);
	ncwrap(ncmpi_inq_varid(nc.id,"y",&varid),__LINE__,par.myrank);
  ncwrap(ncmpi_get_var_double_all(nc.id,varid,data.data()),__LINE__,par.myrank);
	// find the southmost corner, yll
	Kokkos::parallel_reduce( dom.ny_glob , KOKKOS_LAMBDA (int ii, double &y) {
       y =  min(y,data(ii));
    } , Kokkos::Min<double>(dom.yll) );

	Kokkos::fence();
	// find the cell vertex, instead of cell center
	dom.xll -= 0.5 * dom.dxConst;
	dom.yll -= 0.5 * dom.dxConst;
	
	if(par.masterproc){
    std::cout << GOK << "dx = " << dom.dx() << std::endl;
	  std::cout << GOK << "Domain extent : (" << dom.xll << ", " << dom.yll << ") (" << dom.xll + dom.dx() * dom.nx_glob << ", " << dom.yll + dom.dx() * dom.ny_glob<< ")" << std::endl;
    std::cout << GOK << "Coordinates read from " << nc.fname << std::endl;
  }
	return 1;
}

int readNetCDFvariable(const Parallel &par, Domain &dom, State &state, ncStream &nc, std::string vname) {
	#if SERGHEI_DEBUG_INPUT_NETCDF
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "\tReading NetCDF variable " << vname << std::endl;
    #endif

  int varid=-1;

  int ncerr=ncmpi_inq_varid(nc.id,vname.c_str(),&varid);
	if(ncerr != NC_NOERR){
    if(ncerr == NC_ENOTVAR){  // if variable not found, return and deal with it one level up
      return NC_ENOTVAR;
      }
    else{
      ncwrap(ncerr,__LINE__,par.myrank);
    }
  }
  
  std::string attName = "_FillValue";
  float nodata;
  ncwrap(ncmpi_get_att_float(nc.id,varid,attName.c_str(),&nodata),__LINE__,par.myrank);
  if(par.masterproc) std::cout << BDASH << "No data value for variable " << vname << " is " << nodata << std::endl;
	
	
	#if SERGHEI_DEBUG_INPUT_NETCDF
		std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "\tTarget NetCDF variable " << vname << " has index " << varid << std::endl;
    #endif

  int ndim=-1;
  ncwrap(ncmpi_inq_varndims(nc.id, varid, &ndim),__LINE__,par.myrank);

  MPI_Offset st[ndim], ct[ndim];
  if(ndim > 2){
    // This is the approach to read time-dependent varibles
    st[0] = 0; st[1] = par.j_beg; st[2] = par.i_beg;
    ct[0] = 1; ct[1] = dom.ny   ; ct[2] = dom.nx   ;
  }
  else{
    // Read a time-independent variable
    st[0] = par.j_beg; st[1] = par.i_beg;
    ct[0] = dom.ny   ; ct[1] = dom.nx   ;
  }
  
  // we use a buffer view, because of the order of coordinates
  realArr data = realArr("ncdata",dom.nCell);
  ncwrap(ncmpi_get_vara_double_all(nc.id,varid,st,ct,data.data()),__LINE__,par.myrank);
  
  int var = getIOvarID(vname);

  #if SERGHEI_DEBUG_INPUT_NETCDF
  std::cout << "rank = " << par.myrank << "\tst : " << st[0] << "\t" << st[1] << "\t" << st[2]  << std::endl;
  std::cout << "rank = " << par.myrank << "\tct : " << ct[0] << "\t" << ct[1] << "\t" << ct[2]  <<std::endl;
  
  std::cout << GGD << vname << " is internal var " << var << std::endl;
  #endif
  if(var < 0 ){
		std::cerr << RERROR << "Internal variable " << vname << " not found." << std::endl;
		return 0;
	} 
 
  if(var == ioZ){
    Kokkos::parallel_reduce("readNC_z_reduce",dom.nCell, KOKKOS_LAMBDA(int iGlob, int &nValid) {
		int ii = dom.getIndex(iGlob);
	  	state.z(ii) = data(iGlob);
	  	if (state.z(ii) <= nodata) {
	   	 state.isnodata(ii) = true;
	  	} else {
	   	 state.isnodata(ii) = false;
		 nValid++;
	  	}
	}, Kokkos::Sum<int>(dom.nCellValid));
  } 
  else{
  Kokkos::parallel_for("netCDFvariable",dom.nCell, KOKKOS_LAMBDA(int iGlob) {

	int ii = dom.getIndex(iGlob);
    if(!std::isnan(data[iGlob]) && !state.isnodata(ii)){// data[iGlob]=0;
        if(var == ioH) state.h(ii) = data[iGlob];
        if(var == ioU) state.hu(ii) = data[iGlob]*state.h(ii);
        if(var == ioV) state.hv(ii) = data[iGlob]*state.h(ii);
				if(var == ioR) state.roughness(ii) = data[iGlob]; 
    }
	});
}


  if(par.masterproc) std::cout << GOK << "NetCDF variable " << GREEN << BOLD << vname << RESET << " read" << std::endl;
  return 1;

}
#endif

};
#endif
