#ifndef _FILEIO_H_
#define _FILEIO_H_

#include <chrono>
#include <ctime>
#include <unistd.h>
#include <cpuid.h>
#include "define.h"
#include "State.h"
#include "SWSourceSink.h"
#include "BC.h"
#include "pnetcdf.h"
#include "mpi.h"
#include "Indexing.h"
#include "DomainIntegrator.h"
#include "tools.h"
#include "ParticleTracking.h"

#include "GwState.h"
#include "GwDomain.h"

#ifndef SERGHEI_NC_MODE
#define SERGHEI_NC_MODE NC_CLOBBER
#endif

#ifndef SERGHEI_WRITE_HZ
#define SERGHEI_WRITE_HZ 0
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
#endif
#if SERGHEI_NC_REAL == NC_FLOAT
	#define ncmpi_put_vara_real_all ncmpi_put_vara_float_all
#endif

#define SERGHEI_NC_ENABLE_MISSING_VALUE 0
#define SERGHEI_NC_MISSING_VALUE "missing_value"

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

/*
#define EMPTY()
#define DEFER(x) x EMPTY()
#define PASTER(x,y) x ## y
#define EVALUATOR(x,y) PASTER(x,y)
#define CHOOSE_BACKEND(name,id) EVALUATOR(name,id)
*/



class FileIO {

protected:

  int ncid;
  int tDim, xDim, yDim, zDim;
  int tVar, xVar, yVar, hVar, hzVar, uVar, vVar, zVar, z3Var, hdVar, wcVar;
  int infVar,infVolVar;
  std::ofstream domainOutputFile;
  std::ofstream SubsurfaceOutputFile;
  std::ofstream logFile;
	int nOut; // expected number of spatial outputs
	#if SERGHEI_MAXFLOOD > 0
	int hMaxVar, momMaxVar, timehMaxVar;
	#endif

private:
  Kokkos::Timer timer;

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

  real outFreq;
  int outFormat;
  int numOut; // spatial output counter

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
		numOut=0;
		if(outFormat==OUT_VTK){
			outputVTKsubsurface(gw, gdom, par, dir);
		}
        else if (outFormat == OUT_NETCDF)   {
            outputInitNETCDFSub(gw, gdom, par, dir);
        }
		numOut++;
	}

    void outputSubsurface(const GwState &gw, GwDomain const &gdom, Parallel const &par, std::string dir){
        numOut--;
        if(outFormat==OUT_VTK){
			outputVTKsubsurface(gw, gdom, par, dir);
		}
        else if (outFormat == OUT_NETCDF)   {
            outputNETCDFSubsurface(gw, gdom, par, dir);
        }
        numOut++;
	}
    #endif

	// Writes spatial fields
	void output(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par,std::string dir){
    timer.reset();

		if(outFormat==OUT_NETCDF){
			outputNETCDF(state,dom,ss,par,dir);
		}
		if(outFormat==OUT_VTK){
			outputVTK(state,dom,ss,par,dir);
		}
		if(outFormat==OUT_BIN){
			outputBIN(state,dom,par,dir);
		}

		numOut++;
    dom.timers.out += timer.seconds();
	}


	// NetCDF initialiser and writer
  void outputInitNETCDF(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par, std::string dir) {
    int dimids[3];
    MPI_Offset st[3], ct[3];
    realArr xCoord = realArr("xCoord",dom.nx);
    realArr yCoord = realArr("yCoord",dom.ny);
    ncArr data   = ncArr("data",dom.ny*dom.nx);
		static char timeUnits[] = "seconds" ;
	std::string filename;

	filename=dir+"output.nc";
    // Create the file
    ncwrap( ncmpi_create( MPI_COMM_WORLD , filename.c_str() , SERGHEI_NC_MODE , MPI_INFO_NULL , &ncid ) , __LINE__ );

    // Create the dimensions
    ncwrap( ncmpi_def_dim( ncid , "t" , (MPI_Offset) NC_UNLIMITED , &tDim ) , __LINE__ );
    ncwrap( ncmpi_def_dim( ncid , "x" , (MPI_Offset) dom.nx_glob  , &xDim ) , __LINE__ );
    ncwrap( ncmpi_def_dim( ncid , "y" , (MPI_Offset) dom.ny_glob  , &yDim ) , __LINE__ );
    // Create the variables
    dimids[0] = tDim;
    ncwrap( ncmpi_def_var( ncid , "t"      , NC_DOUBLE , 1 , dimids , &tVar ) , __LINE__ );
	 	ncwrap( ncmpi_put_att_text (ncid, tVar, "units",strlen(timeUnits), timeUnits), __LINE__ );
    dimids[0] = xDim;
    ncwrap( ncmpi_def_var( ncid , "x"      , NC_DOUBLE , 1 , dimids , &xVar ) , __LINE__ );
    dimids[0] = yDim;
    ncwrap( ncmpi_def_var( ncid , "y"      , NC_DOUBLE , 1 , dimids , &yVar ) , __LINE__ );
		// time dependend variables
    dimids[0] = tDim; dimids[1] = yDim; dimids[2] = xDim;
    ncwrap( ncmpi_def_var( ncid , "h" , SERGHEI_NC_REAL , 3 , dimids , &hVar  ) , __LINE__ );
#if SERGHEI_WRITE_HZ
    ncwrap( ncmpi_def_var( ncid , "h+z" , SERGHEI_NC_REAL , 3 , dimids , &hzVar  ) , __LINE__ );
#endif
    ncwrap( ncmpi_def_var( ncid , "u"      , SERGHEI_NC_REAL , 3 , dimids , &uVar  ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "v"      , SERGHEI_NC_REAL , 3 , dimids , &vVar  ) , __LINE__ );
    if(ss.inf.model){
      ncwrap( ncmpi_def_var( ncid , "inf" , SERGHEI_NC_REAL , 3 , dimids , &infVar  ) , __LINE__ );
      ncwrap( ncmpi_def_var( ncid , "infVol" , SERGHEI_NC_REAL , 3 , dimids , &infVolVar  ) , __LINE__ );
    }

		#if SERGHEI_MAXFLOOD > 0
		int nc_ndims = 2+SERGHEI_MAXFLOOD-1;
		if(nc_ndims==2){
			dimids[0] = yDim;
			dimids[1] = xDim;
		}
		ncwrap( ncmpi_def_var( ncid , "hMax"    , NC_DOUBLE , nc_ndims , dimids , &hMaxVar ) , __LINE__ );
		ncwrap( ncmpi_def_var( ncid , "momMax"    , NC_DOUBLE , nc_ndims , dimids , &momMaxVar ) , __LINE__ );
		ncwrap( ncmpi_def_var( ncid , "timehMax"    , NC_DOUBLE , nc_ndims , dimids , &timehMaxVar ) , __LINE__ );
		#endif

		dimids[0] = yDim; dimids[1] = xDim;
    ncwrap( ncmpi_def_var( ncid , "z"    , SERGHEI_NC_REAL , 2 , dimids , &zVar ) , __LINE__ );

		// define global attributes
		static char title[] = "SERGHEI simulation";
		ncwrap( ncmpi_put_att(ncid, NC_GLOBAL, "title", NC_CHAR,strlen(title)+1,title), __LINE__);

    // End "define" mode
    ncwrap( ncmpi_enddef( ncid ) , __LINE__ );

    // Compute x, y coordinates
    // for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for("compute_grid_coord_x", dom.nx , KOKKOS_LAMBDA(int i) {
      xCoord(i) = dom.xll + ( par.i_beg + i + 0.5) * dom.dx;
    });
    // for (int j=0; j<dom.ny; j++) {
    Kokkos::parallel_for("compute_grid_coord_x", dom.ny , KOKKOS_LAMBDA(int j) {
      yCoord(j) = dom.yll + dom.ny_glob*dom.dx - ( par.j_beg + j + 0.5) * dom.dx;
    });
    Kokkos::fence();

    // Write out x, y coordinates
    st[0] = par.i_beg;
    ct[0] = dom.nx;
    ncwrap( ncmpi_put_vara_double_all( ncid , xVar , st , ct , xCoord.data() ) , __LINE__ );
    st[0] = par.j_beg;
    ct[0] = dom.ny;
    ncwrap( ncmpi_put_vara_double_all( ncid , yVar , st , ct , yCoord.data() ) , __LINE__ );

    st[0] = par.j_beg; st[1] = par.i_beg;
    ct[0] = dom.ny   ; ct[1] = dom.nx   ;

		// write elevation
		Kokkos::parallel_for("ncwrap_z", dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {

	 	int i,j;
		unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index with the extended domain (including halo cells)
      data(iGlob) = state.z(ii);
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_real_all( ncid , zVar  , st , ct , data.data() ) , __LINE__ );

		#if SERGHEI_NC_ENABLE_MISSING_VALUE
		real missing_value = dom.MISSING_VALUE;
  	ncwrap(nc_put_att_float(ncid, zVar, SERGHEI_NC_MISSING_VALUE, SERGHEI_NC_REAL, 1, &missing_value));
		#endif

    writeStateNETCDF(state, dom, ss, par);

    ncwrap( ncmpi_close(ncid) , __LINE__ );


  }

/////////////////////////////////////////////////////////
	// NetCDF writer - requires the initialisation
  void outputNETCDF(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par, std::string dir) {

	 std::string filename;
	 filename=dir+"output.nc";



    // Create the file
    ncwrap( ncmpi_open( MPI_COMM_WORLD , filename.c_str() , NC_WRITE , MPI_INFO_NULL , &ncid ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "h" , &hVar  ) , __LINE__ );
#if SERGHEI_WRITE_HZ
    ncwrap( ncmpi_inq_varid( ncid , "h+z" , &hzVar  ) , __LINE__ );
#endif
    ncwrap( ncmpi_inq_varid( ncid , "u"      , &uVar  ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "v"      , &vVar  ) , __LINE__ );
    if(ss.inf.model){
      ncwrap( ncmpi_inq_varid( ncid , "inf"      , &infVar  ) , __LINE__ );
      ncwrap( ncmpi_inq_varid( ncid , "infVol"      , &infVolVar  ) , __LINE__ );
    }

		#if SERGHEI_MAXFLOOD > 0
    	ncwrap( ncmpi_inq_varid( ncid , "hMax" , &hMaxVar  ) , __LINE__ );
    	ncwrap( ncmpi_inq_varid( ncid , "momMax" , &momMaxVar ) , __LINE__ );
    	ncwrap( ncmpi_inq_varid( ncid , "timehMax" , &timehMaxVar  ) , __LINE__ );
		#endif

	  writeStateNETCDF(state, dom, ss, par);

    ncwrap( ncmpi_close(ncid) , __LINE__ );

  }

	// Builds the NetCDF dataset for the state
  void writeStateNETCDF(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par) {
  	ncArr data = ncArr("data",dom.ny*dom.nx);

    MPI_Offset st[3], ct[3];
	 	double timeIter[numOut+1];

	 	//write t. As the first one is written in the first iteration we should add +1
	 	for (int i=0; i<numOut+1; i++) { timeIter[i] = i*1.0;}

    st[0] = 0;
	 	ct[0] = numOut+1;

    ncwrap( ncmpi_put_vara_double_all( ncid , tVar ,  st , ct , timeIter ) , __LINE__ );

    st[0] = numOut; st[1] = par.j_beg; st[2] = par.i_beg;
    ct[0] = 1     ; ct[1] = dom.ny   ; ct[2] = dom.nx   ;

    // Write out depth
    Kokkos::parallel_for("ncwrap_h", dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
	 	int i,j;
		unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index with the extended domain (including halo cells)
      data(iGlob) = state.h(ii);
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_real_all( ncid , hVar , st , ct , data.data() ) , __LINE__ );



#if SERGHEI_WRITE_HZ
    Kokkos::parallel_for("ncwrap_h+z", dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
	 	 int i,j;
		unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index with the extended domain (including halo cells)
      data(iGlob) = state.h(ii)+state.z(ii);
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_real_all( ncid , hzVar , st , ct , data.data() ) , __LINE__ );
#endif

    // Write out x-velocity

    Kokkos::parallel_for("ncwrap_u", dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
	 	 int i,j;
		unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index with the extended domain (including halo cells)
		if(state.h(ii)>TOL12){
      	data(iGlob) = state.hu(ii)/state.h(ii);
		}else{
			data(iGlob)=0.0;
		}
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_real_all( ncid , uVar , st , ct , data.data() ) , __LINE__ );

    // Write out y-velocity
	 Kokkos::parallel_for("ncwrap_v", dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
	 	int i,j;
		unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index with the extended domain (including halo cells)
		if(state.h(ii)>TOL12){
      	data(iGlob) = state.hv(ii)/state.h(ii);
		}else{
			data(iGlob)=0.0;
		}


    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_real_all( ncid , vVar , st , ct , data.data() ) , __LINE__ );

//    realArr &infVol_p = state.ss->inf.infVol;
//    realArr &infRate_p = state.ss->inf.rate;
    if(ss.inf.model){
      // infiltration rates
	     Kokkos::parallel_for("ncwrap_inf_rate",dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
	 	      int i,j;
		      unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		      int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index with the extended domain (including halo cells)
      	  data(iGlob) = ss.inf.rate(ii);
      	  //data(iGlob) = infRate_p(ii);
        });
        Kokkos::fence();
        ncwrap( ncmpi_put_vara_real_all( ncid , infVar , st , ct , data.data() ) , __LINE__ );
      // accumulated infiltration
	     Kokkos::parallel_for("ncwrap_inf_cum", dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
	 	      int i,j;
		      unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		      int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index with the extended domain (including halo cells)
      	  //data(iGlob) = infVol_p(ii); // GPU safe
          data(iGlob) = ss.inf.infVol(ii); // GPU unsafe
        });
        Kokkos::fence();
        ncwrap( ncmpi_put_vara_real_all( ncid , infVolVar , st , ct , data.data() ) , __LINE__ );
    }

			#if SERGHEI_MAXFLOOD > 0
				#if SERGHEI_MAXFLOOD == 1		// write only at the end, otherwise writes at every time step
				if(numOut == nOut){
    			st[0] = par.j_beg; st[1] = par.i_beg;
    			ct[0] = dom.ny   ; ct[1] = dom.nx   ;
				#endif
	    	Kokkos::parallel_for("ncwrap_hMax",dom.nCellDomain , KOKKOS_LAMBDA(int iGlob) {
	 	  		int i,j;
		    	unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		    	int ii=(hc+j)*(dom.nx+2*hc)+hc+i;//index with the extended domain (including halo cells)
       		data(iGlob) = state.hMax(ii);
      	});
      	Kokkos::fence();
    		ncwrap( ncmpi_put_vara_double_all( ncid , hMaxVar , st , ct , data.data() ) , __LINE__ );

	    	Kokkos::parallel_for("ncwrap_momMax",dom.nCellDomain, KOKKOS_LAMBDA(int iGlob) {
	 	  		int i,j;
		    	unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		    	int ii=(hc+j)*(dom.nx+2*hc)+hc+i;//index with the extended domain (including halo cells)
       		data(iGlob) = state.momentumMax(ii);
      	});
      	Kokkos::fence();
				ncwrap( ncmpi_put_vara_double_all( ncid , momMaxVar , st , ct , data.data() ) , __LINE__ );

	    	Kokkos::parallel_for("ncwrap_timehMax",dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
	 	  		int i,j;
		    	unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		    	int ii=(hc+j)*(dom.nx+2*hc)+hc+i;//index with the extended domain (including halo cells)
       		data(iGlob) = state.time_hMax(ii);
      	});
      	Kokkos::fence();
	  		ncwrap( ncmpi_put_vara_double_all( ncid , timehMaxVar , st , ct , data.data() ) , __LINE__ );
			#if SERGHEI_MAXFLOOD == 1
			}
			#endif
		#endif
}

/////////////////////////////////////////

  //Error reporting routine for the PNetCDF I/O
  void ncwrap( int ierr , int line ) {
    if (ierr != NC_NOERR) {
      std::cerr<< RERROR "NetCDF Error at line: " << line <<"\n" << ncmpi_strerror(ierr) << "\n";
      exit(-1);
    }
  }

	// Writes VTK file
	void outputVTK(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par, std::string dir){

		std::string filename;
		//filename = dir+"result_."+std::to_string(par.myrank)+".vtk";
		filename = dir+"result"+std::to_string(numOut)+".vtk";
	 	real *xCoord;
  	real *yCoord;
		real *data_cpu;

		int ncells=dom.ny*dom.nx;

		int nVars=3; //z, h, (u,v)
		#if SERGHEI_WRITE_HZ
		nVars++;
		#endif
		#if SERGHEI_DEBUG_BOUNDARY
		nVars++;
		#endif
    if(ss.inf.model) nVars = nVars+2;  // 2 more variables: inf, infVol

		realArr data  = realArr("data",nVars*ncells);
		#ifdef __NVCC__
			cudaMallocHost( &data_cpu , nVars*ncells*sizeof(real) );
		#else
			data_cpu = data.data();
		#endif

		xCoord=(real*) malloc((dom.nx+1)*sizeof(real));
		yCoord=(real*) malloc((dom.ny+1)*sizeof(real));

    	// Compute x, y coordinates
    	for (int i=0; i<dom.nx+1; i++) {
   		xCoord[i] = dom.xll + ( par.i_beg + i) * dom.dx;
    	};

	 	for (int j=0; j<dom.ny+1; j++) {
      	yCoord[j] = dom.yll + dom.ny_glob*dom.dx - ( par.j_beg + j) * dom.dx;
    	};

		int nnodes=(dom.nx+1)*(dom.ny+1);
			int of_sw = 3;
			int offset_z = 0;
			int offset_h = ncells;
			int offset_hu = 2*ncells;
			int offset_hv = 3*ncells;
			int of_bc = of_sw;
			#if SERGHEI_DEBUG_BOUNDARY
			of_bc++;
			#endif
			int offset_bc = of_bc*ncells;
			int of_inf = of_bc + 1;
			int offset_infrate = of_inf * ncells;
			int offset_infVol = (of_inf+1) * ncells;
    	Kokkos::parallel_for("vtkwrap_all", ncells , KOKKOS_LAMBDA(int iGlob) {
			int i,j;
			unpackIndices(iGlob,dom.ny,dom.nx,j,i);
			int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index for the extended domain (including halo cells)
			data(iGlob+offset_z) = state.z(ii);
			data(iGlob+offset_h) = state.h(ii);
			data(iGlob+offset_hu) = state.hu(ii);
			data(iGlob+offset_hv) = state.hv(ii);
      if(ss.inf.model){
        data(iGlob+offset_infrate) = ss.inf.rate(ii);
        data(iGlob+offset_infVol) = ss.inf.infVol(ii);
      }
			#if SERGHEI_DEBUG_BOUNDARY
			data(iGlob + offset_bc) = state.isBound(ii);
			#endif
      // WARNING if you implement a new variable, you have to handle the offsets in a general case (yes, you!), to handle the possibility of different variable combinations
    	});
    	Kokkos::fence();

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
				unpackIndices(iGlob,dom.ny,dom.nx,jj,ii);
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
      }

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
    	realArr data   = realArr("data",dom.ny*dom.nx);


		//header
		filename = dir+std::to_string(par.myrank)+"result.hdr";
		std::ofstream fOutStream1(filename);
		if (fOutStream1.is_open()){
			//fOutStream << "# vtk DataFile Version 3.0.\nOutput file " << filename <<"\nASCII\nDATASET UNSTRUCTURED_GRID\n";
			fOutStream1 << "ncols " << dom.nx << std::endl;
			fOutStream1 << "nrows " << dom.ny << std::endl;
			fOutStream1 << "xllcorner " << dom.xll + par.i_beg * dom.dx << std::endl;
			fOutStream1 << "yllcorner " << dom.yll + dom.ny_glob*dom.dx - par.j_beg * dom.dx << std::endl;
			fOutStream1 << "cellsize " << dom.dx << std::endl;
			fOutStream1 << "nodata_value -9999" << std::endl;
			fOutStream1 << "byteorder msbfirst" << std::endl;
			fOutStream1.close();
		}

		Kokkos::fence();



		//z
    	Kokkos::parallel_for("binwrap_init_z", dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
	 		int i,j;
			unpackIndices(iGlob,dom.ny,dom.nx,j,i);
			int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index with the extended domain (including halo cells)
      	data(iGlob) = state.z(ii);
    	});

		Kokkos::fence();




		filename = dir+std::to_string(par.myrank)+"elevation.bin";
		std::ofstream fOutStream2(filename.c_str(), std::ios::binary);
		fOutStream2.write((char*)&data[0], dom.ny*dom.nx * sizeof(real));
		fOutStream2.close();

		Kokkos::fence();


		//h
    	Kokkos::parallel_for("binwrap_init_h", dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
	 		int i,j;
			unpackIndices(iGlob,dom.ny,dom.nx,j,i);
			int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index with the extended domain (including halo cells)
      	data(iGlob) = state.h(ii);
    	});
		Kokkos::fence();



		filename = dir+std::to_string(par.myrank)+"result"+std::to_string(numOut)+".bin";
		std::ofstream fOutStream3(filename.c_str(), std::ios::binary);
		fOutStream3.write((char*)&data[0], dom.ny*dom.nx * sizeof(real));
		fOutStream3.close();

		Kokkos::fence();

		if (par.nranks > 1){
			MPI_Barrier(MPI_COMM_WORLD);
		}

	}

	// Writes a binary raster file
	void outputBIN(const State &state, Domain const &dom, Parallel const &par, std::string dir){

		std::string filename;
    	realArr data   = realArr("data",dom.ny*dom.nx);

    Kokkos::parallel_for("binwrap_h", dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
	 	int i,j;
		unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index with the extended domain (including halo cells)
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
		fOutStream.write((char*)&data[0], dom.ny*dom.nx * sizeof(real));
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



	int _ncellsBC = 0;
	bint.integrate(extbc,dom,1);
	_ncellsBC = bint.ncellsBC;




#if SERGHEI_DEBUG_BOUNDARY
	if (par.masterproc)
	  {
	    std::cerr << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << "boundary cells (reduced): " << _ncellsBC << "\n";
	    std::cerr << GGD  << GRAY << __PRETTY_FUNCTION__ << RESET << "inflow discharge (integrated) " << bint.inflowDischargeG << "\n";
	    std::cerr << GGD  << GRAY << __PRETTY_FUNCTION__ << RESET << "inflow accumulated (integrated) " << bint.inflowAccumulatedG << "\n";
	    std::cerr << GGD  << GRAY << __PRETTY_FUNCTION__ << RESET << "outflow discharge (integrated) " << bint.outflowDischargeG << "\n";
	    std::cerr << GGD  << GRAY << __PRETTY_FUNCTION__ << RESET << "outflow accumulated (integrated) " << bint.outflowAccumulatedG << std::endl;
	  }
#endif



	if (_ncellsBC)	  {
	    domainOutputFile << "BoundaryInflow ";
	    domainOutputFile << "BoundaryInflowAccum ";
	    domainOutputFile << "BoundaryOutflow ";
	    domainOutputFile << "BoundaryOutflowAccum ";
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
    writeTimeSeries (state, dom, par, sint, bint);

    numObs++;

    return 1;

  }


  // Writes time series files
  void writeTimeSeries (const State &state, Domain const &dom, Parallel const &par, surfaceIntegrator &sint, boundaryIntegrator &bint){
    timer.reset();
    if(par.masterproc) writeDomainTimeSeries(state,dom,par,sint,bint);
    numObs++;
    dom.timers.out += timer.seconds();
  }

  void
  writeDomainTimeSeries (State const &state, Domain const &dom,
			 Parallel const &par, surfaceIntegrator &sint,
			 boundaryIntegrator &bint){

    // Write the data
    std::cout.precision(OUTPUT_PRECISION);
    domainOutputFile << std::scientific << dom.etime << " ";
    domainOutputFile << std::scientific << sint.surfaceVolumeG << " ";

    if (bint.ncellsBC){
	   domainOutputFile << std::scientific << bint.inflowDischargeG << " ";
	   domainOutputFile << std::scientific << bint.inflowAccumulatedG << " ";
      domainOutputFile << std::scientific << bint.outflowDischargeG << " ";
	   domainOutputFile << std::scientific << bint.outflowAccumulatedG << " ";
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
      logFile << "nCellDomain : " << dom.nCellDomainGlobal << std::endl;
      logFile << "nCellValid : " << dom.nCellValidGlobal << std::endl;

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
      ratio = dom.timers.integrate / dom.timers.total;
      logFile << "integrateTime : " << dom.timers.integrate << " : " << ratio << std::endl;
      ratio = dom.timers.dt / dom.timers.total;
      logFile << "dtComputeTime : " << dom.timers.dt << " : " << ratio << std::endl;
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
          writeTimerRank(par,dom.timers.dt,"dtComputeTime");
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
              logFile << "CPU: " << CPUBrandString << std::endl;
              #endif
              logFile << "nTasks : " << par.nranks << std::endl;

              // write compilation setup
      // logFile << "\nSERGHEI_GIT_VERSION " << SERGHEI_GIT_VERSION << std::endl;
              logFile << "\n-------------------------\nMODEL COMPONENT SETUP" <<std::endl;
              logFile << "SERGHEI_TOOLS " << SERGHEI_TOOLS << std::endl;
              logFile << "SERGHEI_SUBSURFACE_MODEL " << SERGHEI_SUBSURFACE_MODEL << std::endl;
              logFile << "SERGHEI_PARTICLE_TRACKING " << SERGHEI_PARTICLE_TRACKING << std::endl;
              logFile << "SERGHEI_VEGETATION_MODEL " << SERGHEI_VEGETATION_MODEL << std::endl;
              logFile << "SERGHEI_FRICTION_MODEL " << SERGHEI_FRICTION_MODEL << std::endl;
              logFile << "SERGHEI_MAXFLOOD " << SERGHEI_MAXFLOOD << std::endl;
      logFile << "SERGHEI_REAL " << SERGHEI_REAL << std::endl;
      logFile << "SERGHEI_NC_MODE " << SERGHEI_NC_MODE << std::endl;
      logFile << "SERGHEI_NC_REAL " ;
        if(SERGHEI_NC_REAL == NC_FLOAT) logFile << "NC_FLOAT";
        if(SERGHEI_NC_REAL == NC_DOUBLE) logFile << "NC_DOUBLE";
        logFile << std::endl;
      logFile << "SERGHEI_WRITE_HZ " << SERGHEI_WRITE_HZ << std::endl;

      logFile << "\n-------------------------\nDEBUG FLAGS" <<std::endl;
      logFile << "SERGHEI_DEBUG_PARALLEL_DECOMPOSITION " << SERGHEI_DEBUG_PARALLEL_DECOMPOSITION << std::endl;
      logFile << "SERGHEI_DEBUG_WORKFLOW " << SERGHEI_DEBUG_WORKFLOW << std::endl;
      logFile << "SERGHEI_DEBUG_KOKKOS_SETUP " << SERGHEI_DEBUG_KOKKOS_SETUP << std::endl;
      logFile << "SERGHEI_DEBUG_BOUNDARY " << SERGHEI_DEBUG_BOUNDARY << std::endl;
      logFile << "SERGHEI_DEBUG_DT " << SERGHEI_DEBUG_DT << std::endl;
      logFile << "SERGHEI_DEBUG_TOOLS " << SERGHEI_DEBUG_TOOLS << std::endl;
      logFile << "SERGHEI_DEBUG_MASS_CONS " << SERGHEI_DEBUG_MASS_CONS << std::endl;
      logFile << "SERGHEI_DEBUG_INFILTRATION " << SERGHEI_DEBUG_INFILTRATION << std::endl;
      logFile << "SERGHEI_DEBUG_MPI " << SERGHEI_DEBUG_MPI << std::endl;
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
        realArr data   = realArr("data",gdom.nCellDomain);
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
        nxhalo = gdom.nx + 2*haloc;
        nyhalo = gdom.ny + 2*haloc;
        st[0] = 0;          st[1] = par.j_beg;  st[2] = par.i_beg;
        ct[0] = gdom.nz;    ct[1] = gdom.ny;    ct[2] = gdom.nx;
        Kokkos::parallel_for( gdom.nCellDomain , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            unpackIndices(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (haloc+kk)*nxhalo*nyhalo + (haloc+jj)*nxhalo + ii + haloc;
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
        realArr data = realArr("data",gdom.nCellDomain);
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

        // Kokkos::parallel_for(gdom.nCellDomain, kernel_output<dspace>(gw.h, data, gdom));
        Kokkos::parallel_for( gdom.nCellDomain , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            unpackIndices(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*nxhalo*nyhalo + (hc+jj)*nxhalo + ii + hc;
            data(idom) = gw.h(iGlob,1);
        });
        Kokkos::fence();
        ncwrap( ncmpi_put_vara_double_all( ncid , hdVar , st , ct , data.data() ) , __LINE__ );

        // Kokkos::parallel_for(gdom.nCellDomain, kernel_output<dspace>(gw.wc, data, gdom));
        Kokkos::parallel_for( gdom.nCellDomain , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, iGlob;
            unpackIndices(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
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
          int nxhalo = gdom.nx+2*haloc;
          int nyhalo = gdom.ny+2*haloc;
          unpackIndices(idom,gdom.nz,gdom.ny,gdom.nx,k,j,i);
          int iGlob = (haloc+k)*nxhalo*nyhalo + (haloc+j)*nxhalo + i + haloc;
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
              unpackIndices(iGlob,gdom.nz,gdom.ny,gdom.nx,kk,jj,ii);
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
  int writeSubTimeSeriesIni (const GwState &gw, GwDomain const &gdom, Domain const &dom, Parallel const &par, std::string dir){

        numObs = 0;
        std::string filename = dir + "SubsurfaceTimeSeries.out";
        SubsurfaceOutputFile.open (filename);

        if (SubsurfaceOutputFile.is_open ()){
            // Write the header
            SubsurfaceOutputFile << "Time ";
            SubsurfaceOutputFile << "SubSurfaceVolume ";
            SubsurfaceOutputFile << "ExchangeVolume ";
            SubsurfaceOutputFile << std::endl;
        }
        else
        {std::cerr << RERROR "Could not create domainTimeSeries.out file" << std::endl; return 0;}

        if(par.masterproc)  {writeSubsurfaceTimeSeries (gw, gdom, dom, par);}
        numObs++;
        return 1;
  }



  void writeSubsurfaceTimeSeries (GwState const &gw, GwDomain const &gdom, Domain const &dom, Parallel const &par){
    // Write the data
    std::cout.precision(OUTPUT_PRECISION);
    SubsurfaceOutputFile << std::scientific << dom.etime << " ";
    SubsurfaceOutputFile << std::scientific << gw.Vtot << " ";
    SubsurfaceOutputFile << std::scientific << gw.Vexch << " ";
    SubsurfaceOutputFile << std::endl;

  }
  #endif


};


#endif
