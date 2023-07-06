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

class FileIO {

protected:

  int ncid;
  int tDim, xDim, yDim, zDim;
  int tVar, xVar, yVar, hVar, hzVar, uVar, vVar, zVar, z3Var, hdVar, wcVar;
  int infVar,infVolVar;
  std::ofstream domainOutputFile;
  std::ofstream SubsurfaceOutputFile;
  std::ofstream logFile;

public:

  real outFreq;
  int outFormat;
  int numOut;

  real obsFreq;
  int numObs = 0;
  int nScreen = 100; //if this value is not specified, the information is displayed every 1000 iterations

// Writes spatial fields for initial state
	void outputIni(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par, std::string dir){
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
	}


	// NetCDF initialiser and writer
  void outputInitNETCDF(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par, std::string dir) {
    int dimids[3];
    MPI_Offset st[3], ct[3];
    realArr xCoord = realArr("xCoord",dom.nx);
    realArr yCoord = realArr("yCoord",dom.ny);
    realArr data   = realArr("data",dom.ny*dom.nx);
	static char title[] = "seconds" ;
	std::string filename;

	filename=dir+"output.nc";



    // Create the file
    ncwrap( ncmpi_create( MPI_COMM_WORLD , filename.c_str() , NC_CLOBBER , MPI_INFO_NULL , &ncid ) , __LINE__ );

    // Create the dimensions
    ncwrap( ncmpi_def_dim( ncid , "t" , (MPI_Offset) NC_UNLIMITED , &tDim ) , __LINE__ );
    ncwrap( ncmpi_def_dim( ncid , "x" , (MPI_Offset) dom.nx_glob  , &xDim ) , __LINE__ );
    ncwrap( ncmpi_def_dim( ncid , "y" , (MPI_Offset) dom.ny_glob  , &yDim ) , __LINE__ );
    // Create the variables
    dimids[0] = tDim;
    ncwrap( ncmpi_def_var( ncid , "t"      , NC_DOUBLE , 1 , dimids , &tVar ) , __LINE__ );
	 ncwrap( ncmpi_put_att_text (ncid, tVar, "units",strlen(title), title), __LINE__ );
    dimids[0] = xDim;
    ncwrap( ncmpi_def_var( ncid , "x"      , NC_DOUBLE , 1 , dimids , &xVar ) , __LINE__ );
    dimids[0] = yDim;
    ncwrap( ncmpi_def_var( ncid , "y"      , NC_DOUBLE , 1 , dimids , &yVar ) , __LINE__ );
    dimids[0] = tDim; dimids[1] = yDim; dimids[2] = xDim;
    ncwrap( ncmpi_def_var( ncid , "h" , NC_DOUBLE , 3 , dimids , &hVar  ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "h+z" , NC_DOUBLE , 3 , dimids , &hzVar  ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "u"      , NC_DOUBLE , 3 , dimids , &uVar  ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "v"      , NC_DOUBLE , 3 , dimids , &vVar  ) , __LINE__ );
    if(ss.inf.model){
      ncwrap( ncmpi_def_var( ncid , "inf" , NC_DOUBLE , 3 , dimids , &infVar  ) , __LINE__ );
      ncwrap( ncmpi_def_var( ncid , "infVol" , NC_DOUBLE , 3 , dimids , &infVolVar  ) , __LINE__ );
    }
    dimids[0] = yDim; dimids[1] = xDim;
    ncwrap( ncmpi_def_var( ncid , "z"    , NC_DOUBLE , 2 , dimids , &zVar ) , __LINE__ );


    // End "define" mode
    ncwrap( ncmpi_enddef( ncid ) , __LINE__ );

    // Compute x, y coordinates
    // for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.nx , KOKKOS_LAMBDA(int i) {
      xCoord(i) = dom.xll + ( par.i_beg + i + 0.5) * dom.dx;
    });
    // for (int j=0; j<dom.ny; j++) {
    Kokkos::parallel_for( dom.ny , KOKKOS_LAMBDA(int j) {
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

    // for (int j=0; j<dom.ny; j++) {
    //   for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
	 	int i,j;
		unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index with the extended domain (including halo cells)
      data(iGlob) = state.z(ii);
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_double_all( ncid , zVar  , st , ct , data.data() ) , __LINE__ );

    writeStateNETCDF(state, dom, ss, par);

    ncwrap( ncmpi_close(ncid) , __LINE__ );


  }

	// NetCDF writer - requires the initialisation
  void outputNETCDF(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par, std::string dir) {

	 std::string filename;
	 filename=dir+"output.nc";



    // Create the file
    ncwrap( ncmpi_open( MPI_COMM_WORLD , filename.c_str() , NC_WRITE , MPI_INFO_NULL , &ncid ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "h" , &hVar  ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "h+z" , &hzVar  ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "u"      , &uVar  ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "v"      , &vVar  ) , __LINE__ );
    if(ss.inf.model){
      ncwrap( ncmpi_inq_varid( ncid , "inf"      , &infVar  ) , __LINE__ );
      ncwrap( ncmpi_inq_varid( ncid , "infVol"      , &infVolVar  ) , __LINE__ );
    }

    writeStateNETCDF(state, dom, ss, par);

    ncwrap( ncmpi_close(ncid) , __LINE__ );

  }

	// Builds the NetCDF dataset for the state
  void writeStateNETCDF(const State &state, Domain const &dom, SourceSinkData &ss, Parallel const &par) {
  	realArr data = realArr("data",dom.ny*dom.nx);
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
    // for (int j=0; j<dom.ny; j++) {
    //   for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
	 	int i,j;
		unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index with the extended domain (including halo cells)
      data(iGlob) = state.h(ii);
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_double_all( ncid , hVar , st , ct , data.data() ) , __LINE__ );

    Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
	 	 int i,j;
		unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index with the extended domain (including halo cells)
      data(iGlob) = state.h(ii)+state.z(ii);
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_double_all( ncid , hzVar , st , ct , data.data() ) , __LINE__ );

    // Write out x-velocity

    Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
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
    ncwrap( ncmpi_put_vara_double_all( ncid , uVar , st , ct , data.data() ) , __LINE__ );

    // Write out y-velocity
	 Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
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
    ncwrap( ncmpi_put_vara_double_all( ncid , vVar , st , ct , data.data() ) , __LINE__ );

//    realArr &infVol_p = state.ss->inf.infVol;
//    realArr &infRate_p = state.ss->inf.rate;
    if(ss.inf.model){
      // infiltration rates
	     Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
	 	      int i,j;
		      unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		      int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index with the extended domain (including halo cells)
      	  data(iGlob) = ss.inf.rate(ii);
      	  //data(iGlob) = infRate_p(ii);
        });
        Kokkos::fence();
        ncwrap( ncmpi_put_vara_double_all( ncid , infVar , st , ct , data.data() ) , __LINE__ );
      // accumulated infiltration
	     Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
	 	      int i,j;
		      unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		      int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index with the extended domain (including halo cells)
      	  //data(iGlob) = infVol_p(ii); // GPU safe
          data(iGlob) = ss.inf.infVol(ii); // GPU unsafe
        });
        Kokkos::fence();
        ncwrap( ncmpi_put_vara_double_all( ncid , infVolVar , st , ct , data.data() ) , __LINE__ );
    }


}


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

		int nVars=4; //4 variables to write: z, h, hu, hv

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

    	Kokkos::parallel_for( ncells , KOKKOS_LAMBDA(int iGlob) {
			int i,j;
			unpackIndices(iGlob,dom.ny,dom.nx,j,i);
			int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index for the extended domain (including halo cells)
			data(iGlob) = state.z(ii);
			data(iGlob+ncells) = state.h(ii);
			data(iGlob+2*ncells) = state.hu(ii);
			data(iGlob+3*ncells) = state.hv(ii);
      if(ss.inf.model){
        data(iGlob+4*ncells) = ss.inf.rate(ii);
        data(iGlob+5*ncells) = ss.inf.infVol(ii);
      }
      // WARNING if you implement a new variable, you have to handle the iGlob+4 in a general case (yes, you!), to handle the possibility of different variable combinations
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

			fOutStream << "POINTS " << nnodes  << " double\n";
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
			fOutStream << "SCALARS z double\n";
			fOutStream << "LOOKUP_TABLE default\n";
			for(iGlob=0;iGlob<ncells;iGlob++){
				fOutStream << std::setprecision(9) << data_cpu[iGlob] << "\n";
			}

			fOutStream << "SCALARS h double\n";
			fOutStream << "LOOKUP_TABLE default\n";
			for(iGlob=0;iGlob<ncells;iGlob++){
				fOutStream << std::setprecision(9) << data_cpu[iGlob+ncells]  << "\n";
			}

			fOutStream << "SCALARS h+z double\n";
			fOutStream << "LOOKUP_TABLE default\n";
			for(iGlob=0;iGlob<ncells;iGlob++){
				fOutStream << std::setprecision(9) << data_cpu[iGlob] + data_cpu[iGlob+ncells] << "\n";
			}


			fOutStream << "VECTORS velocity double\n";
			for(iGlob=0;iGlob<ncells;iGlob++){
				if(data_cpu[iGlob+ncells]>TOL12){
					fOutStream << std::setprecision(9) << data_cpu[iGlob+2*ncells]/data_cpu[iGlob+ncells] << " " << data_cpu[iGlob+3*ncells]/data_cpu[iGlob+ncells] << " 0.0\n";
				}else{
					fOutStream << "0.0 0.0 0.0\n" ;
				}
			}

      if(ss.inf.model){
        fOutStream << "SCALARS inf double\n";
        fOutStream << "LOOKUP_TABLE default\n";
        for(iGlob=0;iGlob<ncells;iGlob++){
          fOutStream << std::setprecision(9) << data_cpu[iGlob+4*ncells]  << "\n";
        }

        fOutStream << "SCALARS infVol double\n";
        fOutStream << "LOOKUP_TABLE default\n";
        for(iGlob=0;iGlob<ncells;iGlob++){
          fOutStream << std::setprecision(9) << data_cpu[iGlob+5*ncells]  << "\n";
        }
      }
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
    	Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
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
    	Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
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

    Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
	 	int i,j;
		unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		int ii=(haloc+j)*(dom.nx+2*haloc)+haloc+i;//index with the extended domain (including halo cells)
      data(iGlob) = state.h(ii);
    });

    Kokkos::fence();

		//MPI_Gatherv(&data[0],cur_proc_data_size, MPI_DOUBLE, total_data_arr, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		/*if (rank_ == 0){
			MPI_Gatherv(state.h.get_address_at(0, 0), cur_proc_data_size, MPI_DOUBLE, total_data_arr, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}else{
			MPI_Gatherv(state.h.get_address_at(1, 0), cur_proc_data_size, MPI_DOUBLE, total_data_arr, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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



#if DEBUG_BOUNDARY

	std::cerr << GGD "Hello from rank: " << par.myrank << "/" << par.nranks << "\n";
	for (int i = 0; i < extbc.size(); i++) {
	  std::cerr << GGD "boundary cells: " << extbc[i].ncellsBC << "\n";
	  std::cerr << GGD "inflow discharge: " << extbc[i].inflowDischarge << "\n";
	  std::cerr << GGD "inflow accumulated: " << extbc[i].inflowAccumulated << "\n";
	  std::cerr << GGD "outflow discharge: " << extbc[i].outflowDischarge << "\n";
	  std::cerr << GGD "outflow accumulated: " << extbc[i].outflowAccumulated << "\n";
	}

#endif



	int _ncellsBC = 0;
	bint.integrate(extbc,1);
	_ncellsBC = bint.ncellsBC;




#if DEBUG_BOUNDARY
	if (par.masterproc)
	  {
	    std::cerr << GGD "boundary cells (reduced): " << _ncellsBC << "\n";
	    std::cerr << GGD "inflow discharge (integrated) " << bint.inflowDischargeG << "\n";
	    std::cerr << GGD "inflow accumulated (integrated) " << bint.inflowAccumulatedG << "\n";
	    std::cerr << GGD "outflow discharge (integrated) " << bint.outflowDischargeG << "\n";
	    std::cerr << GGD "outflow accumulated (integrated) " << bint.outflowAccumulatedG << std::endl;
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
  void writeTimeSeries (const State &state, Domain const &dom,
		   Parallel const &par, surfaceIntegrator &sint,
		   boundaryIntegrator &bint){
         if(par.masterproc) writeDomainTimeSeries(state,dom,par,sint,bint);
    numObs++;
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


  void writeLogFile(Domain const &dom, SergheiTimers &timers, std::string dir){
    std::string filename = dir + "log.out";
    logFile.open(filename);
    if (domainOutputFile.is_open()){
      logFile << "DomainArea : " << dom.area << std::endl;
      logFile << "nCellDomain : " << dom.nCellDomain << std::endl;
      real ratio = timers.total/timers.total;
      logFile << "runTime : " << timers.total << " : " << ratio << std::endl;
      real sim = timers.total-timers.Tinit;
      ratio = sim/timers.total;
      logFile << "simTime : " << sim << " : " << ratio << std::endl;
      ratio = timers.Tinit/timers.total;
      logFile << "initTime : " << timers.Tinit << " : " <<ratio << std::endl;
      ratio = timers.Tout/timers.total;
      logFile << "outputTime : " << timers.Tout << " : " << ratio <<std::endl;
      ratio = timers.Tsweflux/timers.total;
      logFile << "sweFluxTime : " << timers.Tsweflux << " : " << ratio << std::endl;
      ratio = timers.Traininf/timers.total;
      logFile << "rainInfTime : " << timers.Traininf << " : " << ratio << std::endl;
      ratio = timers.Tswe/timers.total;
      logFile << "sweNotFluxTime : " << timers.Tswe << " : " << ratio << std::endl;
      ratio = timers.Texchange / timers.total;
      logFile << "exchangeTime : " << timers.Texchange << " : " << ratio << std::endl;
      ratio = timers.Tintegrate / timers.total;
      logFile << "integrateTime : " << timers.Tintegrate << " : " << ratio << std::endl;


      auto nowtime = std::chrono::system_clock::now();
      std::time_t now_time = std::chrono::system_clock::to_time_t(nowtime);

      logFile << "DateTime : " << std::ctime(&now_time) << std::endl;

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
      // uint *CPUInfo = reinterpret_cast<uint*>(CPUBrandString.data());
      // for (uint i=0; i<3; i++) __cpuid(0x80000002+i, CPUInfo[i*4+0], CPUInfo[i*4+1], CPUInfo[i*4+2], CPUInfo[i*4+3]);
      CPUBrandString.assign(CPUBrandString.data()); // correct null terminator
      logFile << "CPU: " << CPUBrandString << std::endl;
      #endif

      // write compilation setup
      logFile << "\n-------------------------\nMODEL COMPONENT SETUP" <<std::endl;
      logFile << "SERGHEI_TOOLS " << SERGHEI_TOOLS << std::endl;
      //logFile << "SERGHEI_SUBSURFACE_MODEL " << SERGHEI_SUBSURFACE_MODEL << std::endl;
      logFile << "SERGHEI_PARTICLE_TRACKING " << SERGHEI_PARTICLE_TRACKING << std::endl;

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
