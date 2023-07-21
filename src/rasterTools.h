#ifndef __RASTER_H__
#define __RASTER_H__

#include <string>
#include "define.h"
#include "Domain.h"
#include "Parallel.h"

#define SERGHEI_DATA_CHECK_POSITIVE 1

int readRasterField(const std::string fNameIn, const Domain &dom, const Parallel &par, realArr data) {
    std::cout << BDASH << "Reading raster file " << fNameIn << std::endl;
    realArr tmpVar=realArr("var", dom.ny_glob*dom.nx_glob );
    std::ifstream fInStream(fNameIn);
    std::string line;

	  int tnx,tny;
	  real txll,tyll,tdx;
	  real nodata;

	  tnx=-999;
	  tny=-999;
	  txll=-999;
	  tyll=-999;
	  tdx=-999;
	  nodata = 123456789;

	  int ndata=dom.ny_glob*dom.nx_glob;
    std::string str;
	  if (fInStream.is_open()){
      std::getline(fInStream,str,' ');
      std::getline(fInStream,str);
      std::stringstream(str) >> tnx;

      std::getline(fInStream,str,' ');
      std::getline(fInStream,str);
      std::stringstream(str) >> tny;

      std::getline(fInStream,str,' ');
      std::getline(fInStream,str);
      std::stringstream(str) >> txll;

      std::getline(fInStream,str,' ');
      std::getline(fInStream,str);
      std::stringstream(str) >> tyll;

      std::getline(fInStream,str,' ');
      std::getline(fInStream,str);
      std::stringstream(str) >> tdx;

      std::getline(fInStream,str,' ');
      std::getline(fInStream,str);
      std::stringstream(str) >> nodata;

		  //compare the values t* with the DEM raster size, to check if we are using the same values, otherwise error
		  if(dom.ny_glob !=tny || dom.nx_glob !=tnx || dom.xll !=txll || dom.yll !=tyll || dom.dx !=tdx){
			  if(par.masterproc){
				  std::cerr << RERROR "Parameters in " << fNameIn << " don't match DEM file parameters." << std::endl;
        }
			  return 0;
		  }

		  real tmp;

      for (int ii=0; ii<ndata; ii++) {
			  if (!fInStream.fail() && !fInStream.eof()){
				  fInStream >> tmp;
          /*
				  if(tmp<0.0){
				    if(par.masterproc){
					    std::cerr << tmp << std::endl;
					    std::cerr<< RERROR "There are some negative depth values. Unable to continue\n";
					return 0;
					}
				}
          */
				tmpVar(ii)=tmp;
			  }else{
		      if(par.masterproc){
            std::cerr<< RERROR "Error reading " << fNameIn << ". Not enough data" << std::endl;
				    return 0;
          }
			  }
		  }
		  fInStream.close();
    }

    Kokkos::parallel_for("generic_raster_data", dom.ny*dom.nx , KOKKOS_LAMBDA (int iGlob) {
 		  int i,j;
		  unpackIndices(iGlob,dom.ny,dom.nx,j,i);
		  int ii1=(hc+j)*(dom.nx+2*hc)+hc+i;//index for the extended domain (including halo cells)
		  int ii2=(par.j_beg+j)*(dom.nx_glob)+par.i_beg+i;//index for the subdomain (par.j_beg+j,par.i_beg+i)
		  data(ii1)=tmpVar(ii2);
	  });

 	  return 1;
  }

#endif

