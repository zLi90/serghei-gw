
#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "define.h"
#include "Parallel.h"
#include "geometry.h"
#include "Indexing.h"

KOKKOS_INLINE_FUNCTION int getHaloExtension(const int i, const int j, const int nx){
  return( (hc+j)*(nx+2*hc)+hc+i ); //index for the extended domain (including halo cells)
};

class Domain {

public:

  real cfl;
  double simLength;

  double etime;

  real dt;

  SergheiTimers mutable timers;

  //raster variables
  int nx_glob;
  int ny_glob;
  int nx; // physical number of cells in x-direction
  int ny; // physical number of cells in y-direction
  int ncells;
  int nCellDomain;
  int nCellValid;
  real xll; // southwest corner x-coordinate
  real yll; // southwest corner y-coordinate
  real dx;  // resolution

  int iE,iW,iS,iN; //flag to see if the subdomain touch with either a East, West, South or North boundaries

  // global (reduced) variables
  real areaGlobal;
  int nCellValidGlobal;
  int nCellDomainGlobal;

  // other variables
  int BCtype;
  int isRain = 0;

  int nIter;
  int countIterDt;

  real area;
  int id;  // subdomain ID

  geometry::point extent[2];

  void domainArea(){
       area = cellArea()*nCellValid;
   };

  KOKKOS_INLINE_FUNCTION geometry::point getCellCenter(int i, int j){
    geometry::point p;
    p(_X) = i*dx + extent[0](_X);
    p(_Y) = j*dx + extent[0](_Y);
    return(p);
  }

  KOKKOS_INLINE_FUNCTION geometry::point getCellCenter(int iGlob){
    int i,j;
    unpackIndices(iGlob,ny,nx,j,i);
    return(getCellCenter(i,j));
  }

  KOKKOS_INLINE_FUNCTION real cellArea(){
    return(dx*dx);  // WARNING UCM
  }

  KOKKOS_INLINE_FUNCTION void getMatrixIndicesForPoint(const geometry::point &p, int &i, int &j) const{
    // get i,j coordinates of the cell containing the point
    // this only works on a raster-order grid, where j is zero at NORTH boundary
    i = floor((p(_X)-extent[0](_X))/dx);
    j = floor((extent[1](_Y)-p(_Y))/dx);
  };

  KOKKOS_INLINE_FUNCTION int getCellForPoint(const geometry::point &p) const{
    int i,j,iGlob;
    getMatrixIndicesForPoint(p,i,j);
    if(i < 0 || i >= nx || j < 0 || j >= ny){
      iGlob = -1;
    }else{
      iGlob = packIndices(ny,nx,j,i);
    }
    // note: domain decomposition is handled by the i,j coords
    return(iGlob);
  };

  KOKKOS_INLINE_FUNCTION int getIndex(int iGlob) const{
    int i,j;
    unpackIndices(iGlob,ny,nx,j,i);
    int ii=(hc+j)*(nx+2*hc)+hc+i; //index for the extended domain (including halo cells)
    return(ii);
  };

  KOKKOS_INLINE_FUNCTION int getIndexForPoint(const geometry::point &p) const{
    int i,j,ii;
    getMatrixIndicesForPoint(p,i,j);
    if(i < 0 || i >= nx || j < 0 || j >= ny){
      ii = -1;
    }else{
      ii  = getHaloExtension(i,j,nx);
    }
    return(ii);
  };

void initialise() {
   	// Initialize the time
   	etime = 0;
   	nIter = 0;
		countIterDt=0;

    // physical cells onlys
    nCellDomain = nx*ny; // WARNING UCM
    // physical cells + halo cells
		ncells=(ny+2*hc)*(nx+2*hc);
  };

void getStatistics(){
    domainArea();
    MPI_Allreduce(&area, &areaGlobal, 1, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nCellValid, &nCellValidGlobal, 1, MPI_INT , MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nCellDomain, &nCellDomainGlobal, 1, MPI_INT , MPI_SUM, MPI_COMM_WORLD);
};

  int buildDomainDecomposition(Parallel &par) {

    int ierr=1;

    if (par.nranks != par.nproc_x*par.nproc_y) {
      std::cerr << RERROR "ERROR: nproc_x*nproc_y != nranks" << std::endl;
      std::cerr << RERROR << par.nproc_x << " " << par.nproc_y << " " << par.nranks << std::endl;
      exit(-1);
    }

    //Get my x and y process grid ID
    par.px = par.myrank % par.nproc_x;
    par.py = par.myrank / par.nproc_x;

    // name subdomain
    id = par.myrank;

    //Get my beginning and ending global indices
    double nper;
    nper = ((double) nx_glob)/par.nproc_x;
    par.i_beg = (long) round( nper* par.px    );
    par.i_end = (long) round( nper*(par.px+1) )-1;
    nper = ((double) ny_glob)/par.nproc_y;
    par.j_beg = (long) round( nper* par.py    );
    par.j_end = (long) round( nper*(par.py+1) )-1;

    // Determine my number of grid cells
    nx = par.i_end - par.i_beg + 1;
    ny = par.j_end - par.j_beg + 1;

    // Determine my extent, point 0 is SW, point 1 is NE (standard cartesian)
    extent[0](_X) = xll + par.i_beg*dx;
    extent[0](_Y) = yll + ny_glob*dx - (par.j_end+1)*dx;
    extent[1](_X) = xll + (par.i_end+1)*dx;
    extent[1](_Y) = yll + ny_glob*dx - (par.j_beg)*dx;


    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        int pxloc = par.px+i-1;
        if (pxloc < 0            ) pxloc = pxloc + par.nproc_x;
        if (pxloc > par.nproc_x-1) pxloc = pxloc - par.nproc_x;
        int pyloc = par.py+j-1;
        if (pyloc < 0            ) pyloc = pyloc + par.nproc_y;
        if (pyloc > par.nproc_y-1) pyloc = pyloc - par.nproc_y;
        par.neigh(j,i) = pyloc * par.nproc_x + pxloc;
      }
    }

    // Debug output for the parallel decomposition
    #if DEBUG_PARALLEL_DECOMPOSITION
      for (int rr=0; rr < par.nranks; rr++) {
        if (rr == par.myrank) {
          std::cerr << GGD "Hello! My Rank is: " << par.myrank << "\n";
          std::cerr << GGD "My proc grid ID is: " << par.px << " , " << par.py << "\n";
          std::cerr << GGD "I have: " << nx << " x " << ny << " grid cells" << "\n";
          std::cerr << GGD "I start at index: " << par.i_beg << " x " << par.j_beg << "\n";
          std::cerr << GGD << "My extent points are " ;
          std::cerr << "(" << extent[0](_X) << "," << extent[0](_Y) << ") ";
          std::cerr << "(" << extent[1](_X) << "," << extent[1](_Y) << ") " << std::endl;
          std::cerr << GGD "My neighbor matrix is:\n";
          for (int j = 2; j >= 0; j--) {
            for (int i = 0; i < 3; i++) {
              std::cerr << std::setw(6) << par.neigh(j,i) << " ";
            }
            std::cerr << "\n";
          }
          std::cerr << "\n";
          }
          ierr = MPI_Barrier(MPI_COMM_WORLD);
      }
      ierr = MPI_Barrier(MPI_COMM_WORLD);
    #endif
	   return ierr;
	};

};

KOKKOS_INLINE_FUNCTION int getIndex(const int iGlob, Domain const &dom){
  int i,j;
  unpackIndices(iGlob,dom.ny,dom.nx,j,i);
  return ( getHaloExtension(i,j,dom.nx));
//  int ii=(hc+j)*(dom.nx+2*hc)+hc+i; //index for the extended domain (including halo cells)
//  return(ii);
};



#endif
