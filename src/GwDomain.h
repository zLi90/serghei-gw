#ifndef _GW_DOMAIN_H_
#define _GW_DOMAIN_H_

#include "define.h"
#include "Domain.h"

class GwDomain : public Domain {

public:
    // Time stepping options
    real dt_init, dt_max, dt, dtOld;
    // Subsurface domain dimensions
    real thickH, topZ, zll, dx, dy, dz_multiplier;
    int nz, nz_glob, nhalo, nxhc, nyhc, nzhc;
    int nCellSw, nCellSwMem;
    // Domain properties
    int nSoilID, aev, hmin;
    int hasRoot, hasET;
    // Numerical scheme
    int gw_scheme;
    bool async;
    // Kokkos views
    realArr z, dz, sinx, cosx, siny, cosy, rainRate, evapRate, etpmRate;
    intArr isnodata;

    // other variables
    int BCtype;
    realArr globalBuffer;
    geometry::point extent[2];

    KOKKOS_INLINE_FUNCTION void unpackIndicesGw(int const iGlob, int nz, int ny, int nx, int &k, int &j, int &i) const{
      unpackIndicesUniformGrid(iGlob,nz,ny,nx,k,j,i);
    };
    KOKKOS_INLINE_FUNCTION void unpackIndices(int const iGlob, int &k, int &j, int &i) const{
        unpackIndicesUniformGrid(iGlob, nz, ny, nx, k, j, i);
    };

    KOKKOS_INLINE_FUNCTION void unpackIndicesHalo(int const iGlob, int &k, int &j, int &i) const{
        unpackIndicesUniformGrid(iGlob, nzhc, nyhc, nxhc, k, j, i);
    };

    KOKKOS_INLINE_FUNCTION int getHaloExtension(const int i, const int j, const int k) const {
        return( (hc+k)*(nx+2*hc)*(ny+2*hc) + (hc+j)*(nx+2*hc) + hc+i ); //index for the extended domain (including halo cells)
    };

    KOKKOS_INLINE_FUNCTION int getSubdomainExtension(const Parallel &par, const int i, const int j, const int k) const{
        return( k*nx_glob*ny_glob + (par.j_beg+j)*nx_glob + par.i_beg+i ); //index for the subdomain (par.j_beg+j,par.i_beg+i)
    };
    // Initialize surface domain
    // void initialise() {
    //     // physical cells onlys
    //     #if SERGHEI_MESH_UNIFORM
    //     nCell = nx*ny*nz;
    //     nCellMem = (ny+2*hc)*(nx+2*hc)*(nz+2*hc);
    //     nCellGlobal = nx_glob * ny_glob * nz_glob;
    //     #endif
    //     globalBuffer = realArr("globalBuffer", nCellGlobal);
    // };

    void fetchFieldFromGlobalBuffer(const Parallel &par, realArr &data){
        Kokkos::parallel_for("fetch_from_global_buffer", nCell , KOKKOS_CLASS_LAMBDA (int iGlob) {
            int i,j,k;
            unpackIndices(iGlob,k,j,i);
            int ii1 = getHaloExtension(i,j,k);
            int ii2 = getSubdomainExtension(par,i,j,k);
            data(ii1) = globalBuffer(ii2);
        });
    }

};


#endif
