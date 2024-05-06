/* -*- mode: c++; c-default-style: "linux" -*- */

/**
 * Subsurface model container
 **/

#ifndef _SUBSURFACE_H_
#define _SUBSURFACE_H_

#ifndef SERGHEI_SUBSURFACE_MODEL
#define SERGHEI_SUBSURFACE_MODEL 0
#endif

#if SERGHEI_SUBSURFACE_MODEL

#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include <set>

// boundary types, mimics the SWE ones, but can be completely different
#define SUB_BC_PERIODIC 1
#define SUB_BC_REFLECTIVE 2
#define SUB_BC_TRANSMISSIVE 3

real hydraulic_conductivity(real h, real alpha, real n, real Ks);

class SubsurfaceState
{

public:

  // GW variables
  realArr h;  // piezometric head
  realArr qx; // unit discharge in x-direction
  realArr qy; // unit discharge in y-direction
  realArr qz; // unit discharge in z-direction
  realArr wc; // water content

  // hydraulic conductivity
  realArr kx;
  realArr ky;
  realArr kz;

  // impermeable base
  realArr z;

  // vertical grid resolution
  realArr dz;

  //
  // reasoning
  // ---------
  // soil variables are stored in a table, soilID maps the cell to
  // van Genuchten parameters in the table soilTable
  intArr soilID;
  realArr vgTable; // van Genuchten parameters: Ks, Phi, ThetaS, ThetaR, n, m, alpha
  int nVGparam; /* number of van Genuchten parameters */
  //

  // delta fluxes
  realArr dgw0; // 4 variables (h, qx, qy, qz)
  realArr dgw1; // 4 variables (h, qx, qy, qz)
  realArr dgw2; // 4 variables (h, qx, qy, qz)

  boolArr isnodata;

  inline void
  allocate (DomainSubsurface &dom)
  {
    nVGparam = 7;

    h = realArr ("h", dom.nCellMem);
    z = realArr ("z", dom.nCellMem);
    wc = realArr ("wc", dom.nCellMem);
    kx = realArr ("kx", dom.nCellMem);
    ky = realArr ("ky", dom.nCellMem);
    kz = realArr ("kz", dom.nCellMem);

    dz = realArr ("dz", dom.nCellMem);

    soilID = intArr ("soilID", dom.nCellMem);
    isnodata = boolArr ("isnodata", dom.nCellMem);

    dgw0 = realArr ("dgw0", 4 * dom.nCellMem);
    dgw1 = realArr ("dgw1", 4 * dom.nCellMem);
    dgw2 = realArr ("dgw2", 4 * dom.nCellMem);
    vgTable = realArr ("vg", nVGparam * dom.nsoilID);
  }

};

class SubsurfaceModel {

  public:

  std::string model;        // defines the model used for subsurface flow
  std::string initialMode;  // defines the initial model
  real initialValue = 0;
  std::set<std::string> models = {"richards","linearDiffusion"};
  std::set<std::string> initialModes = {"h","theta","saturated","file-h","file-theta"};


  // calculate hydraulic conductivity on cell faces as the arithmic mean
  inline void computeFaceConductivity(SubsurfaceState &statesub, DomainSubsurface &domsub) {

      Kokkos::parallel_for( domsub.nCellDomain , KOKKOS_LAMBDA (int iGlob) {
        int i, j, k, ii, ii2, ii3;
        real alpha, n, Ks, Kcell;
        unpackIndices(iGlob,domsub.nz,domsub.ny,domsub.nx,k,j,i);
        ii = (hc+k)*(domsub.nx+2*hc)*(domsub.ny+2*hc) + (hc+j)*(domsub.nx+2*hc) + i + hc;
        ii3 = statesub.soilID(ii) * statesub.nVGparam;
        n = statesub.vgTable(ii3 + 4);
        alpha = statesub.vgTable(ii3 + 6);

        // As of 20210417, Ksat is assumed isotropic
        Ks = statesub.vgTable(ii3);
        Kcell = hydraulic_conductivity(statesub.h(ii), alpha, n, Ks);
        // x direction
        if (i == domsub.nx-1)   {statesub.kx(ii) = Kcell;}
        else {
            ii2 = ii + 1;
            statesub.kx(ii) = 0.5 * (Kcell + hydraulic_conductivity(statesub.h(ii2), alpha, n, Ks));
        }
        if (i == 0) {
            ii2 = ii - 1;
            statesub.kx(ii2) = Kcell;
        }
        // y direction
        if (j == domsub.ny-1)   {statesub.ky(ii) = Kcell;}
        else {
            ii2 = ii + domsub.nx + 2*hc;
            statesub.ky(ii) = 0.5 * (Kcell + hydraulic_conductivity(statesub.h(ii2), alpha, n, Ks));
        }
        if (j == 0) {
            ii2 = ii - domsub.nx - 2*hc;
            statesub.ky(ii2) = Kcell;
        }
        // z direction
        if (k == domsub.nz-1)   {statesub.kz(ii) = Kcell;}
        else {
            ii2 = ii + (domsub.nx+2*hc)*(domsub.ny+2*hc);
            statesub.kz(ii) = 0.5 * (Kcell + hydraulic_conductivity(statesub.h(ii2), alpha, n, Ks));
        }
        if (k == 0) {
            ii2 = ii - (domsub.nx+2*hc)*(domsub.ny+2*hc);
            statesub.kz(ii2) = Kcell;
        }
     });

  }

};

real wc_to_h(real wc, real alpha, real n, real wcs, real wcr) {
    real h, m, eps=1e-2;
    if (n < 1.0)    {std::cerr<< RERROR "van Genuchten n should be > 1!!!\n";}

    m = 1.0 - 1.0 / n;
    if (wc - wcr < eps) {wc = wcr + eps;}
    if (wc < wcs)
    {h = -(1.0/alpha) * pow(pow((wcs - wcr)/(wc - wcr),(1.0/m)) - 1.0,(1.0/n));}
    else
    {h = 0.0;}

    return h;
}

real h_to_wc(real h, real alpha, real n, real wcs, real wcr) {
    real wc, m, s;
    if (n < 1.0)    {std::cerr<< RERROR "van Genuchten n should be > 1!!!\n";}
    m = 1.0 - 1.0 / n;
    s = pow(1.0 + pow(fabs(alpha*h), n), -m);
    // calculate water content
    if (h >= 0.0)  {wc = wcs;}
    else    {wc = wcr + (wcs-wcr) * s;}
    // check to make sure wcr<=wc<=wcs
    if (wc > wcs)    {wc = wcs;}
    else if (wc < wcr)   {wc = wcr;}

    return wc;
}

real hydraulic_conductivity(real h, real alpha, real n, real Ks) {
    real m, s, Kout;
    m = 1.0 - 1.0 / n;
    s = pow(1.0 + pow(fabs(alpha*h), n), -m);
    Kout = Ks * pow(s,0.5) * pow(1-pow(1-pow(s,1.0/m),m), 2.0);
    if (h > 0.0 || Kout > Ks) {Kout = Ks;}
    return Kout;
}

real dwcdh(real h, real alpha, real n, real wcs, real wcr) {
    real deno, nume, m, Ch;
    m = 1.0 - 1.0 / n;
    nume = alpha*n*m * (wcs-wcr) * pow(fabs(alpha*h),n-1);
    deno = pow((1.0 + pow(fabs(alpha*h),n)), m+1);
    Ch = nume / deno;
    if (h > 0)  {Ch = 0.0;}
    return Ch;
}


#endif

#endif
