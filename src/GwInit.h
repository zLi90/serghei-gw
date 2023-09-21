/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _GW_INIT_H_
#define _GW_INIT_H_

#include "define.h"
#include "Domain.h"
#include "FileIO.h"
#include "GwDomain.h"
#include "GwState.h"
#include "State.h"
#include "Parallel.h"
#include "Parser.h"
#include "SWSourceSink.h"


class GwInit : public Initializer{

// Parse lines, same as Parser : ParserLine
class PsLn{
    public:
    std::string line;
    std::string key ;
    std::stringstream value;
    void lowercase()    {std::for_each(line.begin(), line.end(), [](char & c) {c = ::tolower(c);});}

    void parse()    {
        // make sure key and value are clean (in case of reuse)
        key.clear();
        value.clear();
        // If the line isn't empty and doesn't begin with a comment specifier, split it based on the colon
        if (!line.empty() && line.find("//",0) != 0) {
            uint splitloc = line.find(':',0);
            key   = line.substr(0,splitloc);
            key.erase(std::remove(key.begin(), key.end(), ' '), key.end());
            key.erase(std::remove(key.begin(), key.end(), '\t'), key.end());
            std::string val = line.substr(splitloc+1,line.length()-splitloc);
            size_t splitter = val.find("//",0);
            std::string strloc;
            if (splitter != std::string::npos)  {strloc = val.substr(0,splitter);}
            else {strloc = val;}
            value.clear();
            value.str(strloc);
        }
    }

    void print()    {std::cout << "line: " << line << "\tkey: " << key << "\tvalue: " << value.str() << std::endl;}
    };

public:

    int initialize_gw(GwState &gw, GwDomain &gdom, State &state, Domain &dom, GwBC &gbc, GwMPI &gmpi,
        Parallel &par, FileIO &io, SourceSinkData &ss, std::string inFolder, std::string outFolder) {
        int flag = -1;
        int ii, jj, kk, idx, iGlob, iGlobSW;
        real hdiff, dist;
        // Read subsurface input file
        std::string fNameIn = inFolder + "subsurface.input";
        if (!readGwFile(fNameIn, gdom, par))    {
            std::cerr << GOK << " Reading in subsurface dimensions failed." << std::endl;   return 0;
        }
        // Assumes no decomposition in the vertical direction
        gdom.nx = dom.nx;
        gdom.ny = dom.ny;
        gdom.nx_glob = dom.nx_glob;
        gdom.ny_glob = dom.ny_glob;
        gdom.nz = gdom.nz_glob;
        gdom.nxhc = gdom.nx + 2*hc;
        gdom.nyhc = gdom.ny + 2*hc;
        gdom.nzhc = gdom.nz + 2*hc;
        gdom.dx = dom.dx;
        gdom.dy = dom.dx;
        gdom.xll = dom.xll;
        gdom.yll = dom.yll;
        gdom.zll = 0.0;
        gdom.isRain = dom.isRain;
        gdom.hmin = state.hmin;
        // allocate domain
        gdom.etime = 0;
        gdom.nCellDomain = dom.nx * dom.ny * gdom.nz;
        gdom.ncells = gdom.nxhc*gdom.nyhc*(gdom.nz+2*hc);
        gdom.nhalo = 2*(gdom.nxhc*gdom.nyhc + (gdom.nxhc)*(gdom.nz+2*hc) + (gdom.nyhc)*(gdom.nz+2*hc));
        gdom.z = realArr("z", gdom.ncells);
        gdom.dz = realArr("dz", gdom.ncells);
        gdom.sinx = realArr("sinx", gdom.ncells);
        gdom.cosx = realArr("cosx", gdom.ncells);
        gdom.siny = realArr("siny", gdom.ncells);
        gdom.cosy = realArr("cosy", gdom.ncells);
        gdom.qrain = realArr("qrain", dom.ncells);
        gdom.isnodata = intArr("nodata", gdom.ncells);
        // allocate subsurface state variable
        gw.allocate(gdom);
        gmpi.allocate(gdom);
        // hpair stores the id of the halo cell, and the corresponding interior cell
        // format : hpair(idx, 0) = global id of the halo cell
        //          hpair(idx, 1) = global id of the interior cell
        //          hpair(idx, 2) = direction of the boundary: -1,1,-2,2,-3,3
        gdom.hpair = intArr2("hpair", gdom.nhalo, 3);
        idx = 0;
        for (iGlob = 0; iGlob < gdom.ncells; iGlob++) {
            unpackIndices(iGlob, gdom.nzhc, gdom.nyhc, gdom.nxhc, kk, jj, ii);
            if (ii == 0 & jj > 0 & jj < gdom.ny+1 & kk > 0 & kk < gdom.nz+1)    {
                gdom.hpair(idx,0) = iGlob;   gdom.hpair(idx,1) = iGlob+1; gdom.hpair(idx,2) = -1; idx += 1;
            }
            else if (ii == gdom.nx+hc & jj > 0 & jj < gdom.ny+1 & kk > 0 & kk < gdom.nz+1)    {
                gdom.hpair(idx,0) = iGlob;   gdom.hpair(idx,1) = iGlob-1; gdom.hpair(idx,2) = 1; idx += 1;
            }
            if (jj == 0 & ii > 0 & ii < gdom.nx+1 & kk > 0 & kk < gdom.nz+1)    {
                gdom.hpair(idx,0) = iGlob;   gdom.hpair(idx,1) = iGlob+gdom.nxhc; gdom.hpair(idx,2) = -2; idx += 1;
            }
            else if (jj == gdom.ny+hc & ii > 0 & ii < gdom.nx+1 & kk > 0 & kk < gdom.nz+1)    {
                gdom.hpair(idx,0) = iGlob;   gdom.hpair(idx,1) = iGlob-gdom.nxhc; gdom.hpair(idx,2) = 2; idx += 1;
            }
            if (kk == 0 & ii > 0 & ii < gdom.nx+1 & jj > 0 & jj < gdom.ny+1)    {
                gdom.hpair(idx,0) = iGlob;   gdom.hpair(idx,1) = iGlob+gdom.nxhc*gdom.nyhc; gdom.hpair(idx,2) = -3; idx += 1;
            }
            else if (kk == gdom.nz+hc & ii > 0 & ii < gdom.nx+1 & jj > 0 & jj < gdom.ny+1)    {
                gdom.hpair(idx,0) = iGlob;   gdom.hpair(idx,1) = iGlob-gdom.nxhc*gdom.nyhc; gdom.hpair(idx,2) = 3; idx += 1;
            }
        }
        // get z and dz for interior cells
        for (iGlob = 0; iGlob < gdom.ncells; iGlob++) {
            unpackIndices(iGlob, gdom.nzhc, gdom.nyhc, gdom.nxhc, kk, jj, ii);
            iGlobSW = packIndices(gdom.nyhc, gdom.nxhc, jj, ii);
            //if (state.z(iGlobSW) <= gdom.bottomZ) {std::cerr<< RERROR "GwDomain bottom must be lower than DEM!\n";}
            //gdom.dz(iGlob) = (state.z(iGlobSW) - gdom.bottomZ) / gdom.nz_glob;
            gdom.dz(iGlob) = gdom.thickH / gdom.nz_glob;
            gdom.z(iGlob) = state.z(iGlobSW) - (kk-hc+0.5)*gdom.dz(iGlob);
            // no data cells
            if (state.isnodata(iGlobSW) == 1)   {gdom.isnodata(iGlob) == 1;}
        }
        for (iGlob = 0; iGlob < gdom.ncells; iGlob++) {
            unpackIndices(iGlob, gdom.nzhc, gdom.nyhc, gdom.nxhc, kk, jj, ii);
            if (ii == 0)    {gdom.z(iGlob) = gdom.z(iGlob+1);}
            else if (ii == gdom.nx+1) {gdom.z(iGlob) = gdom.z(iGlob-1);}
            else if (jj == 0)   {gdom.z(iGlob) = gdom.z(iGlob+gdom.nxhc);}
            else if (jj == gdom.ny+1) {gdom.z(iGlob) = gdom.z(iGlob-gdom.nxhc);}
        }
        gmpi.mpi_sendrecv1(gdom.z, gdom, par);
        gmpi.mpi_sendrecv1(gdom.dz, gdom, par);
        // get angles for terrain-following domain
        for (iGlob = 0; iGlob < gdom.ncells; iGlob++)   {
            unpackIndices(iGlob, gdom.nzhc, gdom.nyhc, gdom.nxhc, kk, jj, ii);
            // x direction
            if (ii == 0 || ii >= gdom.nx)    {
                gdom.sinx(iGlob) = 0.0;
                gdom.cosx(iGlob) = 1.0;
            }
            else    {
                hdiff = myfabs(gdom.z(iGlob+1) - gdom.z(iGlob));
                dist = mysqrt(mypow(hdiff,2.0) + mypow(gdom.dx,2.0));
                if (gdom.z(iGlob+1) > gdom.z(iGlob))    {gdom.sinx(iGlob) = hdiff / dist;}
                else    {gdom.sinx(iGlob) = -hdiff / dist;}
                gdom.cosx(iGlob) = gdom.dx / dist;
            }
            // y direction
            if (jj == 0 || jj >= gdom.ny)    {
                gdom.siny(iGlob) = 0.0;
                gdom.cosy(iGlob) = 1.0;
            }
            else    {
                hdiff = myfabs(gdom.z(iGlob+gdom.nxhc) - gdom.z(iGlob));
                dist = mysqrt(mypow(hdiff,2.0) + mypow(gdom.dy,2.0));
                if (gdom.z(iGlob+gdom.nxhc) > gdom.z(iGlob))    {gdom.siny(iGlob) = hdiff / dist;}
                else    {gdom.siny(iGlob) = -hdiff / dist;}
                gdom.cosy(iGlob) = gdom.dy / dist;
            }
        }
        // read VG parameters
        fNameIn = inFolder + "vg.input";
        if (!readVGParameters(fNameIn, gw, gdom, par))   {
            if (par.masterproc) {
                std::cerr << RERROR << " Unable to read van Genuchten parameters." << std::endl;
                return 0;
            }
        }
        // read soil ID
        fNameIn = inFolder + "soilID.input";
        if (!readSoilID(fNameIn, gw, gdom, par)) {
            if (par.masterproc) {
                std::cerr << RERROR "Unable to read soilID from soilID.input" << std::endl;
                return 0;
            }
        }
        // read boundary conditions
        fNameIn = inFolder + "subbc.input";
        if (!readSubBCFile(fNameIn, gw, gdom, gbc, par)) {
            if (par.masterproc) {
                std::cerr << RERROR "Unable to read subsurface BC" << std::endl;
                return 0;
            }
        }
        gbc.topBC = intArr("topbc", dom.ncells);
        // read initial conditions
        fNameIn = inFolder;
        if (!setGwState(fNameIn, gw, gdom, state, gbc, par, io)) {
            if (par.masterproc) {
                std::cerr << RERROR "Unable to read subsurface IC" << std::endl;
                return 0;
            }
        }
        Kokkos::parallel_for(gdom.ncells, KOKKOS_LAMBDA(int iGlob) {
            gw.h(iGlob,0) = gw.h(iGlob,1);  gw.wc(iGlob,0) = gw.wc(iGlob,1);
        });
        gmpi.mpi_sendrecv(gw.h, gdom, par);
        gmpi.mpi_sendrecv(gw.wc, gdom, par);
        // Enforce boundary condition on lateral boundaries
        Kokkos::parallel_for( gdom.nCellDomain , KOKKOS_LAMBDA(int idom) {
            int ii, jj, kk, ii2, ii3, iGlob, ivg;
            real n, alpha, wcr, wcs;
            unpackIndices(idom, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            ii2 = kk*gdom.ny_glob + (par.j_beg+jj);
            ii3 = kk*gdom.nx_glob + (par.i_beg+ii);
            ivg = gw.soilID(iGlob) * gw.nVGparam;
            wcs = gw.vgTable(ivg+2);    wcr = gw.vgTable(ivg+3);
            n = gw.vgTable(ivg+4);  alpha = gw.vgTable(ivg+6);
            if (ii == 0 && par.px == 0 && gbc.bctypeXM != SUB_BC_NOFLOW)    {
                gw.h(iGlob-1,1) = gw.hbcX(ii2,0);
                gw.wc(iGlob-1,1) = h2wc(gw.h(iGlob-1,1), alpha, n, wcs, wcr);
            }
            else if (ii == gdom.nx-1 && par.px == par.nproc_x-1 && gbc.bctypeXP != SUB_BC_NOFLOW)    {
                gw.h(iGlob+1,1) = gw.hbcX(ii2,1);
                gw.wc(iGlob+1,1) = h2wc(gw.h(iGlob+1,1), alpha, n, wcs, wcr);
            }
            if (jj == 0 && par.py == 0 && gbc.bctypeYM != SUB_BC_NOFLOW)    {
                gw.h(iGlob-gdom.nxhc,1) = gw.hbcY(ii3,0);
                gw.wc(iGlob-gdom.nxhc,1) = h2wc(gw.h(iGlob-gdom.nxhc,1), alpha, n, wcs, wcr);
            }
            else if (jj == gdom.ny-1 && par.py == par.nproc_y-1 && gbc.bctypeYP != SUB_BC_NOFLOW)    {
                gw.h(iGlob+gdom.nxhc,1) = gw.hbcY(ii3,1);
                gw.wc(iGlob+gdom.nxhc,1) = h2wc(gw.h(iGlob+gdom.nxhc,1), alpha, n, wcs, wcr);
            }
        });

        // rainfall
        if (dom.isRain) {
            Kokkos::parallel_for( dom.ncells , KOKKOS_LAMBDA(int idom) {
                gdom.qrain(idom) = ss.rainRate(idom);
            });
        }
        io.outputIniSub(gw, gdom, par, outFolder);
        flag = 1;
        return flag;
    }


    /*
        Read input file for the subsurface
    */
    int readGwFile(std::string fNameIn, GwDomain &gdom, Parallel &par) {
        // Initialize all read-in values to -999
        gdom.nz_glob = -999;
        gdom.thickH = -999;
        par.nproc_z = -999;
        std::string strAux;
        // Read in colon-separated key: value file line by line
        std::ifstream fInStream(fNameIn);
        std::string line;
        PsLn pline;
        if (fInStream.is_open()) {
            while (std::getline(fInStream, line)) {
                pline.line = line;
                pline.parse();
                /* Rational for naming:
                ndepth  : since we have ncols and nrows
                parNz   : consistent with parNx and parNy
                thickH : conveys the intent more clearly */
                if (!pline.key.empty())  {
                    if      ( !strcmp( "ndepth" , pline.key.c_str() ) ) { pline.value >> gdom.nz_glob; }
                    else if ( !strcmp( "parNz" , pline.key.c_str() ) ) { pline.value >> par.nproc_z; }
                    else if ( !strcmp( "height"    , pline.key.c_str() ) ) { pline.value >> gdom.thickH; }
                    else if ( !strcmp( "dt_init"    , pline.key.c_str() ) ) { pline.value >> gdom.dt_init; }
                    else if ( !strcmp( "dt_max"    , pline.key.c_str() ) ) { pline.value >> gdom.dt_max; }
                    else if (!strcmp("nSoilID", pline.key.c_str()))   {pline.value >> gdom.nSoilID;}
                    else if ( !strcmp( "gw_scheme"    , pline.key.c_str() ) ) { pline.value >> gdom.gw_scheme; }
                    else if ( !strcmp( "aev"    , pline.key.c_str() ) ) { pline.value >> gdom.aev; }
                }
            }
        }
        else    {
            if (par.masterproc) {
                std::cerr<< RERROR "Unable to open " << fNameIn << "\n";    return 0;
            }
        }
        // Test to make sure all values were initialized
        if (gdom.nz_glob   == -999) { if (par.masterproc) std::cerr << RERROR "key " << "ndepth" << " not set."; exit(-1); }
        if (gdom.thickH   == -999) { if (par.masterproc) std::cerr << RERROR "key " << "height" << " not set."; exit(-1); }
        if (par.nproc_z   == -999) { if (par.masterproc) std::cerr << RERROR "key " << "parNz" << " not set."; exit(-1); }
        if (gdom.dt_init   == -999) { if (par.masterproc) std::cerr << RERROR "key " << "dt_init" << " not set."; exit(-1); }
        if (gdom.dt_max    == -999) { if (par.masterproc) std::cerr << RERROR "key " << "dt_max" << " not set."; exit(-1); }
        if (gdom.gw_scheme    == -999) { if (par.masterproc) std::cerr << RERROR "key " << "gw_scheme" << " not set."; exit(-1); }
        if (gdom.aev    == -999) { if (par.masterproc) std::cerr << RERROR "key " << "aev" << " not set."; exit(-1); }
      // Print out the values
        if (par.masterproc) {
            std::cerr << BDASH "ndepth  : "  << gdom.nz_glob    << "\n";
            std::cerr << BDASH "height : "  << gdom.thickH    << "\n";
            std::cerr << BDASH "parNz   : "  << par.nproc_z    << "\n";
        }
        if (par.masterproc)   {std::cerr<< GOK "Subsurface parameters read\n";}
        return 1;
    }

    /* io (Wed Mar 10 12:46:00 PST 2021): implemented the read in
     functionality for the van Genuchten parameters.

     vgTable structure
     =================

                        parameter index
                ------------------------------->
        s
        o |  Ks1 | Phi1 | ThS1 | ThR1 | n1 | m1 | a1
        i |  Ks2 | Phi2 | ThS2 | ThR2 | n2 | m2 | a2
        l |  Ks3 | Phi3 | ThS3 | ThR3 | n3 | m3 | a3
        i v
        d

    */
    int readVGParameters(std::string fNameIn, GwState &gw, GwDomain &gdom, Parallel &par)  {
        std::ifstream fInStream(fNameIn);
        std::string line;
        PsLn pline;
        int flag = 0;
        if (fInStream.is_open())    {
            while (std::getline(fInStream, line))   {
                pline.line = line;
        	    pline.parse();
                // std::cerr<< "key " << pline.key.c_str() << "\n";
                if (!pline.key.empty()) {
                    // VG alpha
                    if (!strcmp("alpha", pline.key.c_str())) {
            		    std::string tail;
            		    std::string head;
            		    pline.value >> tail;
            		    uint splitloc = tail.find(';', 0);
            		    for (int i = 0; i < gdom.nSoilID - 1; i ++)    {
                            head = tail.substr(0, splitloc);
                            tail = tail.substr(splitloc + 1, tail.length() - splitloc);
                            gw.vgTable(packIndices(gdom.nSoilID, gw.nVGparam, i, 6)) = std::stof(head);
                        }
            		    gw.vgTable(packIndices(gdom.nSoilID, gw.nVGparam, gdom.nSoilID - 1, 6)) = std::stof(tail);
                    }
                    // VG n
                    else if(!strcmp("n", pline.key.c_str()))    {
            		    std::string tail;
            		    std::string head;
            		    pline.value >> tail;
            		    uint splitloc = tail.find(';', 0);
            		    for (int i = 0; i < gdom.nSoilID - 1; i ++)    {
                            head = tail.substr(0, splitloc);
                            tail = tail.substr(splitloc + 1, tail.length() - splitloc);
                            gw.vgTable(packIndices(gdom.nSoilID, gw.nVGparam, i, 4)) = std::stof(head);
                        }
            		    gw.vgTable(packIndices(gdom.nSoilID, gw.nVGparam, gdom.nSoilID - 1, 4)) = std::stof(tail);
                    }
                    // VG Ks
                    else if(!strcmp("Ks", pline.key.c_str()))   {
            		    std::string tail;
            		    std::string head;
            		    pline.value >> tail;
            		    uint splitloc = tail.find(';', 0);
            		    for (int i = 0; i < gdom.nSoilID - 1; i ++)    {
                            head = tail.substr(0, splitloc);
                            tail = tail.substr(splitloc + 1, tail.length() - splitloc);
                            gw.vgTable(packIndices(gdom.nSoilID, gw.nVGparam, i, 0)) = std::stof(head);
                        }
            		    gw.vgTable(packIndices(gdom.nSoilID, gw.nVGparam, gdom.nSoilID - 1, 0)) = std::stof(tail);
                    }
                    // VG porosity
                    else if(!strcmp("Phi", pline.key.c_str()))  {
            		    std::string tail;
            		    std::string head;
            		    pline.value >> tail;
            		    uint splitloc = tail.find(';', 0);
            		    for (int i = 0; i < gdom.nSoilID - 1; i ++)    {
                            head = tail.substr(0, splitloc);
                            tail = tail.substr(splitloc + 1, tail.length() - splitloc);
                            gw.vgTable(packIndices(gdom.nSoilID, gw.nVGparam, i, 1)) = std::stof(head);
                        }
                        gw.vgTable(packIndices(gdom.nSoilID, gw.nVGparam, gdom.nSoilID - 1, 1)) = std::stof(tail);
                    }
                    // VG wcs
                    else if(!strcmp("ThetaR", pline.key.c_str()))   {
            		    std::string tail;
            		    std::string head;
            		    pline.value >> tail;
            		    uint splitloc = tail.find(';', 0);
            		    for (int i = 0; i < gdom.nSoilID - 1; i ++)    {
                            head = tail.substr(0, splitloc);
                            tail = tail.substr(splitloc + 1, tail.length() - splitloc);
                            gw.vgTable(packIndices(gdom.nSoilID, gw.nVGparam, i, 3)) = std::stof(head);
                        }
            		    gw.vgTable(packIndices(gdom.nSoilID, gw.nVGparam, gdom.nSoilID - 1, 3)) = std::stof(tail);
                    }
                    // VG wcr
            		else if(!strcmp("ThetaS", pline.key.c_str())) {
            		    std::string tail;
            		    std::string head;
            		    pline.value >> tail;
            		    uint splitloc = tail.find(';', 0);
            		    for (int i = 0; i < gdom.nSoilID - 1; i ++)    {
                            head = tail.substr(0, splitloc);
                            tail = tail.substr(splitloc + 1, tail.length() - splitloc);
                            gw.vgTable(packIndices(gdom.nSoilID, gw.nVGparam, i, 2)) = std::stof(head);
                        }
                        gw.vgTable(packIndices(gdom.nSoilID, gw.nVGparam, gdom.nSoilID - 1, 2)) = std::stof(tail);
                    }
            		else  {
                        if (par.masterproc) {
                            std::cerr << RERROR "key " << pline.key << " not understood in file " << fNameIn << std::endl;
                        }
                        flag = -1;
                    }
                }
            }
        }
        if (par.masterproc) {std::cerr<< GOK "van Genuchten parameters read\n";}
        flag = 1;
        return flag;
    }

    /*
        Read soil ID
    */
    int readSoilID(std::string fNameIn, GwState &gw, GwDomain &gdom, Parallel &par)    {
        std::ifstream fInStream(fNameIn);
        std::string line, str;
        int ii, jj, kk, ii2, idx, iGlob, tmp, n_soil, ndata=gdom.ny_glob*gdom.nx_glob*gdom.nz_glob;
        intArr tmpVar=intArr("var", ndata);
        // read soil ID file
        if (fInStream.is_open())    {
            std::getline(fInStream,str,' ');
            std::getline(fInStream,str);
            std::stringstream(str) >> n_soil;
            //compare the values t* with the DEM file just to check if we are using the same values, otherwise error
            if (n_soil != gdom.nSoilID)    {
                if (par.masterproc)  {std::cerr<< RERROR "Number of soilID != that specified in the VG Table. Unable to continue\n";}
                if (par.masterproc)  {std::cerr << BDASH "n_soil: " 	<< n_soil 	<< ", in the Table : " << gw.nVGparam <<"\n";}
            	return 0;
            }
            if (n_soil == 1)    {for (int ii=0; ii<ndata; ii++) {tmpVar(ii) = 0;}}
            else {
                for (int ii=0; ii<ndata; ii++) {
          			if (!fInStream.fail() && !fInStream.eof()) {
          				fInStream >> tmp;
          				tmpVar(ii)=tmp;
          			}
                    else {if(par.masterproc){std::cerr<< RERROR "Error reading soil ID file. Not enough data\n";  return 0;}}
          		}
            }
      		fInStream.close();
      	}
        else {
            if (gdom.nSoilID == 1)  {
                for (int ii=0; ii<ndata; ii++) {tmpVar(ii) = 0;}
            }
            else {
                if(par.masterproc){std::cerr<< RERROR "Unable to open soil ID file!\n";  return 0;}
            }
        }
        // assign soil ID to cells
        for (idx = 0; idx < gdom.nCellDomain; idx++)  {
            unpackIndices(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
            // global index for this rank
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            // global index for the entire domain
            ii2 = kk*gdom.nx_glob*gdom.ny_glob + (par.j_beg+jj)*(gdom.nx_glob)+par.i_beg+ii;
            gw.soilID(iGlob) = tmpVar(ii2);
        }
        // halo cells
        for (idx = 0; idx < gdom.nhalo; idx++)  {
            gw.soilID(gdom.hpair(idx,0)) = gw.soilID(gdom.hpair(idx,1));
        }
        if (par.masterproc)   {std::cerr<<GOK "Soil ID  set\n";}
        return 1;
    }

    /*
        Read boundary condition for the subsurface
    */
    int readSubBCFile(std::string fNameIn, GwState &gw, GwDomain &gdom, GwBC &gbc, Parallel &par) {
        std::ifstream fInStream(fNameIn);
        std::string line;
        PsLn pline;
        if (fInStream.is_open()) {
            while (std::getline(fInStream, line)) {
                pline.line = line;
                pline.parse();
                if(!pline.key.empty()){
                    if (!strcmp("bctype-xp", pline.key.c_str())){pline.value >> gbc.bctypeXP;}
                    else if (!strcmp("bctype-xm", pline.key.c_str())){pline.value >> gbc.bctypeXM;}
                    else if (!strcmp("bctype-yp", pline.key.c_str())){pline.value >> gbc.bctypeYP;}
                    else if (!strcmp("bctype-ym", pline.key.c_str())){pline.value >> gbc.bctypeYM;}
                    else if (!strcmp("bctype-zp", pline.key.c_str())){pline.value >> gbc.bctypeZP;}
                    else if (!strcmp("bctype-zm", pline.key.c_str())){pline.value >> gbc.bctypeZM;}
                    else if (!strcmp("hbc-xp", pline.key.c_str())){pline.value >> gbc.hbcXP;}
                    else if (!strcmp("hbc-xm", pline.key.c_str())){pline.value >> gbc.hbcXM;}
                    else if (!strcmp("hbc-yp", pline.key.c_str())){pline.value >> gbc.hbcYP;}
                    else if (!strcmp("hbc-ym", pline.key.c_str())){pline.value >> gbc.hbcYM;}
                    else if (!strcmp("hbc-zp", pline.key.c_str())){pline.value >> gbc.hbcZP;}
                    else if (!strcmp("hbc-zm", pline.key.c_str())){pline.value >> gbc.hbcZM;}
                    else if (!strcmp("qbc-xp", pline.key.c_str())){pline.value >> gbc.qbcXP;}
                    else if (!strcmp("qbc-xm", pline.key.c_str())){pline.value >> gbc.qbcXM;}
                    else if (!strcmp("qbc-yp", pline.key.c_str())){pline.value >> gbc.qbcYP;}
                    else if (!strcmp("qbc-ym", pline.key.c_str())){pline.value >> gbc.qbcYM;}
                    else if (!strcmp("qbc-zp", pline.key.c_str())){pline.value >> gbc.qbcZP;}
                    else if (!strcmp("qbc-zm", pline.key.c_str())){pline.value >> gbc.qbcZM;}
                    else {
                        if(par.masterproc)
                        {std::cerr << RERROR << "In subbc.input: Key " << pline.key << " not understood." << std::endl;   return 0;}
                    }
                }
            }
            fInStream.close();
        }
        else {if (par.masterproc)    {std::cerr << YEXC << "subbc.input not found. Default boundaries used." << std::endl; return 0;}}
        if (par.masterproc) {std::cerr<< GOK "Subsurface BC set\n";}
        return 1;
    }

    // read subsurface initial conditions
    int setGwState(std::string inFolder, GwState &gw, GwDomain &gdom, State &state, GwBC &gbc, Parallel &par, FileIO &io){
        int ii, jj, kk, idx, ivg, iGlob, iGlobSW;
        real wcs, wcr, n, alpha;
        std::ifstream fInStream(inFolder + "subsurface.input");
        std::string line;
        PsLn pline;
        //SubsurfaceModel sub;
        std::string tempStr;
        gw.initialMode="saturated";
        // read initial mode and value
        if (fInStream.is_open()){
            while (std::getline(fInStream, line)) {
                pline.line = line;
                pline.lowercase();
                pline.parse();
                // If the line was valid and a key is stored
                if(!pline.key.empty()){
                    // Match the key, and store the value
                    if(!strcmp("initialmode",pline.key.c_str())){ pline.value >> gw.initialMode;}
                    if(!strcmp("initialvalue",pline.key.c_str())){ pline.value >> gw.initialValue ; }
                }
            }
        }
        else {
            if (par.masterproc){
                std::cerr << RERROR "File " << inFolder + "subsurface.input" << " not found" << std::endl; return 0;
            }
        }
        // initialize the primary variables
        for (iGlob = 0; iGlob < gdom.ncells; iGlob++)   {
            gw.wc(iGlob,0) = 0.0;   gw.wc(iGlob,1) = 0.0;
            gw.h(iGlob,0) = 0.0;    gw.h(iGlob,1) = 0.0;    gw.wc(iGlob,2) = 0.0;
        }
        // read initial condition from file
        if(!gw.initialMode.compare("file-h")) {
            tempStr = "head.input";
            readGwICFile(tempStr, inFolder, gw, gdom, par);
        }
        else if (!gw.initialMode.compare("file-theta")){
            tempStr = "theta.input";
            readGwICFile(tempStr, inFolder, gw, gdom, par);
        }
        else if (!gw.initialMode.compare("file-wt"))    {
            tempStr = "wt.input";
            readGwICFile(tempStr, inFolder, gw, gdom, par);
        }
        else if (!gw.initialMode.compare("saturated")){
            for (iGlob = 0; iGlob < gdom.ncells; iGlob++)   {
                unpackIndices(iGlob, gdom.nzhc, gdom.nyhc, gdom.nxhc, kk, jj, ii);
                iGlobSW = packIndices(gdom.nyhc, gdom.nxhc, jj, ii);
                idx = gw.soilID(iGlob) * gw.nVGparam;
                wcs = gw.vgTable(idx + 2);
                gw.wc(iGlob,0) = wcs;   gw.wc(iGlob,1) = wcs;
                gw.h(iGlob,1) = state.h(iGlobSW) + (kk-0.5)*gdom.dz(iGlob);
                gw.h(iGlob,0) = gw.h(iGlob,1);
            }
        }
        else {
            if (par.masterproc) {
                std::cerr << RERROR "Initial mode must be file-h, file-theta or saturated!" << std::endl; return 0;
            }
        }
        if (par.masterproc){std::cerr<<GOK "Subsurface initial condition set" << std::endl;}

        // set up boundary conditions
        for (idx = 0; idx < gdom.nhalo; idx++)  {
            // get soil parameters
            ivg = gw.soilID(gdom.hpair(idx,1)) * gw.nVGparam;
            wcs = gw.vgTable(ivg + 2);  wcr = gw.vgTable(ivg + 3);
            n = gw.vgTable(ivg + 4);    alpha = gw.vgTable(ivg + 6);
            // apply boundary conditions
            if (gdom.hpair(idx,2) == -1)    {
                gw.h(gdom.hpair(idx,0),1) = gbc.hbcXM;
            }
            else if (gdom.hpair(idx,2) == 1)    {
                gw.h(gdom.hpair(idx,0),1) = gbc.hbcXP;
            }
            else if (gdom.hpair(idx,2) == -2)    {
                gw.h(gdom.hpair(idx,0),1) = gbc.hbcYM;
            }
            else if (gdom.hpair(idx,2) == 2)    {
                gw.h(gdom.hpair(idx,0),1) = gbc.hbcYP;
            }
            else if (gdom.hpair(idx,2) == -3)    {
                switch (gbc.bctypeZM) {
                    case SUB_BC_H_SWE:
                        //unpackIndices(icell, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
                        //i2d = (hc+jj)*gdom.nxhc + ii + hc;
                        //h_h(gdom.hpair(idx,0),1) = state.h(i2d);
                        gw.h(gdom.hpair(idx,0),1) = gbc.hbcZM;
                    case SUB_BC_Q_CONST:
                        gw.h(gdom.hpair(idx,0),1) = gw.h(gdom.hpair(idx,1),1);
                    default:
                        gw.h(gdom.hpair(idx,0),1) = gbc.hbcZM;
                }
            }
            else if (gdom.hpair(idx,2) == 3)    {
                gw.h(gdom.hpair(idx,0),1) = gbc.hbcZP;
            }
            gw.wc(gdom.hpair(idx,0),1) = h2wc(gw.h(gdom.hpair(idx,0),1), alpha, n, wcs, wcr);
        }
        // Setting spatially-distributed boundary conditions
        if (gbc.bctypeXM == 4)  {
            tempStr = "hbcxm.input";
            readGwBCFileX(tempStr, inFolder, gw, gdom, par, -1);
        }
        if (gbc.bctypeXP == 4)  {
            tempStr = "hbcxp.input";
            readGwBCFileX(tempStr, inFolder, gw, gdom, par, 1);
        }
        if (gbc.bctypeYM == 4)  {
            tempStr = "hbcym.input";
            readGwBCFileY(tempStr, inFolder, gw, gdom, par, -1);
        }
        if (gbc.bctypeYP == 4)  {
            tempStr = "hbcyp.input";
            readGwBCFileY(tempStr, inFolder, gw, gdom, par, 1);
        }
        if (gbc.bctypeZM == 4)  {
            tempStr = "hbczm.input";
            readGwBCFileZ(tempStr, inFolder, gw, gdom, par);
        }
        // gw.h.modify<dualDbl::host_mirror_space> ();
        // gw.wc.modify<dualDbl::host_mirror_space> ();
        if (par.masterproc){std::cerr<<GOK "Subsurface boundary condition set" << std::endl;}
        return 1;
    }

    /*
        Read initial condition from file
    */
    int readGwICFile(std::string fNameIn, std::string fDirIn, GwState &gw, GwDomain &gdom, Parallel &par)   {
        std::string fname = fDirIn + fNameIn;
        std::ifstream fInStream(fname);
        std::string line;
        int tnx, tny, iGlob, iGlobSW, idx, ivg, ii, jj, kk, ii2;
        real tmp, wcs, wcr, n, alpha, nodata_value ;
        int ndata;
        if (!strcmp(fNameIn.c_str(), "wt.input")) {ndata = gdom.nx_glob*gdom.ny_glob;}
        else {ndata = gdom.nx_glob*gdom.ny_glob*gdom.nz_glob;}
        realArr tmpVar = realArr("var", ndata);
        std::string str;
     	if (fInStream.is_open()) {
            std::getline(fInStream,str,' ');
            std::getline(fInStream,str);
            std::stringstream(str) >> tnx;
            std::getline(fInStream,str,' ');
            std::getline(fInStream,str);
            std::stringstream(str) >> tny;
            std::getline(fInStream,str,' ');
            std::getline(fInStream,str);
            std::stringstream(str) >> nodata_value;
     		//compare the values t* with the DEM file just to check if we are using the same values, otherwise error
     		if (gdom.ny_glob !=tny || gdom.nx_glob !=tnx) {
                if (par.masterproc) {
                    std::cerr<< RERROR "IC file parameters don't match DEM parameters. Unable to continue\n";
                    if (par.masterproc) {
                        std::cerr << BDASH "nx_glob: " 	<< gdom.nx_glob 	<< tnx <<"\n";
                        std::cerr << BDASH "ny_glob: "<< gdom.ny_glob 	<< tny << "\n";
                    }
                    return 0;
                }
            }
            // read and store data into a temporary view
            for (int ii = 0; ii < ndata; ii++)  {
                if (!fInStream.fail() && !fInStream.eof()) {
                    fInStream >> tmp;
                    tmpVar(ii)=tmp;
                }
                else {
                    if (par.masterproc) {std::cerr<< RERROR "Error reading IC file. Not enough data\n";  return 0;}
                }
            }
            fInStream.close();
        }
        // Copy data into head or water content
        if (!strcmp(fNameIn.c_str(), "head.input")) {
            for (idx = 0; idx < gdom.nCellDomain; idx++)    {
                unpackIndices(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
                // get global index
                iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                // get soil parameters
                ivg = gw.soilID(iGlob) * gw.nVGparam;
                wcs = gw.vgTable(ivg+2);
                wcr = gw.vgTable(ivg+3);
                n = gw.vgTable(ivg+4);
                alpha = gw.vgTable(ivg+6);
                // get head and water content
                ii2 = kk*gdom.nx_glob*gdom.ny_glob + (par.j_beg+jj)*(gdom.nx_glob)+par.i_beg+ii;
                gw.h(iGlob,1) = tmpVar(ii2);
                gw.wc(iGlob,1) = h2wc(gw.h(iGlob,1), alpha, n, wcs, wcr);
                gw.h(iGlob,0) = gw.h(iGlob,1);
                gw.wc(iGlob,0) = gw.wc(iGlob,1);
            }
        }
        else if (!strcmp(fNameIn.c_str(), "theta.input")) {
            for (idx = 0; idx < gdom.nCellDomain; idx++)    {
                unpackIndices(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
                // get global index
                iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                // get soil parameters
                ivg = gw.soilID(iGlob) * gw.nVGparam;
                wcs = gw.vgTable(ivg+2);
                wcr = gw.vgTable(ivg+3);
                n = gw.vgTable(ivg+4);
                alpha = gw.vgTable(ivg+6);
                // get head and water content
                ii2 = kk*gdom.nx_glob*gdom.ny_glob + (par.j_beg+jj)*(gdom.nx_glob)+par.i_beg+ii;
                gw.wc(iGlob,1) = tmpVar(ii2);
                gw.h(iGlob,1) = wc2h(gw.wc(iGlob,1), alpha, n, wcs, wcr);
                gw.h(iGlob,0) = gw.h(iGlob,1);
                gw.wc(iGlob,0) = gw.wc(iGlob,1);
            }
        }
        else if (!strcmp(fNameIn.c_str(), "wt.input")) {
            for (idx = 0; idx < gdom.nCellDomain; idx++)    {
                unpackIndices(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
                // get global index
                iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                // get soil parameters
                ivg = gw.soilID(iGlob) * gw.nVGparam;
                wcs = gw.vgTable(ivg+2);
                wcr = gw.vgTable(ivg+3);
                n = gw.vgTable(ivg+4);
                alpha = gw.vgTable(ivg+6);
                // get head and water content
                ii2 = (par.j_beg+jj)*(gdom.nx_glob)+par.i_beg+ii;
                if (tmpVar(ii2) == nodata_value)    {
                    gw.h(iGlob,1) = 0.0;    gw.wc(iGlob,1) = wcs;
                }
                else {
                    gw.h(iGlob,1) = tmpVar(ii2) - gdom.z(iGlob);
                    gw.wc(iGlob,1) = h2wc(gw.h(iGlob,1), alpha, n, wcs, wcr);
                }
                gw.h(iGlob,0) = gw.h(iGlob,1);
                gw.wc(iGlob,0) = gw.wc(iGlob,1);
            }
        }
        else    {
    	 	if (par.masterproc) {std::cerr<< RERROR "Error reading head/theta IC file. File name might be wrong.\n"; return 0;}
        }
        if (par.masterproc)   {std::cerr<<GOK "Subsurface head/water content set\n";}
        return 1;
    }


    /*
        Read spatially-distributed boundary conditions
        (As of 2023-06-16, Only support ZM, XM, XP faces)
    */
    int readGwBCFileZ(std::string fNameIn, std::string fDirIn, GwState &gw, GwDomain &gdom, Parallel &par)   {
        std::string fname = fDirIn + fNameIn;
        std::ifstream fInStream(fname);
        std::string line;
        int tnx, tny, iGlob, idx, ivg, ii, jj, kk, ii2;
        real tmp, wcs, wcr, n, alpha ;
        int ndata = gdom.nx_glob*gdom.ny_glob;
        realArr tmpVar = realArr("var", ndata);
        std::string str;
     	if (fInStream.is_open()) {
            std::getline(fInStream,str,' ');
            std::getline(fInStream,str);
            std::stringstream(str) >> tnx;
            std::getline(fInStream,str,' ');
            std::getline(fInStream,str);
            std::stringstream(str) >> tny;
     		//compare the values t* with the DEM file just to check if we are using the same values, otherwise error
     		if (gdom.ny_glob !=tny || gdom.nx_glob !=tnx) {
                if (par.masterproc) {
                    std::cerr<< RERROR "BC (hbczm) file parameters don't match DEM parameters. Unable to continue\n";
                    if (par.masterproc) {
                        std::cerr << BDASH "nx_glob: " 	<< gdom.nx_glob 	<< tnx <<"\n";
                        std::cerr << BDASH "ny_glob: "<< gdom.ny_glob 	<< tny << "\n";
                    }
                    return 0;
                }
            }
            // read and store data into a temporary view
            for (int ii = 0; ii < ndata; ii++)  {
                if (!fInStream.fail() && !fInStream.eof()) {
                    fInStream >> tmp;
                    tmpVar(ii)=tmp;
                }
                else {
                    if (par.masterproc) {std::cerr<< RERROR "Error reading BC (hbczm) file. Not enough data\n";  return 0;}
                }
            }
            fInStream.close();
        }
        // Copy data into head or water content
        if (!strcmp(fNameIn.c_str(), "hbczm.input")) {
            for (idx = 0; idx < gdom.nCellDomain; idx++)    {
                unpackIndices(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
                if (kk == 0)    {
                    // get global index
                    iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                    // get soil parameters
                    ivg = gw.soilID(iGlob) * gw.nVGparam;
                    wcs = gw.vgTable(ivg+2);    wcr = gw.vgTable(ivg+3);
                    n = gw.vgTable(ivg+4);  alpha = gw.vgTable(ivg+6);
                    // get index of the halo cell
                    iGlob -= gdom.nxhc*gdom.nyhc;
                    // get head and water content
                    ii2 = (par.j_beg+jj)*(gdom.nx_glob)+par.i_beg+ii;
                    gw.h(iGlob,1) = tmpVar(ii2);
                    gw.wc(iGlob,1) = h2wc(gw.h(iGlob,1), alpha, n, wcs, wcr);
                    gw.h(iGlob,0) = gw.h(iGlob,1);
                    gw.wc(iGlob,0) = gw.wc(iGlob,1);
                }
            }
        }
        else    {
    	 	if (par.masterproc) {std::cerr<< RERROR "Error reading head BC (hbczm) file. File name might be wrong.\n"; return 0;}
        }
        if (par.masterproc)   {std::cerr<<GOK "Subsurface head BC set\n";}
        return 1;
    }


    int readGwBCFileX(std::string fNameIn, std::string fDirIn, GwState &gw, GwDomain &gdom, Parallel &par, int dir)   {
        std::string fname = fDirIn + fNameIn;
        std::ifstream fInStream(fname);
        std::string line;
        int tnz, tny, iGlob, idx, ivg, ii, jj, kk, ii2;
        real tmp, wcs, wcr, n, alpha ;
        int ndata = gdom.nz_glob*gdom.ny_glob;
        realArr tmpVar = realArr("var", ndata);
        std::string str;
     	if (fInStream.is_open()) {
            std::getline(fInStream,str,' ');
            std::getline(fInStream,str);
            std::stringstream(str) >> tny;
            std::getline(fInStream,str,' ');
            std::getline(fInStream,str);
            std::stringstream(str) >> tnz;
     		//compare the values t* with the DEM file just to check if we are using the same values, otherwise error
     		if (gdom.ny_glob !=tny || gdom.nz_glob !=tnz) {
                if (par.masterproc) {
                    std::cerr<< RERROR "BC (hbcx) file parameters don't match DEM parameters. Unable to continue\n";
                    if (par.masterproc) {
                        std::cerr << BDASH "ny_glob: " 	<< gdom.ny_glob 	<< tny <<"\n";
                        std::cerr << BDASH "nz_glob: "<< gdom.nz_glob 	<< tnz << "\n";
                    }
                    return 0;
                }
            }
            // read and store data into a temporary view
            for (int ii = 0; ii < ndata; ii++)  {
                if (!fInStream.fail() && !fInStream.eof()) {
                    fInStream >> tmp;
                    tmpVar(ii)=tmp;
                }
                else {
                    if (par.masterproc) {std::cerr<< RERROR "Error reading BC (hbcx) file. Not enough data\n";  return 0;}
                }
            }
            fInStream.close();
        }
        // Copy data into head or water content
        if (!strcmp(fNameIn.c_str(), "hbcxp.input") && dir == 1) {
            for (idx = 0; idx < gdom.nCellDomain; idx++)    {
                unpackIndices(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
                if (ii == gdom.nx-1)    {
                    // get global index
                    iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                    // get soil parameters
                    ivg = gw.soilID(iGlob) * gw.nVGparam;
                    wcs = gw.vgTable(ivg+2);    wcr = gw.vgTable(ivg+3);
                    n = gw.vgTable(ivg+4);  alpha = gw.vgTable(ivg+6);
                    // get index of the halo cell
                    iGlob += 1;
                    // get head and water content
                    ii2 = kk*gdom.ny_glob + (par.j_beg+jj);
                    gw.h(iGlob,1) = tmpVar(ii2);
                    gw.wc(iGlob,1) = h2wc(gw.h(iGlob,1), alpha, n, wcs, wcr);
                    gw.h(iGlob,0) = gw.h(iGlob,1);
                    gw.wc(iGlob,0) = gw.wc(iGlob,1);

                    gw.hbcX(ii2,1) = tmpVar(ii2);
                }
            }
        }
        else if (!strcmp(fNameIn.c_str(), "hbcxm.input") && dir == -1)  {
            for (idx = 0; idx < gdom.nCellDomain; idx++)    {
                unpackIndices(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
                if (ii == 0)    {
                    // get global index
                    iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                    // get soil parameters
                    ivg = gw.soilID(iGlob) * gw.nVGparam;
                    wcs = gw.vgTable(ivg+2);    wcr = gw.vgTable(ivg+3);
                    n = gw.vgTable(ivg+4);  alpha = gw.vgTable(ivg+6);
                    // get index of the halo cell
                    iGlob -= 1;
                    // get head and water content
                    ii2 = kk*gdom.ny_glob + (par.j_beg+jj);
                    gw.h(iGlob,1) = tmpVar(ii2);
                    gw.wc(iGlob,1) = h2wc(gw.h(iGlob,1), alpha, n, wcs, wcr);
                    gw.h(iGlob,0) = gw.h(iGlob,1);
                    gw.wc(iGlob,0) = gw.wc(iGlob,1);

                    gw.hbcX(ii2,0) = tmpVar(ii2);
                }
            }
        }
        else    {
    	 	if (par.masterproc) {std::cerr<< RERROR "Error reading head BC (hbcx) file. File name might be wrong.\n"; return 0;}
        }
        if (par.masterproc)   {std::cerr<<GOK "Subsurface head BC set\n";}
        return 1;
    }

    int readGwBCFileY(std::string fNameIn, std::string fDirIn, GwState &gw, GwDomain &gdom, Parallel &par, int dir)   {
        std::string fname = fDirIn + fNameIn;
        std::ifstream fInStream(fname);
        std::string line;
        int tnz, tnx, iGlob, idx, ivg, ii, jj, kk, ii2;
        real tmp, wcs, wcr, n, alpha ;
        int ndata = gdom.nz_glob*gdom.nx_glob;
        realArr tmpVar = realArr("var", ndata);
        std::string str;
     	if (fInStream.is_open()) {
            std::getline(fInStream,str,' ');
            std::getline(fInStream,str);
            std::stringstream(str) >> tnx;
            std::getline(fInStream,str,' ');
            std::getline(fInStream,str);
            std::stringstream(str) >> tnz;
     		//compare the values t* with the DEM file just to check if we are using the same values, otherwise error
     		if (gdom.nx_glob !=tnx || gdom.nz_glob !=tnz) {
                if (par.masterproc) {
                    std::cerr<< RERROR "BC (hbcy) file parameters don't match DEM parameters. Unable to continue\n";
                    if (par.masterproc) {
                        std::cerr << BDASH "nx_glob: " 	<< gdom.nx_glob 	<< tnx <<"\n";
                        std::cerr << BDASH "nz_glob: "<< gdom.nz_glob 	<< tnz << "\n";
                    }
                    return 0;
                }
            }
            // read and store data into a temporary view
            for (int ii = 0; ii < ndata; ii++)  {
                if (!fInStream.fail() && !fInStream.eof()) {
                    fInStream >> tmp;
                    tmpVar(ii)=tmp;
                }
                else {
                    if (par.masterproc) {std::cerr<< RERROR "Error reading BC (hbcy) file. Not enough data\n";  return 0;}
                }
            }
            fInStream.close();
        }
        // Copy data into head or water content
        if (!strcmp(fNameIn.c_str(), "hbcyp.input") && dir == 1) {
            for (idx = 0; idx < gdom.nCellDomain; idx++)    {
                unpackIndices(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
                if (jj == gdom.ny-1)    {
                    // get global index
                    iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                    // get soil parameters
                    ivg = gw.soilID(iGlob) * gw.nVGparam;
                    wcs = gw.vgTable(ivg+2);    wcr = gw.vgTable(ivg+3);
                    n = gw.vgTable(ivg+4);  alpha = gw.vgTable(ivg+6);
                    // get index of the halo cell
                    iGlob += gdom.nxhc;
                    // get head and water content
                    ii2 = kk*gdom.nx_glob + (par.i_beg+ii);
                    gw.h(iGlob,1) = tmpVar(ii2);
                    gw.wc(iGlob,1) = h2wc(gw.h(iGlob,1), alpha, n, wcs, wcr);
                    gw.h(iGlob,0) = gw.h(iGlob,1);
                    gw.wc(iGlob,0) = gw.wc(iGlob,1);

                    gw.hbcY(ii2,1) = tmpVar(ii2);
                }
            }
        }
        else if (!strcmp(fNameIn.c_str(), "hbcym.input") && dir == -1)  {
            for (idx = 0; idx < gdom.nCellDomain; idx++)    {
                unpackIndices(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
                if (jj == 0)    {
                    // get global index
                    iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
                    // get soil parameters
                    ivg = gw.soilID(iGlob) * gw.nVGparam;
                    wcs = gw.vgTable(ivg+2);    wcr = gw.vgTable(ivg+3);
                    n = gw.vgTable(ivg+4);  alpha = gw.vgTable(ivg+6);
                    // get index of the halo cell
                    iGlob -= gdom.nxhc;
                    // get head and water content
                    ii2 = kk*gdom.nx_glob + (par.i_beg+ii);
                    gw.h(iGlob,1) = tmpVar(ii2);
                    gw.wc(iGlob,1) = h2wc(gw.h(iGlob,1), alpha, n, wcs, wcr);
                    gw.h(iGlob,0) = gw.h(iGlob,1);
                    gw.wc(iGlob,0) = gw.wc(iGlob,1);

                    gw.hbcY(ii2,0) = tmpVar(ii2);
                }
            }
        }
        else    {
    	 	if (par.masterproc) {std::cerr<< RERROR "Error reading head BC (hbcy) file. File name might be wrong.\n"; return 0;}
        }
        if (par.masterproc)   {std::cerr<<GOK "Subsurface head BC set\n";}
        return 1;
    }


};

#endif
