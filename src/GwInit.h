/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _GW_INIT_H_
#define _GW_INIT_H_

#include "define.h"
#include "Domain.h"
#include "FileIO.h"
#include "GwBC.h"
#include "GwDomain.h"
#include "GwState.h"
#include "State.h"
#include "Parallel.h"
#include "Parser.h"
#include "SourceSink.h"


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

    int initialize_gw(GwState &gw, GwDomain &gdom, State &state, Domain &dom, SubsurfaceBoundaries &gbc, GwMPI &gmpi, GwIntegrator &gint,
        Parallel &par, FileIO &io, SourceSink &ss, std::string inFolder, std::string outFolder) {
        int flag = -1;
        int ii, jj, kk, idx, iGlob, iGlobSW;
        real hdiff, dist, dz_base;
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
        gdom.dx = dom.dxConst;
        gdom.dy = dom.dxConst;
        gdom.xll = dom.xll;
        gdom.yll = dom.yll;
        gdom.zll = 0.0;
        gdom.hmin = state.hmin;
        // allocate domain
        gdom.etime = 0.0;
        gdom.nCell = dom.nx * dom.ny * gdom.nz;
        gdom.nCellMem = gdom.nxhc*gdom.nyhc*gdom.nzhc;
        gdom.nhalo = 2*(gdom.nxhc*gdom.nyhc + gdom.nxhc*gdom.nzhc + gdom.nyhc*gdom.nzhc);
        gdom.nCellSw = dom.nCell;
        gdom.nCellSwMem = dom.nCellMem;
        gdom.z = realArr("z", gdom.nCellMem);
        gdom.dz = realArr("dz", gdom.nCellMem);
        gdom.sinx = realArr("sinx", gdom.nCellMem);
        gdom.cosx = realArr("cosx", gdom.nCellMem);
        gdom.siny = realArr("siny", gdom.nCellMem);
        gdom.cosy = realArr("cosy", gdom.nCellMem);
        gdom.rainRate = realArr("rain", dom.nCellMem);
        gdom.evapRate = realArr("evap", dom.nCellMem);
        gdom.isnodata = intArr("nodata", gdom.nCellMem);
        // allocate subsurface state variable
        gw.allocate(gdom);
        gmpi.allocate(gdom);
        // get z and dz for interior cells
        for (iGlob = 0; iGlob < gdom.nCellMem; iGlob++) {
            gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
            iGlobSW = packIndicesUniformGrid(gdom.nyhc, gdom.nxhc, jj, ii);
            dz_base = gdom.thickH / gdom.nz_glob;
            // Note that when dz_multiplier > 1, the actual domain height will be > gdom.thickH
            gdom.dz(iGlob) = dz_base * mypow(gdom.dz_multiplier, kk);
            gdom.z(iGlob) = state.z(iGlobSW) - (kk-hc+0.5)*gdom.dz(iGlob);
            // no data cells
            if (state.isnodata(iGlobSW) == 1)   {gdom.isnodata(iGlob) == 1;}
        }
        gmpi.mpi_sendrecv1(gdom.z, gdom, par);
        gmpi.mpi_sendrecv1(gdom.dz, gdom, par);
        // get angles for terrain-following domain
        for (iGlob = 0; iGlob < gdom.nCellMem; iGlob++)   {
            gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
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
        // read rainfall
        gdom.isRain = dom.isRain;
        gdom.isEvap = dom.isEvap;
        #if SERGHEI_SWE_MODEL
        if (gdom.isRain) {
            Kokkos::deep_copy(gdom.rainRate, ss.swss.rainRate);
        }
        if (gdom.isEvap) {
            Kokkos::deep_copy(gdom.evapRate, ss.swss.evapRate);
        }
        #endif
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
        fNameIn = inFolder + "gwbc.input";
        if (!readGwBCFile(fNameIn, gdom, gbc, par, gw)) {
            if (par.masterproc) {
                std::cerr << RERROR "Unable to read subsurface BC" << std::endl;
                return 0;
            }
        }
        // read subsurface source sink terms
        fNameIn = inFolder + "gwss.input";
        if (!readGwSSFile(fNameIn, gw, gdom, ss, par))   {
            if (par.masterproc) {
                std::cerr << RERROR << " Unable to read subsurface source/sink." << std::endl;
                return 0;
            }
        }
        // read initial conditions
        fNameIn = inFolder;
        if (!setGwState(fNameIn, gw, gdom, state, gbc, par, io)) {
            if (par.masterproc) {
                std::cerr << RERROR "Unable to read subsurface IC" << std::endl;
                return 0;
            }
        }
        Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
            gw.h(iGlob,0) = gw.h(iGlob,1);  gw.wc(iGlob,0) = gw.wc(iGlob,1);
        });
        gmpi.mpi_sendrecv(gw.h, gdom, par);
        gmpi.mpi_sendrecv(gw.wc, gdom, par);
        for (int k = 0; k < gbc.gwbc.size(); k++) {
            gbc.gwbc[k].applyHBC(gw, gdom, par);
        }
		// scalar transport 
		#if SERGHEI_SUBSURFACE_TRANSPORT
        // read subsurface source sink terms
        fNameIn = inFolder + "scalar.input";
        if (!readScalarFile(fNameIn, gw, gdom, par))   {
            if (par.masterproc) {
                std::cerr << RERROR << " Unable to read scalar transport file." << std::endl;
                return 0;
            }
        }
        // read initial conditions
        fNameIn = inFolder;
        if (!setScalarField(fNameIn, gw, gdom, par, io)) {
            if (par.masterproc) {
                std::cerr << RERROR "Unable to read scalar field IC" << std::endl;
                return 0;
            }
        }
		#endif
		// integrator
        gint.initialize(gw, gdom, gbc.gwbc, ss.gwss);
        gint.integrate(gw, gdom, gbc.gwbc, ss.gwss);
        // output
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
        gdom.dz_multiplier = -999;
        gdom.dt_init = -999;
        gdom.dt_max = -999;
        gdom.nSoilID = -999;
        gdom.gw_scheme = -999;
        gdom.aev = -999;
        gdom.async = -999;
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
                    else if ( !strcmp( "height"    , pline.key.c_str() ) ) { pline.value >> gdom.thickH; }
                    else if ( !strcmp( "dz_multiplier"    , pline.key.c_str() ) ) { pline.value >> gdom.dz_multiplier; }
                    else if ( !strcmp( "dt_init"    , pline.key.c_str() ) ) { pline.value >> gdom.dt_init; }
                    else if ( !strcmp( "dt_max"    , pline.key.c_str() ) ) { pline.value >> gdom.dt_max; }
                    else if ( !strcmp( "nSoilID", pline.key.c_str()))   {pline.value >> gdom.nSoilID;}
                    else if ( !strcmp( "gw_scheme"    , pline.key.c_str() ) ) { pline.value >> gdom.gw_scheme; }
                    else if ( !strcmp( "aev"    , pline.key.c_str() ) ) { pline.value >> gdom.aev; }
                    else if ( !strcmp( "async"    , pline.key.c_str() ) ) { pline.value >> gdom.async; }
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
        if (gdom.dz_multiplier   == -999) { if (par.masterproc) std::cerr << RERROR "key " << "dz_multiplier" << " not set."; exit(-1); }
        if (gdom.thickH   == -999) { if (par.masterproc) std::cerr << RERROR "key " << "height" << " not set."; exit(-1); }
        if (gdom.dt_init   == -999) { if (par.masterproc) std::cerr << RERROR "key " << "dt_init" << " not set."; exit(-1); }
        if (gdom.dt_max    == -999) { if (par.masterproc) std::cerr << RERROR "key " << "dt_max" << " not set."; exit(-1); }
        if (gdom.gw_scheme    == -999) { if (par.masterproc) std::cerr << RERROR "key " << "gw_scheme" << " not set."; exit(-1); }
        if (gdom.aev    == -999) { if (par.masterproc) std::cerr << RERROR "key " << "aev" << " not set."; exit(-1); }
        if (gdom.async    == -999) { if (par.masterproc) std::cerr << RERROR "key " << "async" << " not set."; exit(-1); }
      // Print out the values
        if (par.masterproc) {
            std::cerr << BDASH "Number of grids (nz)  : "  << gdom.nz_glob    << "\n";
            std::cerr << BDASH "Domain thickness : "  << gdom.thickH    << "\n";
            std::cerr << BDASH "Maximum dt   : "  << gdom.dt_max    << "\n";
            std::cerr << BDASH "Asynchronous SW-GW coupling   : "  << gdom.async    << "\n";
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
                            gw.vgTable(packIndicesUniformGrid(gdom.nSoilID, gw.nVGparam, i, 6)) = std::stof(head);
                        }
            		    gw.vgTable(packIndicesUniformGrid(gdom.nSoilID, gw.nVGparam, gdom.nSoilID - 1, 6)) = std::stof(tail);
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
                            gw.vgTable(packIndicesUniformGrid(gdom.nSoilID, gw.nVGparam, i, 4)) = std::stof(head);
                        }
            		    gw.vgTable(packIndicesUniformGrid(gdom.nSoilID, gw.nVGparam, gdom.nSoilID - 1, 4)) = std::stof(tail);
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
                            gw.vgTable(packIndicesUniformGrid(gdom.nSoilID, gw.nVGparam, i, 0)) = std::stof(head);
                        }
            		    gw.vgTable(packIndicesUniformGrid(gdom.nSoilID, gw.nVGparam, gdom.nSoilID - 1, 0)) = std::stof(tail);
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
                            gw.vgTable(packIndicesUniformGrid(gdom.nSoilID, gw.nVGparam, i, 1)) = std::stof(head);
                        }
                        gw.vgTable(packIndicesUniformGrid(gdom.nSoilID, gw.nVGparam, gdom.nSoilID - 1, 1)) = std::stof(tail);
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
                            gw.vgTable(packIndicesUniformGrid(gdom.nSoilID, gw.nVGparam, i, 3)) = std::stof(head);
                        }
            		    gw.vgTable(packIndicesUniformGrid(gdom.nSoilID, gw.nVGparam, gdom.nSoilID - 1, 3)) = std::stof(tail);
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
                            gw.vgTable(packIndicesUniformGrid(gdom.nSoilID, gw.nVGparam, i, 2)) = std::stof(head);
                        }
                        gw.vgTable(packIndicesUniformGrid(gdom.nSoilID, gw.nVGparam, gdom.nSoilID - 1, 2)) = std::stof(tail);
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
        for (idx = 0; idx < gdom.nCell; idx++)  {
            gdom.unpackIndices(idx, kk, jj, ii);
            // global index for this rank
            iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
            // global index for the entire domain
            ii2 = kk*gdom.nx_glob*gdom.ny_glob + (par.j_beg+jj)*(gdom.nx_glob)+par.i_beg+ii;
            gw.soilID(iGlob) = tmpVar(ii2);
        }
        // halo cells
        // for (idx = 0; idx < gdom.nhalo; idx++)  {
        //     gw.soilID(gdom.hpair(idx,0)) = gw.soilID(gdom.hpair(idx,1));
        // }
        if (par.masterproc)   {std::cerr<<GOK "Soil ID  set\n";}
        return 1;
    }

    /*
        Read boundary condition for the subsurface
    */
    int readGwBCFile(std::string fNameIn, GwDomain &gdom, SubsurfaceBoundaries &gbc, Parallel &par, GwState &gw) {
        std::ifstream fInStream(fNameIn);
        std::string dir;
        std::vector<std::string> polygonFile;
        std::vector<std::string> fullPathPoly;
        std::vector<std::string> tsFile, bcFile;
        std::string line;
        PsLn pline;
        int nPoly, nPoly3D;
        dir = fNameIn.substr(0, fNameIn.length() - 10); // 10 chars equivalent to "gwbc.input" to get the dir
        int bccount = 0, bccountFound = 0, ibc = -2, ndata = 0, hasbcfile;
        real val;
        // Read the gwbc.input file
        if (fInStream.is_open()) {
            while (std::getline(fInStream, line)) {
                pline.line = line;
                pline.parse();
                if(!pline.key.empty()){
                    // we should read the number of boundaries here
    	  			if (!strcmp("bccount", pline.key.c_str()))    {
                        pline.value >> bccount;
                        bccountFound = 1;
                        if (bccount < 1) {
                            std::cout << YEXC << "subbc.input indicates zero external boundaries." << std::endl;
                            return 1;
                        }
                        gbc.gwbc.resize(bccount);
                        gbc.id.resize(bccount);
                        polygonFile.resize(bccount);
    	    			fullPathPoly.resize(bccount);
    	    			tsFile.resize(bccount);
                        bcFile.resize(bccount);
    	    			ibc++;	// ibc should be set to -1
    	  			}
                    else if (!strcmp("id", pline.key.c_str())){
                        ibc++;
                        hasbcfile = 0;
                        if (bccount > 0 && ibc >= 0) {pline.value >> gbc.id[ibc];}
                    }
                    else if(!strcmp("bctype", pline.key.c_str()) && ibc >=0 ) {
                        if(bccount > 0) pline.value >> gbc.gwbc[ibc].bctype;
                    }
                    else if(!strcmp("polygon", pline.key.c_str()) && ibc >=0 ) {
                        if(bccount > 0) pline.value >> polygonFile[ibc];
                    }
                    else if (!strcmp("direction", pline.key.c_str()) && ibc >=0 ) {
                        if (bccount > 0) {
                            // pline.value >> gbc.gwbc[ibc].normalx >> gbc.gwbc[ibc].normaly >> gbc.gwbc[ibc].normalz;
                            pline.value >> gbc.gwbc[ibc].direction;
                            if (gbc.gwbc[ibc].direction == 1 || gbc.gwbc[ibc].direction == 2)   {
                                ndata = gdom.ny_glob * gdom.nz;
                            }
                            else if (gbc.gwbc[ibc].direction == 3 || gbc.gwbc[ibc].direction == 4)  {
                                ndata = gdom.nx_glob * gdom.nz;
                            }
                            else if (gbc.gwbc[ibc].direction == 5 || gbc.gwbc[ibc].direction == 6)  {
                                ndata = gdom.nx_glob * gdom.ny_glob;
                            }
                            else {
                                std::cerr << RERROR << "In gwbc.input: direction must be 1, 2, 3, 4, 5 or 6 " << std::endl;
                            }
                        }
                    }
                    else if(!strcmp("bcvals", pline.key.c_str()) && ibc >=0 ){
                        if (bccount > 0){
                            if (ndata <= 0) {
                                std::cerr << RERROR << "In gwbc.input: direction should on top of bcvals " << std::endl;
                            }
                            gbc.gwbc[ibc].bcvals = realArr("bcvals", ndata);
                            pline.value >> val;
                            for (int idx = 0; idx < ndata; idx++)  {
                                gbc.gwbc[ibc].bcvals(idx) = val;
                            }
                        }
                    }
                    else if(!strcmp("timeseries",pline.key.c_str()) && ibc >=0){
                        pline.value >> tsFile[ibc];
                    }
                    else if(!strcmp("bcfile",pline.key.c_str()) && ibc >=0){
                        pline.value >> bcFile[ibc];
                        hasbcfile = 1;
                        if(bccount > 0 ){
                            if (ndata <= 0) {
                                std::cerr << RERROR << "In gwbc.input: direction should on top of bcfile " << std::endl;
                            }
                            gbc.gwbc[ibc].bcvals = realArr("bcvals", ndata);
                        }
                    }
                    else if (ibc >= 0){
                        if (par.masterproc){
                            std::cerr << RERROR << "In gwbc.input: Key " << pline.key << " not understood." << std::endl;
                            return 0;
                        }
                    }
                    if(ibc < -1){
                        if(par.masterproc){
                            std::cerr << RERROR << "No boundaries defined in subbc.input, number of boundaries not defined, or 'id' key not found." << std::endl;
                            return 0;
                        }
                    }
                }
            } //end while
            fInStream.close();
            if(!bccountFound){
                std::cerr << RERROR << "Number of boundaries not defined in gwbc.input. Please define 'bccount'" << std::endl;
                return 0;
            }
            else {
                if (par.masterproc) {std::cerr<< GOK "Subsurface BC set\n";}
            }
        }
        else {
    	 	if (par.masterproc)   {std::cerr << YEXC << "gwbc.input not found. Default boundaries used." << std::endl;}
        }
        // Read polygon file
        for (int k = 0; k < polygonFile.size(); k ++) {
            fullPathPoly[k] = dir + polygonFile[k];
            std::ifstream fPoly(fullPathPoly[k]);
            // read in kth polygon
            if (fPoly.is_open()) {
                fPoly.ignore(256,' ');
                fPoly >> nPoly;
                realArr xPoly=realArr( "xPoly" , nPoly );
                realArr yPoly=realArr( "yPoly" , nPoly );
                for (int i=0; i<nPoly; i++) {
                    if (!fPoly.fail() && !fPoly.eof()) {
                        fPoly >> xPoly(i) >> yPoly(i);
                        #if SERGHEI_DEBUG_BOUNDARY
                          std::cout << GGD << "subbc polygon " << k << ". Point " << i << "/" << nPoly << "\t" << xPoly(i) << "\t" << yPoly(i) << std::endl;
                        #endif
                    }
                    else {
                        if(par.masterproc){
                            std::cerr<< RERROR "Error reading subsurface boundary polygon file " << k << ": " << fullPathPoly[k] << std::endl;
                            return 0;
                        }
                    }
                }
                if(!gbc.gwbc[k].find_bcells(gw, gbc.id[k], gdom, par, nPoly, xPoly, yPoly)) return 0;
            }
            else{
                if (par.masterproc) {
                    std::cerr<< RERROR "Polygon file " << k << ": " << fullPathPoly[k] << " not found." << std::endl;
            	    return 0;
    		    }
            }
            fPoly.close();
        } // end for read in of the kth polygon

        // Read time series boundary conditions
        for (int k = 0; k < tsFile.size(); k++) {
            std::string fname = dir + tsFile[k];
            std::ifstream fts(fname);
            int ndatat=0, readts=0;
            // Read ts file if a ts boundary exists
  			switch (gbc.gwbc[k].bctype) {
                case SUB_BC_Q_T:
                case SUB_BC_H_T:
                case SUB_BC_WT_T:
                    readts = 1;
  					break;
			}
			if (readts) {
				if(fts.is_open()) {
                    fts.ignore(256,' ');
                    fts >> ndatat;
                    if (ndatat > 0) {gbc.gwbc[k].ts.initialise(ndatat);}
                    for (int i = 0; i < ndatat; i++) {
                        if (!fts.fail() && !fts.eof()) {
                            fts >> gbc.gwbc[k].ts.time(i) >> gbc.gwbc[k].ts.value(i);
                        }
                        else {
                            if(par.masterproc){
                                std::cerr<< RERROR "Error reading timeseries file for boundary " << k << ": " << tsFile[k] << std::endl;
                                return 0;
                            }
                        }
  					} // end for ndata
  					fts.close();
  				}
                else {
					if (par.masterproc) {
						std::cerr << RERROR "Error opening timeseries file " << fname << std::endl;   return 0;
					}
  				}
  			}
        } // end for timeseries files
        // Read spatially distributed boundary conditions
        for (int k = 0; k < bcFile.size(); k++) {
            std::string fname = dir + bcFile[k];
            std::ifstream fbc(fname);
            int nx, ny, nz, readbc=0;
            // Read bc file if a bc boundary exists
            if (bcFile[k].length() > 0) {
                switch (gbc.gwbc[k].bctype) {
                    case SUB_BC_Q_CONST:
                        if (gdom.isRain)    {break;}
                        else {readbc = 1;}
                    case SUB_BC_H_CONST:
                    case SUB_BC_WT_CONST:
                        readbc = 1;
      					break;
    			}
            }
            // struct stat buffer;
            if (hasbcfile && readbc && fbc.good()) {
                // get total data size should be read
				if(fbc.is_open()) {
                    std::cout << GOK << "Reading subsurface boundary file : " << fname << std::endl;
                    fbc.ignore(256,' ');
                    fbc >> nx >> ny >> nz;
                    if (nx * ny * nz != ndata)  {
                        std::cerr << RERROR << nx << ny << nz << ndata << std::endl;
                        std::cerr << RERROR "Error reading bc data file " << bcFile[k] << " nx*ny*nz != ndata! " << std::endl;   return 0;
                    }
                    if (ndata > 0)  {
                        int idx = 0;
                        for (int kk = 0; kk < nz; kk++) {
                            for (int jj = 0; jj < ny; jj++) {
                                for (int ii = 0; ii < nx; ii++) {
                                    if (!fbc.fail() && !fbc.eof()) {
                                        fbc >> gbc.gwbc[k].bcvals(idx);
                                        idx += 1;
                                    }
                                    else    {
                                        if(par.masterproc)  {
                                            std::cerr<< RERROR "Error reading bc file for boundary " << k << ": " << bcFile[k] << std::endl;
                                            return 0;
                                        }
                                    }
                                }
                            }
                        }
                    }
  					fbc.close();
  				}
                else {
					if (par.masterproc) {
						std::cerr << RERROR "Error opening bc file " << fname << std::endl;   return 0;
					}
  				}
  			}
        }
        if (par.masterproc) std::cout << GOK << "Subsurface boundary file parsed and boundaries set" << std::endl;
        return 1;
    }



    /*
        Read source sink terms for the subsurface
    */
    int readGwSSFile(std::string fNameIn, GwState &gw, GwDomain &gdom, SourceSink &ss, Parallel &par) {
        std::ifstream fInStream(fNameIn);
        std::string dir;
        std::vector<std::string> polygonFile;
        std::vector<std::string> fullPathPoly;
        std::vector<std::string> tsFile;
        std::string line;
        PsLn pline;
        int nPoly;
        dir = fNameIn.substr(0, fNameIn.length() - 10); // 10 chars equivalent to "gwbc.input" to get the dir
        int sscount = 0, sscountFound = 0, iss = -2, readts = 0, hasssfile;
        real val;
        // Read the gwbc.input file
        if (fInStream.is_open()) {
            while (std::getline(fInStream, line)) {
                pline.line = line;
                pline.parse();
                if(!pline.key.empty()){
                    // we should read the number of boundaries here
    	  			if (!strcmp("sscount", pline.key.c_str()))    {
                        pline.value >> sscount;
                        sscountFound = 1;
                        if (sscount < 1) {
                            std::cout << YEXC << "gwss.input indicates zero source sinks." << std::endl;
                            return 1;
                        }
                        ss.gwss.resize(sscount);
                        ss.id.resize(sscount);
                        polygonFile.resize(sscount);
    	    			fullPathPoly.resize(sscount);
    	    			tsFile.resize(sscount);
    	    			iss++;	// ibc should be set to -1
    	  			}
                    else if (!strcmp("id", pline.key.c_str())){
                        iss++;
                        hasssfile = 0;
                        if (sscount > 0 && iss >= 0) {pline.value >> ss.id[iss];}
                    }
                    else if(!strcmp("sstype", pline.key.c_str()) && iss >=0 ) {
                        if(sscount > 0) pline.value >> ss.gwss[iss].sstype;
                    }
                    else if(!strcmp("polygon", pline.key.c_str()) && iss >=0 ) {
                        if(sscount > 0) pline.value >> polygonFile[iss];
                    }
                    else if(!strcmp("direction", pline.key.c_str()) && iss >=0 ) {
                        if(sscount > 0) pline.value >> ss.gwss[iss].direction;
                    }
                    else if(!strcmp("ssvals", pline.key.c_str()) && iss >=0 ){
                        if (sscount > 0){
                            ss.gwss[iss].ssvals = realArr("ssvals", 1);
                            pline.value >> val;
                            ss.gwss[iss].ssvals(0) = val;
                        }
                    }
                    else if(!strcmp("shapecorrection", pline.key.c_str()) && iss >=0 ) {
                        if(sscount > 0) pline.value >> ss.gwss[iss].Cpipe;
                    }
                    else if(!strcmp("timeseries",pline.key.c_str()) && iss >=0){
                        pline.value >> tsFile[iss];
                        readts = 1;
                    }
                    else if (iss >= 0){
                        if (par.masterproc){
                            std::cerr << RERROR << "In gwss.input: Key " << pline.key << " not understood." << std::endl;
                            return 0;
                        }
                    }
                    if(iss < -1){
                        if(par.masterproc){
                            std::cerr << RERROR << "No Source Sink defined in gwss.input, or 'id' key not found." << std::endl;
                            return 0;
                        }
                    }
                }
            } //end while
            fInStream.close();
            if(!sscountFound){
                std::cerr << RERROR << "Number of Source Sinks not defined in gwss.input. Please define 'sscount'" << std::endl;
                return 0;
            }
            else {
                if (par.masterproc) {std::cerr<< GOK "Subsurface Source/Sink set\n";}
            }
        }
        else {
    	 	if (par.masterproc)   {std::cerr << YEXC << "gwss.input not found. No source/sinks used." << std::endl;}
        }
        // Read polygon 3D file for source/sink terms
        for (int k = 0; k < polygonFile.size(); k ++) {
            fullPathPoly[k] = dir + polygonFile[k];
            std::ifstream fPoly(fullPathPoly[k]);
            // read in kth polygon
            if (fPoly.is_open()) {
                fPoly.ignore(256,' ');
                fPoly >> nPoly;
                realArr xPoly=realArr( "xPoly" , nPoly );
                realArr yPoly=realArr( "yPoly" , nPoly );
                realArr zPoly=realArr( "zPoly" , nPoly );
                for (int i=0; i<nPoly; i++) {
                    if (!fPoly.fail() && !fPoly.eof()) {
                        fPoly >> xPoly(i) >> yPoly(i) >> zPoly(i);
                    }
                    else {
                        if(par.masterproc){
                            std::cerr<< RERROR "Error reading subsurface internal polygon file " << k << ": " << fullPathPoly[k] << std::endl;
                            return 0;
                        }
                    }
                }
                if(!ss.gwss[k].find_icells(gw, ss.id[k], gdom, par, nPoly, xPoly, yPoly, zPoly)) return 0;
            }
            else{
                if (par.masterproc) {
                    std::cerr<< RERROR "Polygon file " << k << ": " << fullPathPoly[k] << " not found." << std::endl;
            	    return 0;
    		    }
            }
            fPoly.close();
        } // end for read in of the kth polygon
        // Read time series boundary conditions
        for (int k = 0; k < tsFile.size(); k++) {
            std::string fname = dir + tsFile[k];
            std::ifstream fts(fname);
            int ndatat=0;
            // Read metero data and calculate ET using Penman-Monteith equation
            if (ss.gwss[k].sstype == 0) {
                if (readts) {
    				if(fts.is_open()) {
                        gdom.hasET = 1;
                        real lat, dayoffset, albedo;
                        TimeSeries wind, solar, rhmax, rhmin, tmax, tmin, crop, lai;
                        // number of data
                        fts.ignore(256,' ');
                        fts >> ndatat;
                        // latitude
                        fts.ignore(256,' ');
                        fts >> lat;
                        // day offset
                        fts.ignore(256,' ');
                        fts >> dayoffset;
                        // albedo
                        fts.ignore(256,' ');
                        fts >> albedo;
                        // allocate time series
                        wind.value = realArr ("w",     ndatat);
                        solar.value = realArr ("s",     ndatat);
                        rhmax.value = realArr ("rmax",     ndatat);
                        rhmin.value = realArr ("rmin",     ndatat);
                        tmax.value = realArr ("tmax",     ndatat);
                        tmin.value = realArr ("tmin",     ndatat);
                        crop.value = realArr ("c",     ndatat);
                        lai.value = realArr ("l",     ndatat);
                        if (ndatat > 0) {
                            ss.gwss[k].ts.initialise(ndatat);
                            ss.gwss[k].evap.initialise(ndatat);
                            ss.gwss[k].tran.initialise(ndatat);
                        }
                        for (int i = 0; i < ndatat; i++) {
                            if (!fts.fail() && !fts.eof()) {
                                fts >> ss.gwss[k].ts.time(i);
                                fts >> tmax.value(i);
                                fts >> tmin.value(i);
                                fts >> rhmax.value(i);
                                fts >> rhmin.value(i);
                                fts >> wind.value(i);
                                fts >> solar.value(i);
                                fts >> crop.value(i);
                                fts >> lai.value(i);
                            }
                            else {
                                if(par.masterproc){
                                    std::cerr<< RERROR "Error reading timeseries file for evapotranspiration " << k << ": " << tsFile[k] << std::endl;
                                    return 0;
                                }
                            }
      					} // end for ndata
      					fts.close();
                        // Use PM equation to calculate evapotranspiration
                        //  Note that the ET flux is in m/s
                        // ASSUMPTIONS:
                        //      wind measured at 2m elevation
                        //      atmosphere pressure = 101kPa
                        //      soil radiation is negligible
                        //      albedo = 0.23 (for grass)
                        real gamma = 0.665e-3 * 101.3;
                        real gsc = 0.082;
                        real sigma = 4.903e-9;
                        real tavg, es, ea, e0max, e0min, delta, rn, ra, rso, rns, rnl, rg, day, dr, d, ws, nume, deno;
                        for (int ii = 0; ii < ndatat; ii++)    {
                            // day in year
                            day = dayoffset + ss.gwss[k].ts.time(ii);
                            dr = 1.0 + 0.033 * cos(2.0 * 3.14 * day / 365);
                            d = 0.409 * sin(2.0 * 3.14 * day / 365 - 1.39);
                            ws = acos(-tan(lat)*tan(d));
                            // temperature and humidity
                            tavg = 0.5 * (tmax.value(ii) + tmin.value(ii));
                            e0max = 0.6108 * exp(17.27*tmax.value(ii)/(237.3+tmax.value(ii)));
                            e0min = 0.6108 * exp(17.27*tmin.value(ii)/(237.3+tmin.value(ii)));
                            es = 0.5 * (e0max + e0min);
                            ea = 0.5 * (e0max * rhmax.value(ii) + e0min * rhmin.value(ii));
                            delta = 4098 * (0.6108 * exp(17.27*tavg/(237.3+tavg))) / pow((tavg+237.3),2.0);
                            // net radiation
                            rns = solar.value(ii) * (1 - albedo);
                            ra = (24*60*gsc*dr/3.14)*(ws*sin(lat)*sin(d) + cos(lat)*cos(d)*sin(ws));
                            rso = 0.75*ra;
                            rnl = 0.25*sigma*(pow(tmax.value(ii),4)+pow(tmin.value(ii),4))*(0.34-0.14*sqrt(ea))*(1.35*solar.value(ii)/rso-0.35);
                            rn = rns - rnl;
                            rg = 0.0;
                            // integrate into the PM equation
                            nume = 0.408*delta*(rn-rg) + 900*gamma*wind.value(ii)*(es-ea)/(tavg+273.0);
                            deno = delta + gamma*(1.0+0.34*wind.value(ii));
                            ss.gwss[k].ts.value(ii) = nume / deno;
                            // from mm/d to m/s
                            ss.gwss[k].ts.time(ii) = ss.gwss[k].ts.time(ii) * 86400.0;
                            ss.gwss[k].ts.value(ii) = -ss.gwss[k].ts.value(ii) * crop.value(ii) / 1e3 / 86400.0;
                            // split evaporation and transpiration
                            ss.gwss[k].evap.time(ii) = ss.gwss[k].ts.time(ii);
                            ss.gwss[k].evap.value(ii) = ss.gwss[k].ts.value(ii) * (1.0 - lai.value(ii));
                            ss.gwss[k].tran.time(ii) = ss.gwss[k].ts.time(ii);
                            ss.gwss[k].tran.value(ii) = ss.gwss[k].ts.value(ii) * lai.value(ii);
                            // distributed along the root depth
                            ss.gwss[k].tran.value(ii) = ss.gwss[k].tran.value(ii) / ss.gwss[k].ndepth;
                        }
      				}
                    else {
                        gdom.hasET = 0;
    					if (par.masterproc) {
    						std::cerr << RERROR "Error opening timeseries file " << fname << std::endl;   return 0;
    					}
      				}
      			}
            }
            else {
                if (readts) {
    				if(fts.is_open()) {
                        fts.ignore(256,' ');
                        fts >> ndatat;
                        if (ndatat > 0) {ss.gwss[k].ts.initialise(ndatat);}
                        for (int i = 0; i < ndatat; i++) {
                            if (!fts.fail() && !fts.eof()) {
                                fts >> ss.gwss[k].ts.time(i) >> ss.gwss[k].ts.value(i);
                            }
                            else {
                                if(par.masterproc){
                                    std::cerr<< RERROR "Error reading timeseries file for boundary " << k << ": " << tsFile[k] << std::endl;
                                    return 0;
                                }
                            }
      					} // end for ndata
      					fts.close();
      				}
                    else {
    					if (par.masterproc) {
    						std::cerr << RERROR "Error opening timeseries file " << fname << std::endl;   return 0;
    					}
      				}
      			}
            }
        } // end for timeseries files
        if (par.masterproc) std::cout << GOK << "Subsurface source/sink file parsed set" << std::endl;
        return 1;
    }

    // read subsurface initial conditions
    int setGwState(std::string inFolder, GwState &gw, GwDomain &gdom, State &state, SubsurfaceBoundaries &gbc, Parallel &par, FileIO &io){
        int ii, jj, kk, idx, ivg, iGlob, iGlobSW;
        real wcs, wcr, n, alpha;
        std::ifstream fInStream(inFolder + "subsurface.input");
        std::string line;
        PsLn pline;
        //SubsurfaceModel sub;
        std::string tempStr;
        // gw.initialMode="saturated";
        gw.initialMode = IC_SAT;
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
                }
            }
        }
        else {
            if (par.masterproc){
                std::cerr << RERROR "File " << inFolder + "subsurface.input" << " not found" << std::endl; return 0;
            }
        }
        // initialize the primary variables
        for (iGlob = 0; iGlob < gdom.nCellMem; iGlob++)   {
            gw.wc(iGlob,0) = 0.0;   gw.wc(iGlob,1) = 0.0;
            gw.h(iGlob,0) = 0.0;    gw.h(iGlob,1) = 0.0;    gw.wc(iGlob,2) = 0.0;
        }
        // read initial condition from file
        if (gw.initialMode == IC_H) {
            tempStr = "head.input";
            readGwICFile(tempStr, inFolder, gw, gdom, par);
        }
        else if (gw.initialMode == IC_WC){
            tempStr = "theta.input";
            readGwICFile(tempStr, inFolder, gw, gdom, par);
        }
        else if (gw.initialMode == IC_WT){
            tempStr = "wt.input";
            readGwICFile(tempStr, inFolder, gw, gdom, par);
        }
        else if (gw.initialMode == IC_SAT){
            for (iGlob = 0; iGlob < gdom.nCellMem; iGlob++)   {
                gdom.unpackIndicesGw(iGlob, gdom.nzhc, gdom.nyhc, gdom.nxhc, kk, jj, ii);
                iGlobSW = packIndicesUniformGrid(gdom.nyhc, gdom.nxhc, jj, ii);
                idx = gw.soilID(iGlob) * gw.nVGparam;
                wcs = gw.vgTable(idx + 2);
                gw.wc(iGlob,0) = wcs;   gw.wc(iGlob,1) = wcs;
                gw.h(iGlob,1) = state.h(iGlobSW) + (kk-0.5)*gdom.dz(iGlob);
                gw.h(iGlob,0) = gw.h(iGlob,1);
            }
        }
        else {
            if (par.masterproc) {
                std::cerr << RERROR "Initial mode must be IC_SAT(1), IC_H(2), IC_WC(3) or IC_WT(4)!" << std::endl; return 0;
            }
        }
        if (par.masterproc){std::cerr<<GOK "Subsurface initial condition set" << std::endl;}
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
            for (idx = 0; idx < gdom.nCell; idx++)    {
                gdom.unpackIndices(idx, kk, jj, ii);
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
            for (idx = 0; idx < gdom.nCell; idx++)    {
                gdom.unpackIndices(idx, kk, jj, ii);
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
            for (idx = 0; idx < gdom.nCell; idx++)    {
                gdom.unpackIndices(idx, kk, jj, ii);
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

};

#endif
