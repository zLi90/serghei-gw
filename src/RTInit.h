/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _RT_INIT_H_
#define _RT_INIT_H_

#include "define.h"
#include "Domain.h"
#include "FileIO.h"
#include "GwDomain.h"
#include "GwState.h"
#include "Parallel.h"
#include "Parser.h"
#include "SourceSink.h"

class RTInit : public Initializer
{

    // Parse lines, same as Parser : ParserLine
    class PsLn
    {
    public:
        std::string line;
        std::string key;
        std::stringstream value;
        void lowercase()
        {
            std::for_each(line.begin(), line.end(), [](char &c)
                          { c = ::tolower(c); });
        }

        void parse()
        {
            // make sure key and value are clean (in case of reuse)
            key.clear();
            value.clear();
            // If the line isn't empty and doesn't begin with a comment specifier, split it based on the colon
            if (!line.empty() && line.find("//", 0) != 0)
            {
                uint splitloc = line.find(':', 0);
                key = line.substr(0, splitloc);
                key.erase(std::remove(key.begin(), key.end(), ' '), key.end());
                key.erase(std::remove(key.begin(), key.end(), '\t'), key.end());
                std::string val = line.substr(splitloc + 1, line.length() - splitloc);
                size_t splitter = val.find("//", 0);
                std::string strloc;
                if (splitter != std::string::npos)
                {
                    strloc = val.substr(0, splitter);
                }
                else
                {
                    strloc = val;
                }
                value.clear();
                value.str(strloc);
            }
        }

        void print() { std::cout << "line: " << line << "\tkey: " << key << "\tvalue: " << value.str() << std::endl; }
    };

public:
    int initialize_rt(RTState &rt, GwState &gw, GwDomain &gdom, GwMPI &gmpi, Parallel &par, FileIO &io, SourceSink &ss, std::string inFolder, std::string outFolder){
        int flag = -1;
        int ii, jj, kk, idx, iGlob, iGlobSW;
        // Read subsurface input file
        std::string fNameIn = inFolder + "transport.input";
        if (!readRTFile(fNameIn, rt, par))
        {
            std::cerr << GOK << " Reading in reactive transport input file failed." << std::endl;
            return 0;
        }
        // allocate subsurface state variable
        rt.allocate(gdom);

        /*+++---------------------------------------
                       修改部分zzb
        ------------------------------------------*/
        // READ BOUNDARY CONDITIONS FOR REACTIVE TRANSPORT
       
        // fNameIn = inFolder + "rtbc.input";
        // if (!readRTBCFile(fNameIn))
        // {
        //     if (par.masterproc)
        //     {
        //         std::cerr << RERROR "Unable to read transport BC" << std::endl;
        //         return 0;
        //     }
        // }

        // // READ SOURCE/SINK TERMS FOR REACTIVE TRANSPORT
        // fNameIn = inFolder + "rtss.input";
        // if (!readRTSSFile(fNameIn))
        // {
        //     if (par.masterproc)
        //     {
        //         std::cerr << RERROR << " Unable to read transport source/sink." << std::endl;
        //         return 0;
        //     }
        // }

        // READ INITIAL CONDITIONS FOR REACTIVE TRANSPORT
        fNameIn = inFolder;
        if (!setRtState(fNameIn, rt, gw, gdom, par, io))
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "Unable to read transport IC" << std::endl;
                return 0;
            }
        }
        /*---------------------------------------
                       修改部分zzb
        ------------------------------------------+++*/

        Kokkos::parallel_for(
            gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) { rt.c(iGlob, 0) = rt.c(iGlob, 1); });
        gmpi.mpi_sendrecv(rt.c, gdom, par);
        /*
        // WRITE OUTPUT FOR REACTIVE TRANSPORT
        io.outputIniRT(outFolder);
        */
        io.outputIniRT(rt, gdom, par, outFolder); //修改部分zzb

        flag = 1;
        return flag;
    }

    /*
        Read input file for the subsurface
    */
    int readRTFile(std::string fNameIn, RTState &rt, Parallel &par)
    {
        // Initialize all read-in values to -999
        rt.n_mass = -999;
        rt.d_base = -999;
        rt.RTinitialMode = -999;
        std::string strAux;
        // Read in colon-separated key: value file line by line
        std::ifstream fInStream(fNameIn);
        std::string line;
        PsLn pline;
        if (fInStream.is_open())
        {
            while (std::getline(fInStream, line))
            {
                pline.line = line;
                pline.parse();
                if (!pline.key.empty())
                {
                    if (!strcmp("n_mass", pline.key.c_str()))
                    {
                        pline.value >> rt.n_mass;
                        std::cout << "rt.n_mass: " << rt.n_mass << std::endl; // 添加打印语句
                    }
                    else if (!strcmp("d_base", pline.key.c_str()))
                    {
                        pline.value >> rt.d_base;
                    }
                    else if (!strcmp("RTinitialMode", pline.key.c_str()))
                    {
                        pline.value >> rt.RTinitialMode;
                        
                    }
                }
            }
        }
        else
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "Unable to open " << fNameIn << "\n";
                return 0;
            }
        }
        // Test to make sure all values were initialized
        if (rt.n_mass == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "key "
                          << "n_mass"
                          << " not set.";
            exit(-1);
        }
        if (rt.d_base == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "key "
                          << "d_base"
                          << " not set.";
            exit(-1);
        }
        if (rt.RTinitialMode == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "key "
                          << "RTinitialMode"
                          << " not set.";
            exit(-1);
        }
        if (par.masterproc)
        {
            std::cerr << GOK "Transport parameters read\n";
        }
        return 1;
    }


/* ----------------20240506修改------------------------*/
    // READ BOUNDARY CONDITIONS file
    /*
    int readRTBCFile(std::string fNameIn, GwDomain &gdom, RTSubsurfaceBoundaries &rtgbc, Parallel &par, RTState &rt) {
        std::ifstream fInStream(fNameIn);
        std::string dir;
        std::vector<std::string> polygonFile;
        std::vector<std::string> fullPathPoly;
        std::vector<std::string> tsFile, bcFile;
        std::string line;
        PsLn pline;
        int nPoly, nPoly3D;
        dir = fNameIn.substr(0, fNameIn.length() - 10); // 10 chars equivalent to "rtgwbc.input" to get the dir
        int rtbccount = 0, bccountFound = 0, ibc = -2, ndata = 0, hasbcfile;
        real val;
        // Read the gwbc.input file
        if (fInStream.is_open()) {
            while (std::getline(fInStream, line)) {
                pline.line = line;
                pline.parse();
                if(!pline.key.empty()){
                    // we should read the number of boundaries here
    	  			if (!strcmp("rtbccount", pline.key.c_str()))    {
                        pline.value >> rtbccount;
                        bccountFound = 1;
                        if (rtbccount < 1) {
                            std::cout << YEXC << "RTsubbc.input indicates zero external boundaries." << std::endl;
                            return 1;
                        }
                        rtgbc.rtgwbc.resize(rtbccount);
                        rtgbc.rtid.resize(rtbccount);
                        polygonFile.resize(rtbccount);
    	    			fullPathPoly.resize(rtbccount);
    	    			tsFile.resize(rtbccount);
                        bcFile.resize(rtbccount);
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
/* ----------------20240506修改------------------------*/


    // READ SOURCE/SINK TERMS FOR REACTIVE TRANSPORT

    // READ INITIAL CONDITIONS FOR REACTIVE TRANSPORT
    int setRtState(std::string inFolder, RTState &rt, GwState &gw, GwDomain &gdom, Parallel &par, FileIO &io){
        int ii, jj, kk, idx, ivg, iGlob, iGlobSW;
        
        std::ifstream fInStream(inFolder + "transport.input");
        std::string line;
        PsLn pline;
        // SubsurfaceModel sub;
        std::string tempStr;
        // gw.initialMode="saturated";
        // rt.RTinitialMode = IC_SAT;//define IC_SAT=1
        // read initial mode and value
        if (fInStream.is_open())
        {
            while (std::getline(fInStream, line))
            {
                pline.line = line;
                pline.lowercase();
                pline.parse();
                // If the line was valid and a key is stored
                if (!pline.key.empty())
                {
                    // Match the key, and store the value
                    if (!strcmp("RTinitialmode", pline.key.c_str()))
                    {
                        pline.value >> rt.RTinitialMode;
                        //std::cout << "RTinitialMode2: " << rt.RTinitialMode << std::endl; // 添加打印语句

                    }
                }
            }
        }
        else
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "File " << inFolder + "transport.input"
                          << " not found" << std::endl;
                return 0;
            }
        }
        // initialize the primary variables
        for (iGlob = 0; iGlob < gdom.nCellMem; iGlob++)
        {

            rt.c(iGlob, 0) = 0.0;
            rt.c(iGlob, 1) = 0.0;
        }
        // read initial condition from file
        // if (rt.RTinitialMode == IC_H)
        //std::cout << "RTinitialMode3: " << rt.RTinitialMode << std::endl; // 添加打印语句
        if (rt.RTinitialMode == IC_H)
        {
            tempStr = "concen.input";
            readRtICFile(tempStr, inFolder,  rt, gdom, par);
        }
        else
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "Initial mode must be RTIC_SAT(1), RTIC_H(2), RTIC_WC(3) or RTIC_WT(4)!" << std::endl;
                return 0;
            }
        }
        if (par.masterproc)
        {
            std::cerr << GOK "Subsurface initial condition set" << std::endl;
        }
        // gw.h.modify<dualDbl::host_mirror_space> ();
        // gw.wc.modify<dualDbl::host_mirror_space> ();
        return 1;
    }
    /*
        Read initial condition from file
    */

    int readRtICFile(std::string fNameIn, std::string fDirIn,  RTState &rt, GwDomain &gdom, Parallel &par)
    {
        std::string fname = fDirIn + fNameIn;
        std::ifstream fInStream(fname);
        std::string line;
        int tnx, tny, iGlob, iGlobSW, idx, ivg, ii, jj, kk, ii2;
        real tmp, nodata_value;
        int ndata = gdom.ny_glob * gdom.nx_glob * gdom.nz_glob;
        intArr tmpVar = intArr("var", ndata);
        std::string str;
        if (fInStream.is_open())
        {
            std::getline(fInStream, str, ' ');
            std::getline(fInStream, str);
            std::stringstream(str) >> tnx;
            std::getline(fInStream, str, ' ');
            std::getline(fInStream, str);
            std::stringstream(str) >> tny;
            std::getline(fInStream, str, ' ');
            std::getline(fInStream, str);
            std::stringstream(str) >> nodata_value;
            // compare the values t* with the DEM file just to check if we are using the same values, otherwise error
            if (gdom.ny_glob != tny || gdom.nx_glob != tnx)
            {
                if (par.masterproc)
                {
                    std::cerr << RERROR "rtIC file parameters don't match DEM parameters. Unable to continue\n";
                    if (par.masterproc)
                    {
                        std::cerr << BDASH "nx_glob: " << gdom.nx_glob << tnx << "\n";
                        std::cerr << BDASH "ny_glob: " << gdom.ny_glob << tny << "\n";
                    }
                    return 0;
                }
            }
            // read and store data into a temporary view
            for (int ii = 0; ii < ndata; ii++)
            {
                if (!fInStream.fail() && !fInStream.eof())
                {
                    fInStream >> tmp;
                    tmpVar(ii) = tmp;
                }
                else
                {
                    if (par.masterproc)
                    {
                        std::cerr << RERROR "Error reading rtIC file. Not enough data\n";
                        return 0;
                    }
                }
            }
            fInStream.close();
        }
        // Copy data into head or water content
        if (!strcmp(fNameIn.c_str(), "concen.input"))
        {
            for (idx = 0; idx < gdom.nCell; idx++)
            {
                gdom.unpackIndices(idx, kk, jj, ii);
                // gdom.unpackIndicesGw(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
                // get global index
                iGlob = (hc + kk) * gdom.nxhc * gdom.nyhc + (hc + jj) * gdom.nxhc + ii + hc;
                // get concentration
                ii2 = kk * gdom.nx_glob * gdom.ny_glob + (par.j_beg + jj) * (gdom.nx_glob) + par.i_beg + ii;
                rt.c(iGlob, 1) = tmpVar(ii2);

                rt.c(iGlob, 0) = rt.c(iGlob, 1);
            }
        }
        else
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "Error reading rt head/theta IC file. File name might be wrong.\n";
                return 0;
            }
        }
        if (par.masterproc)
        {
            std::cerr << GOK "rt Subsurface head/water content set\n";
        }
        return 1;
    }
    /*
        Read boundary condition for the Reactive transport
    */



    /*---------------------------------------
                   修改部分zzb
    ------------------------------------------+++*/
};

#endif
