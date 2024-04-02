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


class RTInit : public Initializer{

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

    int initialize_rt(RTState &rt, GwState &gw, GwDomain &gdom, Parallel &par, FileIO &io, SourceSink &ss, std::string inFolder, std::string outFolder) {
        int flag = -1;
        int ii, jj, kk, idx, iGlob, iGlobSW;
        // Read subsurface input file
        std::string fNameIn = inFolder + "transport.input";
        if (!readRTFile(fNameIn, rt))    {
            std::cerr << GOK << " Reading in reactive transport input file failed." << std::endl;   return 0;
        }
        // allocate subsurface state variable
		rt.allocate(gdom);
		/*
		// READ BOUNDARY CONDITIONS FOR REACTIVE TRANSPORT
        fNameIn = inFolder + "rtbc.input";
        if (!readRTBCFile(fNameIn)) {
            if (par.masterproc) {
                std::cerr << RERROR "Unable to read transport BC" << std::endl;
                return 0;
            }
        }
		*/
		/*
		// READ SOURCE/SINK TERMS FOR REACTIVE TRANSPORT
        fNameIn = inFolder + "rtss.input";
        if (!readRTSSFile(fNameIn))   {
            if (par.masterproc) {
                std::cerr << RERROR << " Unable to read transport source/sink." << std::endl;
                return 0;
            }
        }
		*/
		/*
		// READ INITIAL CONDITIONS FOR REACTIVE TRANSPORT
        fNameIn = inFolder;
        if (!setRTState(fNameIn)) {
            if (par.masterproc) {
                std::cerr << RERROR "Unable to read transport IC" << std::endl;
                return 0;
            }
        }
		*/
        Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {rt.c(iGlob,0) = rt.c(iGlob,1);});
        gmpi.mpi_sendrecv(rt.c, gdom, par);
		/*
		// WRITE OUTPUT FOR REACTIVE TRANSPORT
        io.outputIniRT(outFolder);
		*/
        flag = 1;
        return flag;
    }


    /*
        Read input file for the subsurface
    */
    int readRTFile(std::string fNameIn, RTState &rt) {
        // Initialize all read-in values to -999
		rt.n_mass = -999;
		rt.d_base = -999;		
        std::string strAux;
        // Read in colon-separated key: value file line by line
        std::ifstream fInStream(fNameIn);
        std::string line;
        PsLn pline;
        if (fInStream.is_open()) {
            while (std::getline(fInStream, line)) {
                pline.line = line;
                pline.parse();
                if (!pline.key.empty())  {
                    if      ( !strcmp( "n_mass" , pline.key.c_str() ) ) { pline.value >> rt.n_mass; }
                    else if ( !strcmp( "d_base"    , pline.key.c_str() ) ) { pline.value >> rt.d_base; }
                }
            }
        }
        else    {
            if (par.masterproc) {
                std::cerr<< RERROR "Unable to open " << fNameIn << "\n";    return 0;
            }
        }
        // Test to make sure all values were initialized
        if (rt.n_mass   == -999) { if (par.masterproc) std::cerr << RERROR "key " << "n_mass" << " not set."; exit(-1); }
        if (rt.d_base   == -999) { if (par.masterproc) std::cerr << RERROR "key " << "d_base" << " not set."; exit(-1); }
        if (par.masterproc)   {std::cerr<< GOK "Transport parameters read\n";}
        return 1;
    }

};

#endif
