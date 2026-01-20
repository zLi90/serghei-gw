/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _RT_INIT_SW_H_
#define _RT_INIT_SW_H_

#include "define.h"
#include "Domain.h"
#include "BC.h"
#include "FileIO.h"
#include "State.h"
#include "Parallel.h"
#include "Parser.h"
#include "SourceSink.h"
#include "RTInit.h"
#include "RTStateSW.h"

class RTInitSW : public Initializer  
{
    
public:
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


    int initialize_rt(RTStateSW &rtsw, State &state, Domain &dom, SourceSinkData &ss, ExternalBoundaries &ebc, FileIO &io, Parallel &par, std::string inFolder, std::string outFolder)   
    {

        //read transportsw input file
        std::string fNameIn = inFolder + "transportsw.input";

        if (!readRTFileSW(fNameIn, rtsw, par))
        {
            std::cerr << GOK << " Reading in SW reactive transport  input file failed." << std::endl;
            return 0;
        
        }

        // allocate surface state variable
        rtsw.allocate(dom);
    
        // read boundary conditions
        fNameIn = inFolder + "rtswbc.input";
        if (!readRTBCFileSW(fNameIn, dom, ebc, par, state)) {
            if (par.masterproc) {
                std::cerr << RERROR "Unable to read rtsw BC SW File " << std::endl;
                return 0;
            }
        }
        // set initial state of sw solute
        fNameIn = inFolder;
        if (!setRtStateSW(fNameIn, rtsw, dom, par)) {
            if (par.masterproc) {
                std::cerr << RERROR "Unable to set SW solute IC " << std::endl;
                return 0;
            }
        }        
        // apply boundary conditions
        for (int k = 0; k < ebc.rtbc.size(); k++)   {
            ebc.rtbc[k].applyrtbc(rtsw, dom);
        }
        Kokkos::parallel_for("update", dom.nCellMem, KOKKOS_LAMBDA (int iGlob) {
            // state.c(iGlob) = rtsw.c(iGlob,1);
            rtsw.c(iGlob,0) = rtsw.c(iGlob,1); // update c(0) to c(1)
        });
        
        // write output
        io.outputIniRTSW(rtsw, dom,  par, outFolder);
        
        return 1;
    }
    
    
    
/*
    set initial solute concentration
*/
int setRtStateSW(std::string inFolder, RTStateSW &rtsw, Domain &dom, Parallel &par) 
{
    int ii, jj, kk, idx, ivg, iGlob, iGlobSW;

    // std::ifstream fInStream(inFolder + "transportsw.input");
    std::string line;
    PsLn pline;
    // SubsurfaceModel sub;
    std::string tempStr;
    //read initial condition from file

    if (rtsw.RtInitialModeSW == RTIC_CONST)   {
        for (iGlob = 0; iGlob < dom.nCellMem; iGlob++){
            rtsw.c(iGlob, 0) = rtsw.c0;   rtsw.c(iGlob, 1) = rtsw.c0;

            // printf("iGlob: %d, rtsw.c0: %f, rtsw.c(iGlob, 0): %f, rtsw.c(iGlob, 1): %f\n", iGlob, rtsw.c0, rtsw.c(iGlob, 0), rtsw.c(iGlob, 1));
        }
    }
    else if (rtsw.RtInitialModeSW == RTIC_FILE)   {
        std::string fnameIn = inFolder + "concensw.input";
        realArr tmpc = realArr("tmpc", dom.nCellMem);
        if(!readRasterField(fnameIn, dom, par, tmpc)) return 0;
        for (iGlob = 0; iGlob < dom.nCellMem; iGlob++){
            rtsw.c(iGlob, 0) = tmpc(iGlob);   rtsw.c(iGlob, 1) = tmpc(iGlob);
        }
    }
    else {
        if (par.masterproc) {std::cerr << RERROR "Initial Mode for rtsw SW must be IC_CONST=0 or IC_FILE=1! Now " << rtsw.RtInitialModeSW << std::endl;}
        exit(-1);
    }
    if (par.masterproc) {std::cerr << GOK "Initial solute concentration for SW is set" << std::endl;}
    return 1;
}
/*
    Read SW rtsw transport parameters from file
*/
    /*
        Read input file for the subsurface
    */
    int readRTFileSW(std::string fNameIn, RTStateSW &rtsw, Parallel &par)
    {
        // Initialize all read-in values to -999
        rtsw.RtInitialModeSW = -999;
        // rtsw.ReactionModule = -999;
        rtsw.c0 = -999;
        rtsw.n_mass = -999;
        rtsw.InitialCon = -999;
        // rtsw.diffusion_molecular = -999;
        // rtsw.alpha_T = -999;
        // rtsw.alpha_L = -999;
        // rtsw.rt_scheme = -999;
        // rtsw.Up_Weighting_vplus = -999;
        // rtsw.Up_Weighting_vminus = -999;
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
                        pline.value >> rtsw.n_mass;
                        // std::cout << "rtsw.n_mass: " << rtsw.n_mass << std::endl; // 添加打印语句
                        std::cout << BDASH "rtsw.n_mass: " << rtsw.n_mass << std::endl;

                    }
                    // else if (!strcmp("diffusion_molecular", pline.key.c_str()))
                    // {
                    //     pline.value >> rtsw.diffusion_molecular;
                    // }
                    else if (!strcmp("RtInitialModeSW", pline.key.c_str()))
                    {
                        pline.value >> rtsw.RtInitialModeSW;
                        std::cout << BDASH "RtInitialModeSW: " << rtsw.RtInitialModeSW << std::endl;

                    }
                    
                            
                    else if (!strcmp("InitialCon", pline.key.c_str())) 
                    {
                        pline.value >> rtsw.c0;

                    }
                        
                    //  else if (!strcmp("InitialCon", pline.key.c_str()))
                    // {
                    //     pline.value >> rtsw.c0;
                    //     //std::cout << "rtsw.RT_Solid_initialMode: " << rtsw.RT_Solid_initialMode << std::endl; // 添加打印语句
                    // }                   
                    // else if (!strcmp("alpha_T", pline.key.c_str()))
                    // {
                    //     pline.value >> rtsw.alpha_T;
                    // }
                    // else if (!strcmp("alpha_L", pline.key.c_str()))
                    // {
                    //     pline.value >> rtsw.alpha_L;
                    // }
                    // else if (!strcmp("diffusion_molecular", pline.key.c_str()))
                    // {
                    //     pline.value >> rtsw.diffusion_molecular;
                    // }
                    // else if (!strcmp("rt_scheme", pline.key.c_str()))
                    // {
                    //     pline.value >> rtsw.rt_scheme;
                    // }
                    
                    // else if (!strcmp("Up_Weighting_vplus", pline.key.c_str()))
                    // {
                    //     pline.value >> rtsw.Up_Weighting_vplus;
                    // }
                    // else if (!strcmp("Up_Weighting_vminus", pline.key.c_str()))
                    // {
                    //     pline.value >> rtsw.Up_Weighting_vminus;
                    // }
                    // else if (!strcmp("ReactionModule", pline.key.c_str()))
                    // {
                    //     pline.value >> rtsw.ReactionModule;
                    // }
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
        if (rtsw.n_mass == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "key "
                          << "n_mass"
                          << " not set."
                          << std::endl;
            exit(-1);
        }
        // if (rtsw.diffusion_molecular == -999)
        // {
        //     if (par.masterproc)
        //         std::cerr << RERROR "key "
        //                   << "diffusion_molecular"
        //                   << " not set."
        //                   << std::endl;
        //     exit(-1);
        // }
        
        if (rtsw.RtInitialModeSW == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "key "
                          << "RtInitialModeSW"
                          << " not set."
                          << std::endl;
            exit(-1);
            
        }
        if (rtsw.RtInitialModeSW != -999){
                    if (rtsw.c0 == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "key "
                          << "InitialCon"
                          << " not set."
                          << std::endl;
            exit(-1);
        }
        }

        // if (rtsw.alpha_T == -999)
        // {
        //     if (par.masterproc)
        //         std::cerr << RERROR "key "
        //                   << "alpha_T"
        //                   << " not set."
        //                   << std::endl;
        //     exit(-1);
        // }
        // if (rtsw.alpha_L == -999)
        // {
        //     if (par.masterproc)
        //         std::cerr << RERROR "key "
        //                   << "alpha_L"
        //                   << " not set."
        //                   << std::endl;
        //     exit(-1);
        // }
        // if (rtsw.diffusion_molecular == -999)
        // {
        //     if (par.masterproc)
        //         std::cerr << RERROR "key "
        //                   << "diffusion_molecular"
        //                   << " not set."
        //                   << std::endl;
        //     exit(-1);
        // }
        // if (rtsw.rt_scheme == -999)
        // {
        //     if (par.masterproc)
        //         std::cerr << RERROR "key "
        //                   << "rt_scheme"
        //                   << " not set."
        //                   << std::endl;
        //     exit(-1);
        // }
        // if (rtsw.Up_Weighting_vplus == -999)
        // {
        //     if (par.masterproc)
        //         std::cerr << RERROR "key "
        //                   << "Up_Weighting_vplus"
        //                   << " not set."
        //                   << std::endl;
        //     exit(-1);
        // }
        // if (rtsw.Up_Weighting_vminus == -999)
        // {
        //     if (par.masterproc)
        //         std::cerr << RERROR "key "
        //                   << "Up_Weighting_vminus"
        //                   << " not set."
        //                   << std::endl;
        //     exit(-1);
        // }
        // if (rtsw.ReactionModule == -999)
        // {
        //     if (par.masterproc)
        //         std::cerr << RERROR "key "
        //                   << "ReactionModule"
        //                   << " not set."
        //                   << std::endl;
        //     exit(-1);
        // }
        
        if (par.masterproc)
        {
            std::cerr << GOK "Transport parameters read\n";
        }
        return 1;
    }


/*
    read SW rtsw boundary conditions from file
*/
int readRTBCFileSW(std::string fNameIn, Domain &dom, ExternalBoundaries &ebc, Parallel &par, State &state) 
{
    std::ifstream fInStream(fNameIn);
    std::string dir;
    std::vector<std::string> timeseriesFile;
    std::string line;
    PsLn pline;
    dir = fNameIn.substr(0, fNameIn.length() - 12); // 12 chars equivalent to "rtswbc.input" to get the dir

    int ibc = -1;
    int bccount = 0;
    int bccountFound = 0;
    if(fInStream.is_open()){
        int bccount = ebc.extbc.size();
        ebc.rtbc.resize(bccount);
        timeseriesFile.resize(bccount);
        
        // copy values from ExtBC
    	for (int k = 0; k < ebc.extbc.size(); k++) {
        	ebc.rtbc[k].ncellsBC = ebc.extbc[k].ncellsBC;
        	ebc.rtbc[k].bcells = intArr("bcells",ebc.rtbc[k].ncellsBC);
        	Kokkos::deep_copy(ebc.rtbc[k].bcells, ebc.extbc[k].bcells);
    	}
        
        while (std::getline(fInStream, line)) {
            pline.line = line;
            pline.parse();
		    #if SERGHEI_DEBUG_BOUNDARY
		    pline.print();
		    #endif
            if(!pline.key.empty())  {
                if(!strcmp("id",pline.key.c_str()))    {ibc++;}
                else if (!strcmp("rtbctype", pline.key.c_str()) && ibc >=0 ) {if(bccount > 0) pline.value >> ebc.rtbc[ibc].bctype;}
                else if (!strcmp("direction", pline.key.c_str()) && ibc >=0 ) {if(bccount > 0) pline.value >> ebc.rtbc[ibc].normalx >> ebc.rtbc[ibc].normaly;}
                else if (!strcmp("bcvals", pline.key.c_str()) && ibc >=0 ){
                    ebc.rtbc[ibc].bcvals = realArr("bcvals", 1);
                    pline.value >> ebc.rtbc[ibc].bcvals(0);
                }
                else if(!strcmp("timeseries", pline.key.c_str()) && ibc >=0){
                    pline.value >> timeseriesFile[ibc];
                }
                else if(ibc >= 0){
                    if(par.masterproc){
                        std::cerr << RERROR << "In rtswbc.input: Key " << pline.key << " not understood." << std::endl;
                        return 0;
                    }
                }
            }
        } //end while
        fInStream.close();
    }
    else{
        if(par.masterproc)  std::cerr << YEXC << "rtswbc.input not found. Default boundaries used." << std::endl;
    }
    // normalise BC normal vectors by magnitude
    for (int k = 0; k < ebc.rtbc.size(); k++) {
        real mod = sqrt(ebc.rtbc[k].normalx*ebc.rtbc[k].normalx +ebc.rtbc[k].normaly*ebc.rtbc[k].normaly);
        ebc.rtbc[k].normalx /= mod;
        ebc.rtbc[k].normaly /= mod;
    }
    // read timeseries files
    for (int k = 0; k < timeseriesFile.size(); k++) {
        std::string fname = dir + timeseriesFile[k];
        std::ifstream fbc(fname);
        int ndata=0;
        // read in kth polygon
        int readts=0;
        switch (ebc.rtbc[k].bctype)
          
        {
            case SW_RT_BC_Dirichlet_T:
            case SW_RT_BC_Neumann_T:
            case SW_RT_BC_Cauchy_T:
            readts = 1;
            break;
        }
        // read timeseries file 
        if (readts)   {
            if (fbc.is_open())  {
                fbc.ignore(256,' ');
                fbc >> ndata;
                if (ndata > 0)  {
                    ebc.rtbc[k].hydrograph.initialise(ndata);
                }
				for (int i=0; i<ndata; i++) {
					if (!fbc.fail() && !fbc.eof()) {
						fbc >> ebc.rtbc[k].hydrograph.time(i) >> ebc.rtbc[k].hydrograph.value(i);
					}
					else {
						if (par.masterproc){
							std::cerr<< RERROR "Error reading SW solute timesries file for boundary " << k << ": " << timeseriesFile[k] << std::endl;
							return 0;
						}
					}
				} // end for ndata
				fbc.close();
            }
            else    {
                if (par.masterproc) {
                    std::cerr << RERROR "Error opening timeseries for SW transport boundary file " << fname << std::endl;
                    return 0;
                }
            }
        }
    } // end for hydrograph files
    if (par.masterproc) std::cout << GOK << "rtsw SW boundary file parsed and boundaries set" << std::endl;
    return 1;
}

};
#endif
