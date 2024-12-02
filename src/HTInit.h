/* -*- mode: c++; c-default-style: "linux" -*- */

#ifndef _HT_INIT_H_
#define _HT_INIT_H_

#include "define.h"
#include "Domain.h"
#include "FileIO.h"
#include "GwDomain.h"
#include "GwState.h"
#include "Parallel.h"
#include "Parser.h"
#include "SourceSink.h"
#include "RTBC.h"
#include "HTBC.h"

class HTInit : public Initializer
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
    int initialize_ht(HTState &ht, RTState &rt, GwState &gw, GwDomain &gdom, HTSubsurfaceBoundaries &htgbc, GwMPI &gmpi, Parallel &par, FileIO &io, SourceSink &ss, std::string inFolder, std::string outFolder)
    {
        int flag = -1;
        int ii, jj, kk, idx, iGlob, iGlobSW;
        // Read soil heat input file
        std::string fNameIn = inFolder + "heat.input";
        if (!readHTFile(fNameIn, ht, par))
        {
            std::cerr << GOK << " Reading in heat transport input file failed." << std::endl;
            return 0;
        
        }
        // allocate subsurface state variable
        ht.allocate(gdom);
        // gmpi.allocate(gdom);//20240510添加
        // READ BOUNDARY CONDITIONS FOR REACTIVE TRANSPORT

        fNameIn = inFolder + "htgwbc.input";
        if (!readHTBCFile(fNameIn, gdom, htgbc, par, ht, gw))

        {
            if (par.masterproc)
            {
                std::cerr << RERROR "Unable to read heat BC" << std::endl;
                return 0;
            }
        }

        // 读取溶质物理化学反应相关参数
      /*  if (rt.ReactionModule == 1 || rt.ReactionModule == 2)
        {
            fNameIn = inFolder + "reaction.input";
            if (!readRTReactionFile(fNameIn, rt, par))
            {
                if (par.masterproc)
                {
                    std::cerr << RERROR "Unable to read transport Reaction" << std::endl;
                    return 0;
                }
            }
        }*/
        // READ SOURCE/SINK TERMS FOR REACTIVE TRANSPORT
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
        if (!setHtState(fNameIn, ht, gw, gdom, par, io))
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "Unable to read HTM transport IC" << std::endl;
                return 0;
            }
        }
        Kokkos::parallel_for(
            gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) { 
                ht.T(iGlob, 0) = ht.T(iGlob, 1);                
                });
        gmpi.mpi_sendrecv(ht.T, gdom, par);

        // 施加边界条件
        for (int k = 0; k < htgbc.htgwbc.size(); k++)
        {
            htgbc.htgwbc[k].applyTemperatureBC(ht, gw, gdom, par);
        }
        /*
        // WRITE OUTPUT FOR REACTIVE TRANSPORT
        io.outputIniRT(outFolder);
        */
        io.outputIniHT(ht, gdom, par, outFolder);

        flag = 1;
        return flag;
    }

    /*
        Read input file for the subsurface
    */
    int readHTFile(std::string fNameIn, HTState &ht, Parallel &par)
    {
        // Initialize all read-in values to -999
        ht.HT_Aquifer_initialMode = -999;       
        ht.ht_scheme = -999;
        ht.Cw = -999;
        ht.b1 = -999;
        ht.b2 = -999;
        ht.b3 = -999;
        ht.theta_n = -999;
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
                    if (!strcmp("HT_Aquifer_initialMode", pline.key.c_str()))
                    {
                        pline.value >> ht.HT_Aquifer_initialMode;
                        //std::cout << "ht.HT_Aquifer_initialMode: " << ht.HT_Aquifer_initialMode << std::endl; // 添加打印语句
                    }                                     
                    else if (!strcmp("ht_scheme", pline.key.c_str()))
                    {
                        pline.value >> ht.ht_scheme;
                    }                   
                    else if (!strcmp("Cw", pline.key.c_str()))
                    {
                        pline.value >> ht.Cw;
                    }
                    else if (!strcmp("b1", pline.key.c_str()))
                    {
                        pline.value >> ht.b1;
                    }
                    else if (!strcmp("b2", pline.key.c_str()))
                    {
                        pline.value >> ht.b2;
                    }
                    else if (!strcmp("b3", pline.key.c_str()))
                    {
                        pline.value >> ht.b3;
                    }
                    else if (!strcmp("theta_n", pline.key.c_str()))
                    {
                        pline.value >> ht.theta_n;
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
        if (ht.HT_Aquifer_initialMode == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "key "
                          << "HT_Aquifer_initialMode"
                          << " not set."
                          << std::endl;
            exit(-1);
        }                     
        if (ht.ht_scheme == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "key "
                          << "ht_scheme"
                          << " not set."
                          << std::endl;
            exit(-1);
        }
        if (ht.Cw == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "key "
                          << "Cw"
                          << " not set."
                          << std::endl;
            exit(-1);
        }
        if (ht.b1 == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "key "
                          << "b1"
                          << " not set."
                          << std::endl;
            exit(-1);
        }
        if (ht.b2 == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "key "
                          << "b2"
                          << " not set."
                          << std::endl;
            exit(-1);
        }
        if (ht.b3 == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "key "
                          << "b3"
                          << " not set."
                          << std::endl;
            exit(-1);
        }
        if (ht.theta_n == -999)
        {
            if (par.masterproc)
                std::cerr << RERROR "key "
                          << "theta_n"
                          << " not set."
                          << std::endl;
            exit(-1);
        }
                       
        if (par.masterproc)
        {
            std::cerr << GOK "Heat transport parameters read\n";
        }
        return 1;
    }
    /*
    读取反应相关参数
    */
   /* int readRTReactionFile(std::string fNameIn, RTState &rt, Parallel &par)
    {
        // Initialize all read-in values to -999
        // rt.ReactionModule = -999;
        rt.Kd = -999;
        // rt.lambda = -999;
        rt.AdsorptionDesorptionModel = -999;
        rt.rho_b = -999;
        rt.Kf = -999;
        rt.Nf = -999;
        rt.alpha_D = -999;
        rt.beta_D = -999;
        rt.beta = -999;
        rt.f = -999;
        rt.lambda_1 = -999;
        rt.lambda_2 = -999;
        rt.k1 = -999;
        rt.k2 = -999;
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
                    if (!strcmp("lambda", pline.key.c_str())) // 衰减系数
                    {
                        pline.value >> rt.lambda;
                    }

                    else if (!strcmp("rho_b", pline.key.c_str())) 
                    {
                        pline.value >> rt.rho_b;
                    }
                    
                    
                    else if (!strcmp("AdsorptionDesorptionModel", pline.key.c_str()))
                    {
                        pline.value >> rt.AdsorptionDesorptionModel;
                    }

                    else if (!strcmp("Kd", pline.key.c_str()))
                    {
                        pline.value >> rt.Kd;
                        // std::cout << "---Kd: " << rt.Kd << std::endl; // 添加打印语句
                    }

                    else if (!strcmp("Kf", pline.key.c_str()))
                    {
                        pline.value >> rt.Kf;
                        // std::cout << "---Kf: " << rt.Kf << std::endl; // 添加打印语句
                    }

                    else if (!strcmp("Nf", pline.key.c_str()))
                    {
                        pline.value >> rt.Nf;
                        // std::cout << "---Nf: " << rt.Nf << std::endl; // 添加打印语句
                    }

                    else if (!strcmp("alpha_D", pline.key.c_str()))
                    {
                        pline.value >> rt.alpha_D;
                    }

                    else if (!strcmp("beta_D", pline.key.c_str()))
                    {
                        pline.value >> rt.beta_D;
                    }
                    else if (!strcmp("beta", pline.key.c_str()))
                    {
                        pline.value >> rt.beta;
                    }
                    else if (!strcmp("f", pline.key.c_str()))
                    {
                        pline.value >> rt.f;
                    }
                    else if (!strcmp("lambda_1", pline.key.c_str()))
                    {
                        pline.value >> rt.lambda_1;
                    }
                    else if (!strcmp("lambda_2", pline.key.c_str()))
                    {
                        pline.value >> rt.lambda_2;
                    }
                    else if (!strcmp("k1", pline.key.c_str()))
                    {
                        pline.value >> rt.k1;
                    }
                    else if (!strcmp("k2", pline.key.c_str()))
                    {
                        pline.value >> rt.k2;
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
                        if (rt.AdsorptionDesorptionModel == 1)
                        {
                            std::cout << GOK "AdsorptionDesorptionModel:Linear adsorption mode\n"; 
                                if (rt.Kd == -999)
                                {
                                    if (par.masterproc)
                                        std::cerr << RERROR "key "
                                                << "Kd"
                                                << " not set."
                                                << std::endl;
                                    exit(-1);
                                }
                                                       
                        }
                        else if (rt.AdsorptionDesorptionModel == 2)
                        {
                            std::cout << GOK "RTM AdsorptionDesorptionModel: Freundlich isothermal adsorption \n";
 
                                if (rt.Kf == -999){
                                    if (par.masterproc)
                                        std::cerr << RERROR "key "
                                                  << "Kf"
                                                  << " not set."
                                                  << std::endl;
                                    exit(-1);
                                }
                            
                                else if (rt.Nf == -999){
                                    if (par.masterproc)
                                        std::cerr << RERROR "key"
                                                  << "Nf"
                                                  << " not set."
                                                  << std::endl;
                                    exit(-1);
                                }                               
                            
                        }
                        else if (rt.AdsorptionDesorptionModel == 3)
                        {
                            std::cout << GOK "RTM AdsorptionDesorptionModel: Langmuir isothermal adsorption \n";

                                if (rt.alpha_D == -999){
                                    if (par.masterproc)
                                        std::cerr << RERROR "key"
                                                  << " alpha_D"
                                                  << " not set."
                                                  << std::endl;
                                    exit(-1);
                                }                                   
                                else if (rt.beta_D == -999){
                                    if (par.masterproc)
                                        std::cerr << RERROR "key"
                                                  << "beta_D"
                                                  << " not set."
                                                  << std::endl;
                                    exit(-1);
                                }                                    
                            
                        }

                        else if (rt.AdsorptionDesorptionModel == 4)
                        {
                            std::cout << GOK "RTM AdsorptionDesorptionModel: Nonequilibrium sorption \n";

                                if (rt.beta == -999){
                                    if (par.masterproc)
                                        std::cerr << RERROR "key"
                                                  << " beta"
                                                  << " not set."
                                                  << std::endl;
                                    exit(-1);
                                }
                                else if (rt.f == -999){
                                    if (par.masterproc)
                                        std::cerr << RERROR "key"
                                                  << "f"
                                                  << " not set."
                                                  << std::endl;
                                    exit(-1);
                                } 
                                else if (rt.rho_b == -999){
                                    if (par.masterproc)
                                        std::cerr << RERROR "key"
                                                  << "rho_b"
                                                  << " not set."
                                                  << std::endl;
                                    exit(-1);
                                }                                                                    
                                else if (rt.lambda_1 == -999){
                                    if (par.masterproc)
                                        std::cerr << RERROR "key"
                                                  << "lambda_1"
                                                  << " not set."
                                                  << std::endl;
                                    exit(-1);
                                }  
                                else if (rt.lambda_2 == -999){
                                    if (par.masterproc)
                                        std::cerr << RERROR "key"
                                                  << "lambda_2"
                                                  << " not set."
                                                  << std::endl;
                                    exit(-1);
                                }                                      
                                else if (rt.rho_b == -999){
                                    if (par.masterproc)
                                        std::cerr << RERROR "key"
                                                  << "rho_b"
                                                  << " not set."
                                                  << std::endl;
                                    exit(-1);
                                }  
                                else if (rt.Kd == -999){
                                    if (par.masterproc)
                                        std::cerr << RERROR "key"
                                                  << "Kd"
                                                  << " not set."
                                                  << std::endl;
                                    exit(-1);
                                }
                                    

                        } 
                        if (rt.k1 == -999)
                        {
                                    if (par.masterproc)
                                        std::cerr << RERROR "key"
                                                  << "k1"
                                                  << " not set."
                                                  << std::endl;
                                    exit(-1);
                                }
                        if (rt.k2 == -999)
                        {
                                    if (par.masterproc)
                                        std::cerr << RERROR "key"
                                                  << "k2"
                                                  << " not set."
                                                  << std::endl;
                                    exit(-1);
                                }                         
                    // }

                    // if (!strcmp("Bd", pline.key.c_str()))
                    // {
                    //     pline.value >> rt.Bd;
                    // }

        // Test to make sure all values were initialized

        // if (rt.Bd == -999)
        // {
        //     if (par.masterproc)
        //         std::cerr << RERROR "key "
        //                   << "Bd"
        //                   << " not set.";
        //     exit(-1);
        // }
        // if (rt.Kd == -999)
        // {
        //     if (par.masterproc)
        //         std::cerr << RERROR "key "
        //                   << "Kd"
        //                   << " not set.";
        //     exit(-1);
        // }
        // if (rt.lambda == -999)
        // {
        //     if (par.masterproc)
        //         std::cerr << RERROR "key "
        //                   << "lambda"
        //                   << " not set.";
        //     exit(-1);
        // }

        if (par.masterproc)
        {
            std::cerr << GOK "RTM Reaction parameters read\n";
        }
        return 1;
    }*/

    /*
        READ BOUNDARY CONDITIONS file
    */
    int readHTBCFile(std::string fNameIn, GwDomain &gdom, HTSubsurfaceBoundaries &htgbc, Parallel &par, HTState &ht, GwState &gw)
    {
        std::ifstream fInStream(fNameIn);
        std::string dir;
        std::vector<std::string> polygonFile;
        std::vector<std::string> fullPathPoly;
        std::vector<std::string> tsFile, bcFile;
        std::string line;
        PsLn pline;
        int nPoly, nPoly3D;
        dir = fNameIn.substr(0, fNameIn.length() - 12); // 12 chars equivalent to "htgwbc.input" to get the dir
        int bccount = 0, bccountFound = 0, ibc = -2, ndata = 0, hasbcfile;
        real val;
        // Read the gwbc.input file
        if (fInStream.is_open())
        {
            while (std::getline(fInStream, line))
            {
                pline.line = line;
                pline.parse();
                if (!pline.key.empty())
                {
                    // we should read the number of boundaries here
                    if (!strcmp("bccount", pline.key.c_str()))
                    {
                        // strcmp函数用于比较两个字符串是否相等,相等则返回0
                        pline.value >> bccount;
                        bccountFound = 1;
                        if (bccount < 1)
                        {
                            std::cout << YEXC << "HTM htgwbc.input indicates zero external boundaries." << std::endl;
                            return 1;
                        }
                        htgbc.htgwbc.resize(bccount);
                        htgbc.id.resize(bccount);
                        polygonFile.resize(bccount);
                        fullPathPoly.resize(bccount);
                        tsFile.resize(bccount);
                        bcFile.resize(bccount);
                        ibc++; // ibc should be set to -1
                    }
                    else if (!strcmp("id", pline.key.c_str()))
                    {
                        ibc++; // 变为0,利用ibc作为边界添加不同id的索引，第一个id关键词为索引0，下一个id关键词为索引1
                        hasbcfile = 0;
                        if (bccount > 0 && ibc >= 0)
                        {
                            pline.value >> htgbc.id[ibc];
                        }
                    }
                    else if (!strcmp("htbctype", pline.key.c_str()) && ibc >= 0)
                    {
                        if (bccount > 0)
                            pline.value >> htgbc.htgwbc[ibc].htbctype;
                    }
                    else if (!strcmp("polygon", pline.key.c_str()) && ibc >= 0)
                    {
                        if (bccount > 0)
                            pline.value >> polygonFile[ibc];
                    }
                    else if (!strcmp("direction", pline.key.c_str()) && ibc >= 0)
                    {
                        if (bccount > 0)
                        {
                            // pline.value >> gbc.gwbc[ibc].normalx >> gbc.gwbc[ibc].normaly >> gbc.gwbc[ibc].normalz;
                            pline.value >> htgbc.htgwbc[ibc].direction;
                            if (htgbc.htgwbc[ibc].direction == 1 || htgbc.htgwbc[ibc].direction == 2)
                            {
                                ndata = gdom.ny_glob * gdom.nz;
                            }
                            else if (htgbc.htgwbc[ibc].direction == 3 || htgbc.htgwbc[ibc].direction == 4)
                            {
                                ndata = gdom.nx_glob * gdom.nz;
                            }
                            else if (htgbc.htgwbc[ibc].direction == 5 || htgbc.htgwbc[ibc].direction == 6)
                            {
                                ndata = gdom.nx_glob * gdom.ny_glob;
                            }
                            else
                            {
                                std::cerr << RERROR << "In htgwbc.input: direction must be 1, 2, 3, 4, 5 or 6 " << std::endl;
                            }
                        }
                    }
                    else if (!strcmp("bcvals", pline.key.c_str()) && ibc >= 0)
                    {
                        if (bccount > 0)
                        {
                            if (ndata <= 0)
                            {
                                std::cerr << RERROR << "In htgwbc.input: direction should on top of bcvals " << std::endl;
                            }
                            htgbc.htgwbc[ibc].bcvals = realArr("bcvals", ndata);
                            pline.value >> val;
                            for (int idx = 0; idx < ndata; idx++)
                            {
                                htgbc.htgwbc[ibc].bcvals(idx) = val;
                            }
                        }
                    }
                    else if (!strcmp("timeseries", pline.key.c_str()) && ibc >= 0)
                    {
                        pline.value >> tsFile[ibc];
                    }
                    else if (!strcmp("bcfile", pline.key.c_str()) && ibc >= 0)
                    {
                        pline.value >> bcFile[ibc];
                        hasbcfile = 1;
                        if (bccount > 0)
                        {
                            if (ndata <= 0)
                            {
                                std::cerr << RERROR << "In htgwbc.input: direction should on top of bcfile " << std::endl;
                            }
                            htgbc.htgwbc[ibc].bcvals = realArr("bcvals", ndata);
                        }
                    }
                    else if (ibc >= 0)
                    {
                        if (par.masterproc)
                        {
                            std::cerr << RERROR << "In htgwbc.input: Key " << pline.key << " not understood." << std::endl;
                            return 0;
                        }
                    }
                    if (ibc < -1)
                    {
                        if (par.masterproc)
                        {
                            std::cerr << RERROR << "No boundaries defined in htgwbc.input, number of boundaries not defined, or 'id' key not found." << std::endl;
                            return 0;
                        }
                    }
                }
            } // end while
            fInStream.close();
            if (!bccountFound)
            {
                std::cerr << RERROR << "Number of boundaries not defined in htgwbc.input. Please define 'bccount'" << std::endl;
                return 0;
            }
            else
            {
                if (par.masterproc)
                {
                    std::cerr << GOK "HTM Subsurface BC set\n";
                }
            }
        }
        else
        {
            if (par.masterproc)
            {
                std::cerr << YEXC << "htgwbc.input not found. Default boundaries used." << std::endl;
            }
        }
        // Read polygon file
        for (int k = 0; k < polygonFile.size(); k++)
        {
            fullPathPoly[k] = dir + polygonFile[k];

            std::ifstream fPoly(fullPathPoly[k]);
            // std::cout <<fullPathPoly[k]<< std::endl;
            // read in kth polygon
            if (fPoly.is_open())
            {
                fPoly.ignore(256, ' ');
                fPoly >> nPoly;
                realArr xPoly = realArr("xPoly", nPoly);
                realArr yPoly = realArr("yPoly", nPoly);
                for (int i = 0; i < nPoly; i++)
                {
                    if (!fPoly.fail() && !fPoly.eof())
                    {
                        fPoly >> xPoly(i) >> yPoly(i);
#if SERGHEI_DEBUG_BOUNDARY
                        std::cout << GGD << "HTM subbc polygon " << k << ". Point " << i << "/" << nPoly << "\t" << xPoly(i) << "\t" << yPoly(i) << std::endl;
#endif
                    }
                    else
                    {
                        if (par.masterproc)
                        {
                            std::cerr << RERROR "Error reading HTM subsurface boundary polygon file " << k << ": " << fullPathPoly[k] << std::endl;
                            return 0;
                        }
                    }
                }
                if (!htgbc.htgwbc[k].find_bcells(gw, htgbc.id[k], gdom, par, nPoly, xPoly, yPoly))
                    return 0;
            }
            else
            {
                if (par.masterproc)
                {
                    std::cerr << RERROR "HTM Polygon file " << k << ": " << fullPathPoly[k] << " not found." << std::endl;
                    return 0;
                }
            }
            fPoly.close();

        } // end for read in of the kth polygon

        // Read time series boundary conditions
        for (int k = 0; k < tsFile.size(); k++)
        {
            std::string fname = dir + tsFile[k];
            std::ifstream fts(fname);
            int ndatat = 0, readts = 0;
            // Read ts file if a ts boundary exists
            switch (htgbc.htgwbc[k].htbctype)
            {
            case SUB_HT_BC_Dirichlet_T:
            case SUB_HT_BC_Neumann_T:
            case SUB_HT_BC_Cauchy_T:
                readts = 1;
                break;
            }
            if (readts)
            {
                if (fts.is_open())
                {
                    fts.ignore(256, ' ');
                    fts >> ndatat;
                    if (ndatat > 0)
                    {
                        htgbc.htgwbc[k].ts.initialise(ndatat);
                    }
                    for (int i = 0; i < ndatat; i++)
                    {
                        if (!fts.fail() && !fts.eof())
                        {
                            fts >> htgbc.htgwbc[k].ts.time(i) >> htgbc.htgwbc[k].ts.value(i);
                        }
                        else
                        {
                            if (par.masterproc)
                            {
                                std::cerr << RERROR "Error reading timeseries file for HTM boundary " << k << ": " << tsFile[k] << std::endl;
                                return 0;
                            }
                        }
                    } // end for ndata
                    fts.close();
                }
                else
                {
                    if (par.masterproc)
                    {
                        std::cerr << RERROR "Error opening HTM timeseries file " << fname << std::endl;
                        return 0;
                    }
                }
            }
        } // end for timeseries files
        // Read spatially distributed boundary conditions
        for (int k = 0; k < bcFile.size(); k++)
        {
            std::string fname = dir + bcFile[k];
            std::ifstream fbc(fname);
            int nx, ny, nz, readbc = 0;
            // Read bc file if a bc boundary exists
            if (bcFile[k].length() > 0)
            {
                switch (htgbc.htgwbc[k].htbctype)
                {
                case SUB_HT_BC_Neumann_CONST:
                    if (gdom.isRain)
                    {
                        break;
                    }
                    else
                    {
                        readbc = 1;
                    }
                case SUB_HT_BC_Dirichlet_CONST:
                case SUB_HT_BC_Cauchy_CONST:
                    readbc = 1;
                    break;
                }
            }
            // struct stat buffer;
            if (hasbcfile && readbc && fbc.good())
            {
                // get total data size should be read
                if (fbc.is_open())
                {
                    std::cout << GOK << "Reading HTM subsurface boundary file : " << fname << std::endl;
                    fbc.ignore(256, ' ');
                    fbc >> nx >> ny >> nz;
                    if (nx * ny * nz != ndata)
                    {
                        std::cerr << RERROR << nx << ny << nz << ndata << std::endl;
                        std::cerr << RERROR "Error reading HTM bc data file " << bcFile[k] << " nx*ny*nz != ndata! " << std::endl;
                        return 0;
                    }
                    if (ndata > 0)
                    {
                        int idx = 0;
                        for (int kk = 0; kk < nz; kk++)
                        {
                            for (int jj = 0; jj < ny; jj++)
                            {
                                for (int ii = 0; ii < nx; ii++)
                                {
                                    if (!fbc.fail() && !fbc.eof())
                                    {
                                        fbc >> htgbc.htgwbc[k].bcvals(idx);
                                        idx += 1;
                                    }
                                    else
                                    {
                                        if (par.masterproc)
                                        {
                                            std::cerr << RERROR "Error reading HTM bc file for boundary " << k << ": " << bcFile[k] << std::endl;
                                            return 0;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    fbc.close();
                }
                else
                {
                    if (par.masterproc)
                    {
                        std::cerr << RERROR "Error opening HTM bc file " << fname << std::endl;
                        return 0;
                    }
                }
            }
        }
        if (par.masterproc)
            std::cout << GOK << "HTM Subsurface boundary file parsed and boundaries set" << std::endl;
        return 1;
    }


    // READ SOURCE/SINK TERMS FOR REACTIVE TRANSPORT


    /*
        READ INITIAL CONDITIONS FOR REACTIVE TRANSPORT
    */
/*
    int setRtState(std::string inFolder, RTState &rt, GwState &gw, GwDomain &gdom, Parallel &par, FileIO &io)
    {
        int ii, jj, kk, idx, ivg, iGlob, iGlobSW;

        std::ifstream fInStream(inFolder + "transport.input");
        std::string line;
        PsLn pline;
        // SubsurfaceModel sub;
        std::string tempStr;
        // gw.initialMode="saturated";
        rt.RT_Aquifer_initialMode = IC_H;
        // IC_H =2,如果transport.input文件中没有RTinitialMode这个参数，那么默认为IC_H，即初始化浓度（类似head.input文件格式）
        //  read initial mode and value
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
                        pline.value >> rt.RT_Aquifer_initialMode;
                        // std::cout << "RTinitialMode2: " << rt.RT_Aquifer_initialMode << std::endl; // 添加打印语句
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
        // if (rt.RT_Aquifer_initialMode == IC_H)
        // std::cout << "RTinitialMode3: " << rt.RT_Aquifer_initialMode << std::endl; // 添加打印语句
        if (rt.RT_Aquifer_initialMode == IC_H)
        {
            tempStr = "concen.input";
            readRtICFile(tempStr, inFolder, rt, gdom, par);
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
*/
// 20240814修改
    /*
        READ AQUIFER AND Solid INITIAL CONDITIONS FOR REACTIVE TRANSPORT
    */
    int setHtState(std::string inFolder, HTState &ht, GwState &gw, GwDomain &gdom, Parallel &par, FileIO &io)
    {
        int ii, jj, kk, idx, ivg, iGlob, iGlobSW;

        std::ifstream fInStream(inFolder + "heat.input");
        std::string line;
        PsLn pline;
        // SubsurfaceModel sub;
        std::string tempStr;
        // gw.initialMode="saturated";
        //初始化地下水rt模型中含水层和固相中初始浓度值模式为0
        // rt.RT_Aquifer_initialMode = IC_REACTIVE_TRANSPORT_CON;
        // rt.RT_Solid_initialMode = IC_REACTIVE_TRANSPORT_ZERO;
        //read initial model and value
        // if (fInStream.is_open()){
        //     std::cout << "111111" << std::endl; // 添加打印语句
        //     while (std::getline(fInStream, line)){
        //         pline.line = line;
        //         pline.lowercase();
        //         pline.parse();
        //         // If the line was valid and a key is stored
        //         if (!pline.key.empty()){
        //             // Match the key, and store the value
        //             if (!strcmp("RT_Aquifer_initialMode", pline.key.c_str())){
        //                 pline.value >> rt.RT_Aquifer_initialMode;
        //                 std::cout << "1111RT_Aquifer_initialMode: " << rt.RT_Aquifer_initialMode << std::endl; // 添加打印语句
        //             }
        //             if (!strcmp("RT_Solid_initialMode", pline.key.c_str())){
        //                 pline.value >> rt.RT_Solid_initialMode;
        //             }
        //         }
        //     }
        // }
        // else {
        //     if (par.masterproc){
        //         std::cerr << RERROR "File" << inFolder + "transport.input"
        //                   << " not found" << std::endl;
        //     }
        // }
        // initialize the primary variables
        for (iGlob = 0; iGlob < gdom.nCellMem; iGlob++){
            ht.T(iGlob, 0) = 0.0;
            ht.T(iGlob, 1) = 0.0;            
        }
        //read initial condition from file
        if (ht.HT_Aquifer_initialMode == IC_HEAT_TRANSPORT_TEM){
            tempStr = "temperature.input";
            readHtICFile(tempStr, inFolder, ht, gdom, par);
        }
        else if(ht.HT_Aquifer_initialMode != 0 && ht.HT_Aquifer_initialMode != 1){
            // if (par.masterproc){
                std::cerr << RERROR "HTM Aquifer initial mode must be HTIC_ZERO(1) OR HTIC_TEM(2)!" << std::endl;
            // }
            exit(-1);
        }     

        


        if (par.masterproc){std::cerr<<GOK "Subsurface Heat transport initial condition set" << std::endl;}
        return 1;

    }
// 20240814修改



    /*
        Read rt subsurface initial condition from file
    */

    int readHtICFile(std::string fNameIn, std::string fDirIn, HTState &ht, GwDomain &gdom, Parallel &par)
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
                    std::cerr << RERROR "HTM IC file parameters don't match DEM parameters. Unable to continue\n";
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
                        std::cerr << RERROR "Error reading HTM IC file. Not enough data\n";
                        return 0;
                    }
                }
            }
            fInStream.close();
        }
        else if (par.masterproc)
        {
            std::cerr << RERROR "Unable to open " << fNameIn << "\n";
            exit(-1);
        }
        // Copy data into head or water content
        if (!strcmp(fNameIn.c_str(), "temperature.input"))
        {
            for (idx = 0; idx < gdom.nCell; idx++)
            {
                gdom.unpackIndices(idx, kk, jj, ii);
                // gdom.unpackIndicesGw(idx, gdom.nz, gdom.ny, gdom.nx, kk, jj, ii);
                // get global index
                iGlob = (hc + kk) * gdom.nxhc * gdom.nyhc + (hc + jj) * gdom.nxhc + ii + hc;
                // get concentration
                ii2 = kk * gdom.nx_glob * gdom.ny_glob + (par.j_beg + jj) * (gdom.nx_glob) + par.i_beg + ii;
                ht.T(iGlob, 1) = tmpVar(ii2);               

                ht.T(iGlob, 0) = ht.T(iGlob, 1);
            }
        }
        else
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "Error reading HTM temperature IC file. File name might be wrong.\n";
                return 0;
            }
        }
        if (par.masterproc)
        {
            std::cerr << GOK "HTM Subsurface temperature content set\n";
        }
        return 1;
    }


    /*
        Read boundary condition for the Heat transport
    */

    /*---------------------------------------
                   修改ysl
    ------------------------------------------+++*/
};

#endif
