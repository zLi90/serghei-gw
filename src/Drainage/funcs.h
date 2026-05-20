#ifndef FUNCS_H
#define FUNCS_H

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include "enums.h"
#include "globals.h"
#include "DrainageDywave.h"
#include "../define.h"
#include "xsect.h"
#include "DrainageState.h"

XXsect xxsect;

class ReadFileFuncs{
public:
    real dt;
    #define UCHAR(x) (((x) >= 'a' && (x) <= 'z') ? ((x)&~32) : (x))
 

    inline void project_readInput(std::string in, Node &Tnode, Link &Tlink, Conduit &Tconduit, Outfall &Toutfall)
    {
        input_countObjects(in);
        createObjects(Tnode, Tlink, Tconduit, Toutfall);
        input_readData(in, Tnode, Tlink, Tconduit, Toutfall);
    }

    inline void project_validate(Node& Tnode, Link& Tlink, Conduit& Tconduit, Outfall& Toutfall)
    {
        int i;
        for ( i=0; i<Nobjects[NODE]; i++) Tnode.oldDepth(i) = Tnode.fullDepth(i);
        for (i = 0; i < Nobjects[NODE]; i++) {
        if (Tnode.fullDepth(i) > 0.0 && Tnode.pondedArea(i) > 0.0) {
            Tnode.fullVolume(i) = Tnode.fullDepth(i) * Tnode.pondedArea(i);
            }
        }
        for ( i=0; i<Nobjects[LINK]; i++) link_validate(i, Tnode, Tlink, Toutfall, Tconduit );
    }

    inline void input_countObjects(std::string in)
    {

        std::string fNameIn = in + "Drainage.inp";
        std::ifstream fInStream(fNameIn);  
        std::string line;
        std::string currentSection;


        if (!fInStream.is_open()) {
            std::cerr << "Error opening file.\n";
        }

        while (std::getline(fInStream, line)) {
            if (line.find("[JUNCTIONS]") != std::string::npos) {
            currentSection = "NODES";
        } else if (line.find("[INLETS]") != std::string::npos) {
            currentSection = "NODES";
        } else if (line.find("[CONDUITS]") != std::string::npos) {
            currentSection = "LINKS";
        } else if (line.find(";;") != std::string::npos || line.empty()) {
            continue;
        } else {
            if (currentSection == "NODES") {
                Nobjects[NODE]++;
            }
            else if (currentSection == "LINKS")
            {
                Nobjects[LINK]++;
                Nlinks[CONDUIT]++;
            }
        }
            
        }

        fInStream.close();
    }

    inline int findmatch(char *s, const char *keyword[])
    
    {
    int i = 0;
    while (keyword[i] != NULL)
    {
        if (match(s, keyword[i])) return(i);
        i++;
    }
    return(-1);
    }

    inline int  match(char *str, const char *substr)

    {
        int i,j,k;
        if (!substr[0]) return(0);
        for (k = 0; str[k]; k++)
        {
            if (str[k] != ' ') break;
        }
        for (i = k,j = 0; substr[j]; i++,j++)
        {
            if (!str[i] || UCHAR(str[i]) != UCHAR(substr[j])) return(0);
        }
        return(1);
    }

    
    inline size_t sstrncpy(char *dest, const char *src, size_t n)
    {
        size_t offset = 0;  
        if (n > 0) {
            while (src[offset] != '\0' && offset < n) {
                dest[offset] = src[offset];
                offset++;
            }
        }
        dest[offset] = '\0';
        return std::char_traits<char>::length(dest);  
    }


    inline void createObjects(Node &Tnode, Link &Tlink, Conduit &Tconduit, Outfall &Toutfall)

    {
        //===================Node===========================
        Tnode.xcoor = realArr("N_xcoor", Nobjects[NODE]);
        Tnode.ycoor = realArr("N_ycoor", Nobjects[NODE]);
        Tnode.typee = intArr("N_type", Nobjects[NODE]);
        Tnode.head = realArr("N_head", Nobjects[NODE]);
        Tnode.subIndex = intArr("N_subIndex", Nobjects[NODE]);
        Tnode.rptFlag = realArr("N_rptFlag", Nobjects[NODE]);
        Tnode.invertElev = realArr("N_invertElev", Nobjects[NODE]);
        Tnode.initDepth = realArr("N_initDepth", Nobjects[NODE]);
        Tnode.fullDepth = realArr("N_fullDepth", Nobjects[NODE]);
        Tnode.surDepth = realArr("N_surDepth", Nobjects[NODE]);
        Tnode.pondedArea = realArr("N_pondedArea", Nobjects[NODE]);
        Tnode.degree = realArr("N_degree", Nobjects[NODE]);
        Tnode.inlet = realArr("N_inlet", Nobjects[NODE]);
        Tnode.updated = realArr("N_updated", Nobjects[NODE]);
        Tnode.crownElev = realArr("N_crownElev", Nobjects[NODE]);
        Tnode.inflow = realArr("N_inflow", Nobjects[NODE]);
        Tnode.outflow = realArr("N_outflow", Nobjects[NODE]);
        Tnode.losses = realArr("N_losses", Nobjects[NODE]);
        Tnode.oldVolume = realArr("N_oldVolume", Nobjects[NODE]);
        Tnode.newVolume = realArr("N_newVolume", Nobjects[NODE]);
        Tnode.fullVolume = realArr("N_fullVolume", Nobjects[NODE]);
        Tnode.overflow = realArr("N_overflow", Nobjects[NODE]);
        Tnode.oldDepth = realArr("N_oldDepth", Nobjects[NODE]);
        Tnode.newDepth = realArr("N_newDepth", Nobjects[NODE]);
        Tnode.oldLatFlow = realArr("N_oldLatFlow", Nobjects[NODE]);
        Tnode.newLatFlow = realArr("N_newLatFlow", Nobjects[NODE]);
        Tnode.oldFlowInflow = realArr("N_oldFlowInflow", Nobjects[NODE]);
        Tnode.oldNetInflow = realArr("N_oldNetInflow", Nobjects[NODE]);
        Tnode.apiExtInflow = realArr("N_apiExtInflow", Nobjects[NODE]);
        Tnode.SDinflow = realArr("N_SDinflow", Nobjects[NODE]);
        Tnode.SDoutflow = realArr("N_SDoutflow", Nobjects[NODE]);
        Tnode.converged = realArr("N_converged", Nobjects[NODE]);
        Tnode.newSurfArea = realArr("N_newSurfArea", Nobjects[NODE]);
        Tnode.oldSurfArea = realArr("N_oldSurfArea", Nobjects[NODE]);
        Tnode.sumdqdh = realArr("N_sumdqdh", Nobjects[NODE]);
        Tnode.dYdT = realArr("N_dYdT", Nobjects[NODE]);
        Tnode.sealed = intArr("N_sealed", Nobjects[NODE]);
        Tnode.toNode = intArr("N_toNode", Nobjects[NODE]);
        Tnode.connectedInlet = intArr("N_connectedInlet", Nobjects[NODE]);
        Tnode.numInlet = intArr("N_numInlet", Nobjects[NODE]);
        Tnode.inletIndex = intArr("N_inletIndex", Nobjects[NODE]);
        Tnode.backflow = realArr("N_backflow", Nobjects[NODE]);
        Tnode.backflowRatio = realArr("N_backflowRatio", Nobjects[NODE]);
        // ==========================================


        //===================Link===========================
        Tlink.typee = intArr("L_type", Nobjects[LINK]);
        Tlink.qFull = realArr("L_qFull", Nobjects[LINK]);
        Tlink.subIndex = intArr("L_subIndex", Nobjects[LINK]);
        Tlink.rptFlag = realArr("L_rptFlag", Nobjects[LINK]);
        Tlink.node1 = intArr("L_node1", Nobjects[LINK]);
        Tlink.node2 = intArr("L_node2", Nobjects[LINK]);
        Tlink.offset1 = realArr("L_offset1", Nobjects[LINK]);
        Tlink.offset2 = realArr("L_offset2", Nobjects[LINK]);
        Tlink.yFull = realArr("L_yFull", Nobjects[LINK]);
        Tlink.wMax = realArr("L_wMax", Nobjects[LINK]);
        Tlink.ywMax = realArr("L_ywMax", Nobjects[LINK]);
        Tlink.aFull = realArr("L_aFull", Nobjects[LINK]);
        Tlink.rFull = realArr("L_rFull", Nobjects[LINK]);
        Tlink.sFull = realArr("L_sFull", Nobjects[LINK]);
        Tlink.sMax = realArr("L_sMax", Nobjects[LINK]);
        Tlink.s = realArr("L_s", Nobjects[LINK]);
        Tlink.qc = realArr("L_qc", Nobjects[LINK]);
        Tlink.q0 = realArr("L_q0", Nobjects[LINK]);
        Tlink.qLimit = realArr("L_qLimit", Nobjects[LINK]);
        Tlink.cLossInlet = realArr("L_cLossInlet", Nobjects[LINK]);
        Tlink.cLossOutlet = realArr("L_cLossOutlet", Nobjects[LINK]);
        Tlink.cLossAvg = realArr("L_cLossAvg", Nobjects[LINK]);
        Tlink.seepRate = realArr("L_seepRate", Nobjects[LINK]);
        Tlink.hasFlapGate = realArr("L_hasFlapGate", Nobjects[LINK]);
        Tlink.oldFlow = realArr("L_oldFlow", Nobjects[LINK]);
        Tlink.newFlow = realArr("L_newFlow", Nobjects[LINK]);
        Tlink.oldDepth = realArr("L_oldDepth", Nobjects[LINK]);
        Tlink.newDepth = realArr("L_newDepth", Nobjects[LINK]);
        Tlink.oldVolume = realArr("L_oldVolume", Nobjects[LINK]);
        Tlink.newVolume = realArr("L_newVolume", Nobjects[LINK]);
        Tlink.surfArea1 = realArr("L_surfArea1", Nobjects[LINK]);
        Tlink.surfArea2 = realArr("L_surfArea2", Nobjects[LINK]);
        Tlink.flowClass = intArr("L_flowClass", Nobjects[LINK]);
        Tlink.dqdh = realArr("L_dqdh", Nobjects[LINK]);
        Tlink.setting = realArr("L_setting", Nobjects[LINK]);
        Tlink.targetSetting = realArr("L_targetSetting", Nobjects[LINK]);
        Tlink.timeLastSet = realArr("L_timeLastSet", Nobjects[LINK]);
        Tlink.froude = realArr("L_froude", Nobjects[LINK]);
        Tlink.direction = realArr("L_direction", Nobjects[LINK]);
        Tlink.bypassed = realArr("L_bypassed", Nobjects[LINK]);
        Tlink.normalFlow = realArr("L_normalFlow", Nobjects[LINK]);
        Tlink.inletControl = realArr("L_inletControl", Nobjects[LINK]);
        // ==========================================

        //===================Conduit===========================
        Tconduit.length = realArr("C_length", Nlinks[CONDUIT]);
        Tconduit.modLength = realArr("C_modLength", Nlinks[CONDUIT]);
        Tconduit.roughness = realArr("C_roughness", Nlinks[CONDUIT]);
        Tconduit.slope = realArr("C_slope", Nlinks[CONDUIT]);
        Tconduit.beta = realArr("C_beta", Nlinks[CONDUIT]);
        Tconduit.qMax = realArr("C_qMax", Nlinks[CONDUIT]);
        Tconduit.a1 = realArr("C_a1", Nlinks[CONDUIT]);
        Tconduit.a2 = realArr("C_a2", Nlinks[CONDUIT]);
        Tconduit.q1 = realArr("C_q1", Nlinks[CONDUIT]);
        Tconduit.q2 = realArr("C_q2", Nlinks[CONDUIT]);
        Tconduit.q1Old = realArr("C_q1Old", Nlinks[CONDUIT]);
        Tconduit.q2Old = realArr("C_q2Old", Nlinks[CONDUIT]);
        Tconduit.evapLossRate = realArr("C_evapLossRate", Nlinks[CONDUIT]);
        Tconduit.seepLossRate = realArr("C_seepLossRate", Nlinks[CONDUIT]);
        Tconduit.capacityLimited = realArr("C_capacityLimited", Nlinks[CONDUIT]);
        Tconduit.superCritical = realArr("C_superCritical", Nlinks[CONDUIT]);
        Tconduit.hasLosses = realArr("C_hasLosses", Nlinks[CONDUIT]);
        Tconduit.fullState = realArr("C_fullState", Nlinks[CONDUIT]);
        Tconduit.roughFactor = realArr("C_roughFactor", Nlinks[CONDUIT]);
        // ==========================================



        for (int j = 0; j < Nobjects[LINK]; j++)
        {
            Tlink.typee(j)   = -1;
            Tlink.cLossInlet(j)   = 0.0;
            Tlink.cLossOutlet(j)  = 0.0;
            Tlink.cLossAvg(j)     = 0.0;
            Tlink.hasFlapGate(j)  = FALSE;
        }

    }

    inline void input_readData(std::string in, Node &Tnode, Link &Tlink, Conduit &Tconduit, Outfall &Toutfall)
    {
        std::string fNameIn = in + "Drainage.inp";
        std::ifstream fInStream(fNameIn); 
        std::string line;
        std::string currentSection;

    if (!fInStream.is_open()) {
        std::cerr << "Error opening file.\n";
    }

    while (std::getline(fInStream, line)) {
        if (line.find("[JUNCTIONS]") != std::string::npos) {
            currentSection = "JUNCTIONS";
        } else if (line.find("[CONDUITS]") != std::string::npos) {
            currentSection = "LINKS";
        } else if (line.find("[INLETS]") != std::string::npos) {
            currentSection = "INLETS";
        } else if (line.find(";;") != std::string::npos || line.empty()) {
            continue; 
        } else {
            if (currentSection == "JUNCTIONS") {
                readNodeData(line, Tnode);
            } else if (currentSection == "INLETS") {
                readInletData(line, Tnode);
            }
            else if (currentSection == "LINKS")
            {
                readLinkData(CONDUIT, line, Tlink, Tconduit, Tnode);
            }
        }
    }

    fInStream.close();
    }

    // inline void readNodeData(const std::string& line, Node &Tnode, Outfall &Toutfall)//zn251020
    inline void readNodeData(const std::string& line, Node &Tnode)//zn251020
    {
        int j = Mobjects[NODE];
        std::istringstream iss(line);
        real elevation, maxDepth, initDepth, surDepth, pondedArea, xcoor, ycoor;
        int id, sealed;
        std::string type;
        iss >> id >> elevation >> maxDepth >> initDepth >> surDepth >> pondedArea >> xcoor >> ycoor >> type >> sealed;
        if (type == "outlet"){Tnode.typee(j) = OUTFALL;}
        else if (type == "junction"){Tnode.typee(j) = JUNCTION;}
        Tnode.sealed(j) = sealed;
        Tnode.invertElev(j) = elevation;
        Tnode.crownElev(j)  = elevation;
        Tnode.initDepth(j)  = 0.0;
        Tnode.newVolume(j)  = 0.0;
        Tnode.fullVolume(j) = 0.0;
        Tnode.fullDepth(j)  = 0.0;
        Tnode.surDepth(j)  = 0.0;
        Tnode.inflow(j) = 0.0;
        Tnode.outflow(j) = 0.0;
        Tnode.fullDepth(j) = maxDepth;
        Tnode.initDepth(j) = initDepth;
        Tnode.surDepth(j)  = surDepth;
        Tnode.pondedArea(j) = pondedArea;
        Tnode.xcoor(j) = xcoor;
        Tnode.ycoor(j) = ycoor;
        Tnode.SDinflow(j) = 0.0;
        Tnode.SDoutflow(j) = 0.0;
    
        Mobjects[NODE]++;
    }

    inline void readInletData(const std::string& line, Node &Tnode)//zn251020
    {
        int j = Mobjects[NODE];
        std::istringstream iss(line);
        real elevation, maxDepth, initDepth, surDepth, pondedArea, xcoor, ycoor;
        int id, toNode;
        iss >> id >> elevation >> maxDepth >> initDepth >> surDepth >> pondedArea >> xcoor >> ycoor >> toNode;
        Tnode.typee(j) = INLET;
        Tnode.toNode(j) = toNode;
        Tnode.connectedInlet(Tnode.toNode(j)) = TRUE;
        Tnode.numInlet(Tnode.toNode(j)) += 1;
        Tnode.inletIndex(Tnode.toNode(j)) = j;
        Tnode.backflow(j) = 0.0;
        Tnode.invertElev(j) = elevation;
        Tnode.crownElev(j)  = elevation;
        Tnode.initDepth(j)  = 0.0;
        Tnode.newVolume(j)  = 0.0;
        Tnode.fullVolume(j) = 0.0;
        Tnode.fullDepth(j)  = 0.0;
        Tnode.surDepth(j)  = 0.0;
        Tnode.inflow(j) = 0.0;
        Tnode.outflow(j) = 0.0;
        Tnode.fullDepth(j) = maxDepth;
        Tnode.initDepth(j) = initDepth;
        Tnode.surDepth(j)  = surDepth;
        Tnode.pondedArea(j) = pondedArea;
        Tnode.xcoor(j) = xcoor;
        Tnode.ycoor(j) = ycoor;
        Tnode.SDinflow(j) = 0.0;
        Tnode.SDoutflow(j) = 0.0;

        Mobjects[NODE]++;
    }

    inline void readLinkData(int type, const std::string& line, Link &Tlink, Conduit &Tconduit, Node &Tnode)
    {
        int j = Mobjects[LINK];
        int k = Mlinks[type];
        std::istringstream iss(line);
        real length, roughness, InOffset, OutOffset, Initflow, maxflow, geom1;
        int id, from, to;
        char shape;
        iss >> id >> from >> to >> length >> roughness >> InOffset >> OutOffset >> Initflow >> maxflow >> geom1;
        Tlink.node1(j) = from;
        Tlink.node2(j) = to;
        Tlink.subIndex(j) = k;
        Tlink.offset1(j) = 0.0;
        Tlink.offset2(j) = 0.0;
        Tlink.q0(j)       = 0.0;
        Tlink.qFull(j)    = 0.0;
        Tlink.setting(j)  = 1.0;
        Tlink.targetSetting(j) = 1.0;
        Tlink.hasFlapGate(j) = 0;     
        Tlink.direction(j)   = 1;
        Tconduit.length(k)    = length;
        // printf("the length is: %f\n", Tconduit.length(k));
        Tconduit.modLength(k) = length;
        Tconduit.roughness(k) = roughness;
        Tlink.qLimit(j)       = maxflow;
        Tlink.yFull(j)        =  geom1;
        Tlink.wMax(j)         =  geom1;
        Tlink.aFull(j)        =  PI / 4.0 * geom1 * geom1;
        Tlink.rFull(j)        =  0.2500 * geom1;
        Tlink.sFull(j)        =  Tlink.aFull(j) * pow(Tlink.rFull(j), 2./3.);
        Tlink.sMax(j)         =  1.08 * Tlink.sFull(j);
        Tlink.ywMax(j)        =  0.5 * Tlink.yFull(j);
        Tlink.typee(j) = 0;
        Tnode.crownElev(Tlink.node1(j)) = max(Tnode.crownElev(Tlink.node1(j)),
                                        Tnode.invertElev(Tlink.node1(j)) + Tlink.yFull(j));
        Tnode.crownElev(Tlink.node2(j)) = max(Tnode.crownElev(Tlink.node2(j)),
                                        Tnode.invertElev(Tlink.node2(j)) + Tlink.yFull(j));
        Mobjects[LINK]++;
        Mlinks[type]++;
    }

    void  link_validate(int j, Node& Tnode, Link& Tlink, Outfall& Toutfall, Conduit& Tconduit)

    {
        int   n;
        conduit_validate(j, Tlink.subIndex(j), Tnode, Toutfall, Tlink, Tconduit);
   
        n = Tlink.node1(j);
        if ( Tnode.surDepth(n) > 0.0 )
        {
           Tnode.fullDepth(n) = max(Tnode.fullDepth(n),
                                Tlink.offset1(j) + Tlink.yFull(j));
        }

        n = Tlink.node2(j);
        if ( (Tnode.surDepth(n) > 0.0) &&
                Tlink.typee(n) == CONDUIT )
        {
            Tnode.fullDepth(n) = max(Tnode.fullDepth(n),
                                Tlink.offset2(j) + Tlink.yFull(j));
        }
    }


    inline void  conduit_validate(int j, int k,Node& Tnode, Outfall& Toutfall, Link& Tlink, Conduit& Tconduit)

    {
        double aa;
        double lengthFactor, roughness, slope;
        if ( Tlink.offset1(j) < 0.0 )
        {
            Tlink.offset1(j) = 0.0;
        }
        if ( Tlink.offset2(j) < 0.0 )
        {
            Tlink.offset2(j) = 0.0;
        }

        slope = conduit_getSlope(j, Tnode, Tlink, Tconduit);
        Tconduit.slope(k) = slope;

        if (slope < 0.0 &&
            Tlink.typee(j) != DUMMY )
        {
            conduit_reverse(j, k, Tlink, Tconduit);
        }

        roughness = Tconduit.roughness(k);
        lengthFactor = 1.0;
        if (LengtheningStep > 0.0 &&
            Tlink.typee(j) != DUMMY )
        {
            lengthFactor = conduit_getLengthFactor(j, k, roughness, Tlink, Tconduit);
        }

        if ( lengthFactor != 1.0 )
        {
            Tconduit.modLength(k) = lengthFactor *  Tconduit.length(k);
            slope /= lengthFactor;
            roughness = roughness / sqrt(lengthFactor);
        }

        Tconduit.roughFactor(k) = GRAV * pow(roughness,2);
        Tconduit.beta(k) = PHI * sqrt(fabs(slope)) / roughness;
        Tlink.qFull(j) = Tlink.sFull(j) * Tconduit.beta(k);
        Tconduit.qMax(k) = Tlink.sMax(j) * Tconduit.beta(k);
        aa = Tconduit.beta(k) / sqrt(32.2) *
            pow(Tlink.yFull(j), 0.1666667) * 0.3;
        if ( aa >= 1.0 ) Tconduit.superCritical(k) = TRUE;
        else             Tconduit.superCritical(k) = FALSE;
        Tconduit.hasLosses(k) = FALSE;

    }
    

       inline double conduit_getLengthFactor(int j, int k, double roughness, Link& Tlink, Conduit& Tconduit)

    {
        double ratio;
        double yFull;
        double vFull;
        double tStep;

        yFull = Tlink.yFull(j);
        if (xxsect.xsect_isOpen(Tlink.typee(j)) )
        {
            yFull = Tlink.aFull(j) / xxsect.xsect_getWofY(j, yFull, Tlink);
        }
        vFull = 1.0 / roughness * Tlink.sFull(j) *
                sqrt(fabs(Tconduit.slope(k))) / Tlink.aFull(j);

        if ( LengtheningStep == 0.0 ) tStep = RouteStep;
        else                          tStep = min(RouteStep, LengtheningStep);
        ratio = (sqrt(GRAV*yFull) + vFull) * tStep / Tconduit.length(k) ;
        if ( ratio > 1.0 ) return ratio;
        else return 1.0;
    }

     inline double conduit_getSlope(int j, Node& Tnode, Link& Tlink, Conduit& Tconduit)

    {
        double elev1, elev2, delta, slope;
        int k = Tlink.subIndex(j);
        double length = Tconduit.length(k);
        elev1 = Tlink.offset1(j) + Tnode.invertElev(Tlink.node1(j));
        elev2 = Tlink.offset2(j) + Tnode.invertElev(Tlink.node2(j));
        delta = fabs(elev1 - elev2);
        if ( delta < MIN_DELTA_Z )
        {
            delta = MIN_DELTA_Z;
        }
        if ( delta >= length )
        {
            slope = delta / length;
        }
        else slope = delta / sqrt(pow(length,2) - pow(delta,2));
        if ( elev1 < elev2 ) slope = -slope;
        return slope;
    }

    inline void conduit_reverse(int j, int k, Link& Tlink, Conduit& Tconduit)

    {
        int    i;
        double z;
        double cLoss;
        i = Tlink.node1(j);
        Tlink.node1(j) = Tlink.node2(j);
        Tlink.node2(j) = i;
        z = Tlink.offset1(j);
        Tlink.offset1(j) = Tlink.offset2(j);
        Tlink.offset2(j) = z;
        cLoss = Tlink.cLossInlet(j);
        Tlink.cLossInlet(j) = Tlink.cLossOutlet(j);
        Tlink.cLossOutlet(j) = cLoss;
        Tconduit.slope(k) = -Tconduit.slope(k);
        Tlink.direction(j) *= (signed char)-1;
        Tlink.q0(j) = -Tlink.q0(j);
    }

};


#endif 
