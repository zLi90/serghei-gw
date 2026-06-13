#ifndef FUNCS_H
#define FUNCS_H

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <map>
#include <algorithm>
#include "enums.h"
#include "globals.h"
#include "../define.h"
#include "xsect.h"
#include "DrainageState.h"
#include "pump.h"
#include "storage.h"

XXsect xxsect;

class ReadFileFuncs{
public:
    real dt;
    double simTime = 0.0;

    struct HostCurvePoint { real x; real y; };
    struct HostCurve {
        int id;
        std::vector<HostCurvePoint> pts;
    };
    std::vector<HostCurve> hostCurves_;
    std::vector<std::vector<std::pair<real, real>>> hostRiverStages_;

    #define UCHAR(x) (((x) >= 'a' && (x) <= 'z') ? ((x)&~32) : (x))

    inline void resetDrainageCounters()
    {
        for (int i = 0; i < MAX_OBJ_TYPES; ++i) {
            Nobjects[i] = 0;
            Mobjects[i] = 0;
        }
        for (int i = 0; i < MAX_LINK_TYPES; ++i) {
            Nlinks[i] = 0;
            Mlinks[i] = 0;
        }
        for (int i = 0; i < MAX_NODE_TYPES; ++i) {
            Nnodes[i] = 0;
            Mnodes[i] = 0;
        }
        hostCurves_.clear();
        hostRiverStages_.clear();
    }

    inline void project_readInput(std::string in, Node &Tnode, Link &Tlink, Conduit &Tconduit,
        Outfall &Toutfall, Storage &Tstorage, Pump &Tpump, PumpCurves &TpumpCurves,
        RiverStages &TRiver)
    {
        resetDrainageCounters();
        input_countObjects(in);
        createObjects(Tnode, Tlink, Tconduit, Toutfall, Tstorage, Tpump, TpumpCurves, TRiver);
        input_readData(in, Tnode, Tlink, Tconduit, Toutfall, Tstorage, Tpump);
        finalizePumpCurves(Tpump, TpumpCurves);
        finalizeRiverStages(TRiver);
    }

    inline void node_initState(int j, Node& Tnode, Storage& Tstorage)
    {
        // seed depth/volume from InitDepth before routing init.
        Tnode.oldDepth(j) = Tnode.initDepth(j);
        Tnode.newDepth(j) = Tnode.oldDepth(j);
        Tnode.crownElev(j) = Tnode.invertElev(j);

        if (Tnode.typee(j) == STORAGE) {
            Tnode.fullVolume(j) = DrainageStorage::storage_getVolume(
                Tnode.subIndex(j), Tnode.fullDepth(j), Tstorage);
            Tnode.oldVolume(j) = Tnode.newVolume(j) = DrainageStorage::storage_getVolume(
                Tnode.subIndex(j), Tnode.oldDepth(j), Tstorage);
        } else if (Tnode.fullDepth(j) > 0.0 && Tnode.pondedArea(j) > 0.0) {
            Tnode.fullVolume(j) = Tnode.fullDepth(j) * Tnode.pondedArea(j);
            Tnode.oldVolume(j) = Tnode.newVolume(j) =
                Tnode.fullVolume(j) * (Tnode.oldDepth(j) / Tnode.fullDepth(j));
        } else {
            Tnode.fullVolume(j) = 0.0;
            Tnode.oldVolume(j) = 0.0;
            Tnode.newVolume(j) = 0.0;
        }

        Tnode.oldLatFlow(j) = 0.0;
        Tnode.newLatFlow(j) = 0.0;
        Tnode.apiExtInflow(j) = 0.0;
        Tnode.losses(j) = 0.0;
        Tnode.inflow(j) = 0.0;
        Tnode.outflow(j) = 0.0;
        Tnode.overflow(j) = 0.0;
    }

    inline void project_validate(Node& Tnode, Link& Tlink, Conduit& Tconduit, Outfall& Toutfall,
        Storage& Tstorage, Pump& Tpump, PumpCurves& TpumpCurves)
    {
        int i;
        for (i = 0; i < Nobjects[NODE]; i++) {
            if (Tnode.initDepth(i) > Tnode.fullDepth(i) + Tnode.surDepth(i)) {
                std::cerr << "Warning: node " << i
                          << " has initial depth greater than maximum depth.\n";
            }
            if (Tnode.typee(i) == STORAGE &&
                DrainageStorage::storage_getVolume(
                    Tnode.subIndex(i), Tnode.fullDepth(i), Tstorage) < 0.0) {
                std::cerr << "Warning: storage node " << i << " has negative full volume.\n";
            }
            node_initState(i, Tnode, Tstorage);
        }
        for ( i=0; i<Nobjects[LINK]; i++)
            link_validate(i, Tnode, Tlink, Toutfall, Tconduit, Tpump, TpumpCurves);
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
        } else if (line.find("[STORAGE]") != std::string::npos) {
            currentSection = "STORAGE";
        } else if (line.find("[CONDUITS]") != std::string::npos) {
            currentSection = "LINKS";
        } else if (line.find("[PUMPS]") != std::string::npos) {
            currentSection = "PUMPS";
        } else if (line.find("[CURVES]") != std::string::npos) {
            currentSection = "CURVES";
        } else if (line.find("[OUTFALLS]") != std::string::npos) {
            currentSection = "OUTFALLS";
        } else if (line.find("[RIVERSTAGE]") != std::string::npos) {
            currentSection = "RIVERSTAGE";
        } else if (line.find(";;") != std::string::npos || line.empty()) {
            continue;
        } else {
            if (currentSection == "NODES") {
                Nobjects[NODE]++;
            }
            else if (currentSection == "STORAGE") {
                Nobjects[NODE]++;
                Nnodes[STORAGE]++;
            }
            else if (currentSection == "LINKS")
            {
                Nobjects[LINK]++;
                Nlinks[CONDUIT]++;
            }
            else if (currentSection == "PUMPS")
            {
                Nobjects[LINK]++;
                Nlinks[PUMP]++;
            }
            else if (currentSection == "CURVES" || currentSection == "RIVERSTAGE")
            {
                /* counted while reading */
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


    inline void createObjects(Node &Tnode, Link &Tlink, Conduit &Tconduit, Outfall &Toutfall,
        Storage &Tstorage, Pump &Tpump, PumpCurves &TpumpCurves, RiverStages &TRiver)
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
        Tnode.SDinflowObs = realArr("N_SDinflowObs", Nobjects[NODE]);
        Tnode.SDoutflowObs = realArr("N_SDoutflowObs", Nobjects[NODE]);
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
        Tnode.outfallType = intArr("N_outfallType", Nobjects[NODE]);
        Tnode.outfallStageSeries = intArr("N_outfallStageSeries", Nobjects[NODE]);
        Tnode.outfallFixedStage = realArr("N_outfallFixedStage", Nobjects[NODE]);
        Tnode.outfallHasFlapGate = intArr("N_outfallHasFlapGate", Nobjects[NODE]);
        for (int j = 0; j < Nobjects[NODE]; ++j) {
            Tnode.outfallType(j) = FREE_OUTFALL;
            Tnode.outfallStageSeries(j) = -1;
            Tnode.outfallFixedStage(j) = 0.0;
            Tnode.outfallHasFlapGate(j) = 0;
        }
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

        if (Nnodes[STORAGE] > 0) {
            Tstorage.shape = intArr("S_shape", Nnodes[STORAGE]);
            Tstorage.a0 = realArr("S_a0", Nnodes[STORAGE]);
            Tstorage.a1 = realArr("S_a1", Nnodes[STORAGE]);
            Tstorage.a2 = realArr("S_a2", Nnodes[STORAGE]);
        }

        if (Nlinks[PUMP] > 0) {
            Tpump.pumpCurve = intArr("P_pumpCurve", Nlinks[PUMP]);
            Tpump.type = intArr("P_type", Nlinks[PUMP]);
            Tpump.yOn = realArr("P_yOn", Nlinks[PUMP]);
            Tpump.yOff = realArr("P_yOff", Nlinks[PUMP]);
            Tpump.initSetting = realArr("P_initSetting", Nlinks[PUMP]);
            Tpump.xMin = realArr("P_xMin", Nlinks[PUMP]);
            Tpump.xMax = realArr("P_xMax", Nlinks[PUMP]);
        }

        TpumpCurves.nCurves = 0;
        TpumpCurves.nPts = 0;
        TpumpCurves.startIdx = intArr("PC_start", MAX_PUMP_CURVES);
        TpumpCurves.nPtsPerCurve = intArr("PC_npts", MAX_PUMP_CURVES);
        TpumpCurves.x = realArr("PC_x", MAX_PUMP_CURVE_PTS);
        TpumpCurves.y = realArr("PC_y", MAX_PUMP_CURVE_PTS);

        TRiver.nSeries = 0;
        TRiver.nPts = 0;
        TRiver.startIdx = intArr("RS_start", MAX_RIVER_STAGE_SERIES);
        TRiver.nPtsPerSeries = intArr("RS_npts", MAX_RIVER_STAGE_SERIES);
        TRiver.time = realArr("RS_time", MAX_RIVER_STAGE_PTS);
        TRiver.stage = realArr("RS_stage", MAX_RIVER_STAGE_PTS);



        for (int j = 0; j < Nobjects[LINK]; j++)
        {
            Tlink.typee(j)   = -1;
            Tlink.cLossInlet(j)   = 0.0;
            Tlink.cLossOutlet(j)  = 0.0;
            Tlink.cLossAvg(j)     = 0.0;
            Tlink.hasFlapGate(j)  = FALSE;
        }

    }

    inline void input_readData(std::string in, Node &Tnode, Link &Tlink, Conduit &Tconduit,
        Outfall &Toutfall, Storage &Tstorage, Pump &Tpump)
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
        } else if (line.find("[STORAGE]") != std::string::npos) {
            currentSection = "STORAGE";
        } else if (line.find("[PUMPS]") != std::string::npos) {
            currentSection = "PUMPS";
        } else if (line.find("[CURVES]") != std::string::npos) {
            currentSection = "CURVES";
        } else if (line.find("[OUTFALLS]") != std::string::npos) {
            currentSection = "OUTFALLS";
        } else if (line.find("[RIVERSTAGE]") != std::string::npos) {
            currentSection = "RIVERSTAGE";
        } else if (line.find(";;") != std::string::npos || line.empty()) {
            continue; 
        } else {
            if (currentSection == "JUNCTIONS") {
                readNodeData(line, Tnode);
            } else if (currentSection == "INLETS") {
                readInletData(line, Tnode);
            } else if (currentSection == "STORAGE") {
                readStorageData(line, Tnode, Tstorage);
            }
            else if (currentSection == "LINKS")
            {
                readLinkData(CONDUIT, line, Tlink, Tconduit, Tnode);
            }
            else if (currentSection == "PUMPS")
            {
                readPumpData(line, Tlink, Tpump, Tnode);
            }
            else if (currentSection == "CURVES")
            {
                readCurveData(line);
            }
            else if (currentSection == "OUTFALLS")
            {
                readOutfallData(line, Tnode);
            }
            else if (currentSection == "RIVERSTAGE")
            {
                readRiverStageData(line);
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
        Tnode.outfallType(j) = FREE_OUTFALL;
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
        Tnode.SDinflowObs(j) = 0.0;
        Tnode.SDoutflowObs(j) = 0.0;
    
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
        Tnode.SDinflowObs(j) = 0.0;
        Tnode.SDoutflowObs(j) = 0.0;

        Mobjects[NODE]++;
    }

    inline int parseStorageShape(const std::string& word)
    {
        std::string s = word;
        for (auto& c : s) c = static_cast<char>(UCHAR(c));
        if (s == "FUNCTIONAL") return FUNCTIONAL_SHAPE;
        if (s == "CYLINDRICAL") return CYLINDRICAL_SHAPE;
        return -1;
    }

    inline void readStorageData(const std::string& line, Node &Tnode, Storage &Tstorage)
    {
        int j = Mobjects[NODE];
        int k = Mnodes[STORAGE];
        std::istringstream iss(line);
        real elevation, maxDepth, initDepth, p1, p2, p3, surDepth, pondedArea, xcoor, ycoor;
        int id;
        std::string shapeWord;
        iss >> id >> elevation >> maxDepth >> initDepth >> shapeWord
            >> p1 >> p2 >> p3 >> surDepth >> pondedArea >> xcoor >> ycoor;

        const int shape = parseStorageShape(shapeWord);
        if (shape < 0) {
            std::cerr << "Warning: unknown storage shape '" << shapeWord
                      << "' on node " << id << "; defaulting to FUNCTIONAL.\n";
        }

        Tnode.typee(j) = STORAGE;
        Tnode.subIndex(j) = k;
        Tnode.invertElev(j) = elevation;
        Tnode.crownElev(j) = elevation;
        Tnode.fullDepth(j) = maxDepth;
        Tnode.initDepth(j) = initDepth;
        Tnode.surDepth(j) = surDepth;
        Tnode.pondedArea(j) = pondedArea;
        Tnode.xcoor(j) = xcoor;
        Tnode.ycoor(j) = ycoor;
        Tnode.inflow(j) = 0.0;
        Tnode.outflow(j) = 0.0;
        Tnode.newVolume(j) = 0.0;
        Tnode.fullVolume(j) = 0.0;
        Tnode.SDinflow(j) = 0.0;
        Tnode.SDoutflow(j) = 0.0;
        Tnode.SDinflowObs(j) = 0.0;
        Tnode.SDoutflowObs(j) = 0.0;

        Tstorage.shape(k) = (shape >= 0) ? shape : FUNCTIONAL_SHAPE;
        if (Tstorage.shape(k) == FUNCTIONAL_SHAPE) {
            Tstorage.a1(k) = p1;
            Tstorage.a2(k) = p2;
            Tstorage.a0(k) = p3;
        } else {
            const real A = p1 / 2.0;
            const real B = p2 / 2.0;
            Tstorage.a0(k) = PI * A * B;
            Tstorage.a1(k) = 0.0;
            Tstorage.a2(k) = 0.0;
        }

        Mobjects[NODE]++;
        Mnodes[STORAGE]++;
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
        Tlink.offset1(j) = InOffset;
        Tlink.offset2(j) = OutOffset;
        Tlink.q0(j)       = Initflow;
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

    inline HostCurve* findOrCreateHostCurve(int curveId)
    {
        for (auto& c : hostCurves_) {
            if (c.id == curveId) return &c;
        }
        HostCurve c;
        c.id = curveId;
        hostCurves_.push_back(c);
        return &hostCurves_.back();
    }

    inline void readCurveData(const std::string& line)
    {
        std::istringstream iss(line);
        int curveId;
        std::string typeWord;
        real x, y;
        iss >> curveId >> typeWord >> x >> y;
        if (typeWord != "pump3" && typeWord != "PUMP3") {
            std::cerr << "Warning: only pump3 curves supported; skipping line.\n";
            return;
        }
        HostCurve* c = findOrCreateHostCurve(curveId);
        HostCurvePoint pt;
        pt.x = x;
        pt.y = y;
        c->pts.push_back(pt);
    }

    inline bool isIdealPumpCurveToken(const std::string& tok)
    {
        return tok == "*" || tok == "IDEAL" || tok == "ideal" || tok == "Ideal" ||
               tok == "-1";
    }

    inline void readPumpData(const std::string& line, Link &Tlink, Pump &Tpump, Node &Tnode)
    {
        (void)Tnode;
        int j = Mobjects[LINK];
        int k = Mlinks[PUMP];
        std::istringstream iss(line);
        int id, from, to, initOn = 1;
        std::string curveTok;
        real startup = 0.0, shutoff = 0.0;
        iss >> id >> from >> to >> curveTok;
        if (!(iss >> initOn)) initOn = 1;
        iss >> startup >> shutoff;

        Tlink.node1(j) = from;
        Tlink.node2(j) = to;
        Tlink.subIndex(j) = k;
        Tlink.typee(j) = PUMP;
        Tlink.offset1(j) = 0.0;
        Tlink.offset2(j) = 0.0;
        Tlink.q0(j) = 0.0;
        Tlink.qFull(j) = 0.0;
        Tlink.yFull(j) = 0.0;
        Tlink.setting(j) = initOn ? 1.0 : 0.0;
        Tlink.targetSetting(j) = Tlink.setting(j);
        Tlink.hasFlapGate(j) = 0;
        Tlink.direction(j) = 1;
        Tlink.qLimit(j) = 0.0;
        Tlink.dqdh(j) = 0.0;

        if (isIdealPumpCurveToken(curveTok)) {
            Tpump.type(k) = IDEAL_PUMP;
            Tpump.pumpCurve(k) = -1;
        } else {
            Tpump.type(k) = TYPE3_PUMP;
            try {
                Tpump.pumpCurve(k) = std::stoi(curveTok);
            } catch (...) {
                std::cerr << "Warning: invalid pump curve id '" << curveTok
                          << "'; treating as IDEAL pump.\n";
                Tpump.type(k) = IDEAL_PUMP;
                Tpump.pumpCurve(k) = -1;
            }
        }
        Tpump.initSetting(k) = Tlink.setting(j);
        Tpump.yOn(k) = startup;
        Tpump.yOff(k) = shutoff;
        Tpump.xMin(k) = 0.0;
        Tpump.xMax(k) = 0.0;

        Mobjects[LINK]++;
        Mlinks[PUMP]++;
    }

    inline bool parseFlapGateToken(const std::string& tok)
    {
        if (tok == "1" || tok == "YES" || tok == "yes" || tok == "Yes" ||
            tok == "TRUE" || tok == "true" || tok == "True")
            return true;
        return false;
    }

    inline void readOutfallData(const std::string& line, Node &Tnode)
    {
        std::istringstream iss(line);
        int nodeId;
        std::string mode;
        real param = 0.0;
        std::string gated;
        iss >> nodeId >> mode >> param;
        if (nodeId < 0 || nodeId >= Nobjects[NODE]) {
            std::cerr << "Warning: invalid outfall node id " << nodeId << "\n";
            return;
        }
        if (mode == "free" || mode == "FREE") {
            Tnode.outfallType(nodeId) = FREE_OUTFALL;
        } else if (mode == "fixed" || mode == "FIXED") {
            Tnode.outfallType(nodeId) = FIXED_OUTFALL;
            Tnode.outfallFixedStage(nodeId) = param;
        } else if (mode == "timeseries" || mode == "TIMESERIES") {
            Tnode.outfallType(nodeId) = TIMESERIES_OUTFALL;
            Tnode.outfallStageSeries(nodeId) = static_cast<int>(param);
        } else {
            std::cerr << "Warning: unknown outfall mode '" << mode << "'\n";
        }
        if (iss >> gated) {
            Tnode.outfallHasFlapGate(nodeId) = parseFlapGateToken(gated) ? 1 : 0;
        }
    }

    inline void readRiverStageData(const std::string& line)
    {
        std::istringstream iss(line);
        int seriesId;
        real t, stage;
        iss >> seriesId >> t >> stage;
        if (seriesId < 0) return;
        if (static_cast<int>(hostRiverStages_.size()) <= seriesId)
            hostRiverStages_.resize(seriesId + 1);
        hostRiverStages_[seriesId].push_back(std::make_pair(t, stage));
    }

    inline void finalizePumpCurves(Pump &Tpump, PumpCurves &TpumpCurves)
    {
        std::sort(hostCurves_.begin(), hostCurves_.end(),
            [](const HostCurve& a, const HostCurve& b){ return a.id < b.id; });

        std::map<int,int> idToIndex;
        int ptCount = 0;
        int cCount = 0;
        for (const auto& hostCurve : hostCurves_) {
            if (cCount >= MAX_PUMP_CURVES) break;
            idToIndex[hostCurve.id] = cCount;
            TpumpCurves.startIdx(cCount) = ptCount;
            TpumpCurves.nPtsPerCurve(cCount) = static_cast<int>(hostCurve.pts.size());
            for (const auto& pt : hostCurve.pts) {
                if (ptCount >= MAX_PUMP_CURVE_PTS) break;
                TpumpCurves.x(ptCount) = pt.x;
                TpumpCurves.y(ptCount) = pt.y;
                ++ptCount;
            }
            ++cCount;
        }
        TpumpCurves.nCurves = cCount;
        TpumpCurves.nPts = ptCount;

        for (int k = 0; k < Nlinks[PUMP]; ++k) {
            if (Tpump.type(k) == IDEAL_PUMP) {
                Tpump.pumpCurve(k) = -1;
                continue;
            }
            const int fileId = Tpump.pumpCurve(k);
            auto it = idToIndex.find(fileId);
            if (it != idToIndex.end()) {
                Tpump.pumpCurve(k) = it->second;
            } else {
                std::cerr << "Warning: pump curve id " << fileId << " not found.\n";
                Tpump.pumpCurve(k) = -1;
            }
        }
    }

    inline void finalizeRiverStages(RiverStages &TRiver)
    {
        int seriesCount = static_cast<int>(hostRiverStages_.size());
        if (seriesCount > MAX_RIVER_STAGE_SERIES)
            seriesCount = MAX_RIVER_STAGE_SERIES;
        int ptCount = 0;
        for (int s = 0; s < seriesCount; ++s) {
            TRiver.startIdx(s) = ptCount;
            TRiver.nPtsPerSeries(s) = static_cast<int>(hostRiverStages_[s].size());
            for (const auto& pr : hostRiverStages_[s]) {
                if (ptCount >= MAX_RIVER_STAGE_PTS) break;
                TRiver.time(ptCount) = pr.first;
                TRiver.stage(ptCount) = pr.second;
                ++ptCount;
            }
        }
        TRiver.nSeries = seriesCount;
        TRiver.nPts = ptCount;
    }

    void  link_validate(int j, Node& Tnode, Link& Tlink, Outfall& Toutfall, Conduit& Tconduit,
        Pump& Tpump, PumpCurves& TpumpCurves)
    {
        if (Tlink.typee(j) == PUMP) {
            pump_validate(j, Tlink.subIndex(j), Tnode, Tlink, Tpump, TpumpCurves);
            return;
        }
        conduit_validate(j, Tlink.subIndex(j), Tnode, Toutfall, Tlink, Tconduit);
   
        int   n;
        n = Tlink.node1(j);
        if ( Tnode.typee(n) != STORAGE || Tnode.surDepth(n) > 0.0 )
        {
           Tnode.fullDepth(n) = max(Tnode.fullDepth(n),
                                Tlink.offset1(j) + Tlink.yFull(j));
        }

        n = Tlink.node2(j);
        if ( (Tnode.typee(n) != STORAGE || Tnode.surDepth(n) > 0.0) &&
                Tlink.typee(j) == CONDUIT )
        {
            Tnode.fullDepth(n) = max(Tnode.fullDepth(n),
                                Tlink.offset2(j) + Tlink.yFull(j));
        }
    }

    inline void pump_validate(int j, int k, Node& Tnode, Link& Tlink, Pump& Tpump,
        PumpCurves& TpumpCurves)
    {
        (void)Tnode;
        if (Tpump.type(k) == IDEAL_PUMP) {
            Tlink.qFull(j) = 0.0;
            if (Tpump.yOn(k) > 0.0 && Tpump.yOn(k) <= Tpump.yOff(k))
                std::cerr << "Warning: pump startup depth should exceed shutoff depth.\n";
            return;
        }
        const int cid = Tpump.pumpCurve(k);
        if (cid < 0 || cid >= TpumpCurves.nCurves) {
            std::cerr << "Warning: pump link " << j << " has invalid curve index.\n";
            return;
        }
        const int start = TpumpCurves.startIdx(cid);
        const int npt = TpumpCurves.nPtsPerCurve(cid);
        if (npt <= 0) return;
        real qMax = TpumpCurves.y(start);
        Tpump.xMin(k) = TpumpCurves.x(start);
        Tpump.xMax(k) = TpumpCurves.x(start);
        for (int p = 1; p < npt; ++p) {
            Tpump.xMin(k) = std::min(Tpump.xMin(k), TpumpCurves.x(start + p));
            Tpump.xMax(k) = std::max(Tpump.xMax(k), TpumpCurves.x(start + p));
            qMax = std::max(qMax, TpumpCurves.y(start + p));
        }
        Tlink.qFull(j) = qMax;
        if (Tpump.yOn(k) > 0.0 && Tpump.yOn(k) <= Tpump.yOff(k))
            std::cerr << "Warning: pump startup depth should exceed shutoff depth.\n";
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
