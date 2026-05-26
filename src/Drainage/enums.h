#ifndef ENUMS_H
#define ENUMS_H

//-----------------------------------------------------------------------------
// Reference: EPA SWMM5 (Storm Water Management Model) Version 5.2
// Source: https://github.com/USEPA/Stormwater-Management-Model/blob/develop/src/solver/enums.h
// Note: This file has been simplified for our use case.
// We have omitted considerations for subcatchments, unit conversions, file outputs, etc.
//-----------------------------------------------------------------------------

//-------------------------------------
// Names of major object types
//-------------------------------------
 enum ObjectType {
      NODE,                            // conveyance system node
      LINK,                            // conveyance system link
     };

//-------------------------------------
// Names of Node sub-types
//-------------------------------------

enum NodeType
{
     JUNCTION,
     OUTFALL,
     INLET //zn251020
};

//-------------------------------------
// Names of Link sub-types
//-------------------------------------
// const int MAX_LINK_TYPES 1
 enum LinkType {
      CONDUIT,
      PUMP
      };

 enum PumpCurveKind {
      PUMP3_CURVE = 0
      };

/** Pump link types (aligned with EPA SWMM PumpType). */
enum PumpType {
      TYPE3_PUMP = 2,                  // flow vs. head delivered (PUMP3 curve)
      IDEAL_PUMP = 5                   // outflow = inlet inflow + overflow
      };

//-------------------------------------
// Cross section shape types
//-------------------------------------
 enum XsectType {
      DUMMY,                           // 0      
      CIRCULAR,                        // 1      closed
      FILLED_CIRCULAR,                 // 2      closed
      RECT_CLOSED,                     // 3      closed
      RECT_OPEN,                       // 4       
      TRAPEZOIDAL,                     // 5       
      TRIANGULAR,                      // 6       
      PARABOLIC,                       // 7
      POWERFUNC,                       // 8      
      RECT_TRIANG,                     // 9       
      RECT_ROUND,                      // 10
      MOD_BASKET,                      // 11      
      HORIZ_ELLIPSE,                   // 12     closed
      VERT_ELLIPSE,                    // 13     closed
      ARCH,                            // 14     closed
      EGGSHAPED,                       // 15     closed
      HORSESHOE,                       // 16     closed
      GOTHIC,                          // 17     closed
      CATENARY,                        // 18     closed
      SEMIELLIPTICAL,                  // 19     closed
      BASKETHANDLE,                    // 20     closed
      SEMICIRCULAR,                    // 21     closed
      IRREGULAR,                       // 22
      CUSTOM,                          // 23     closed
      FORCE_MAIN,                      // 24     closed
      STREET_XSECT};                   // 25



//-------------------------------------
// Conduit flow classifications
//-------------------------------------
 enum FlowClassType {
      DRY,                             // dry conduit
      UP_DRY,                          // upstream end is dry
      DN_DRY,                          // downstream end is dry
      SUBCRITICAL,                     // sub-critical flow
      SUPCRITICAL,                     // super-critical flow
      UP_CRITICAL,                     // free-fall at upstream end
      DN_CRITICAL,                     // free-fall at downstream end
      MAX_FLOW_CLASSES,                // number of distinct flow classes
      UP_FULL,                         // upstream end is full
      DN_FULL,                         // downstream end is full
      ALL_FULL};                       // completely full



 enum NormalFlowType {
      SLOPE,                           // based on slope only
      FROUDE,                          // based on Fr only
      BOTH,                            // based on slope & Fr
      NEITHER};

 enum OutfallType {
      FREE_OUTFALL,                    // critical depth outfall condition
      NORMAL_OUTFALL,                  // normal flow depth outfall condition
      FIXED_OUTFALL,                   // fixed depth outfall condition
      TIDAL_OUTFALL,                   // variable tidal stage outfall condition
      TIMESERIES_OUTFALL};             // variable time series outfall depth


#endif 
