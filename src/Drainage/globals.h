//-----------------------------------------------------------------------------
// Reference: EPA SWMM5 (Storm Water Management Model) Version 5.2
// Source: https://github.com/USEPA/Stormwater-Management-Model/blob/develop/src/solver/globals.h
// Note: This file has been simplified for our use case.
// We have omitted considerations for subcatchments, unit conversions, file outputs, etc.
//-----------------------------------------------------------------------------

#ifndef GLOBALS_H
#define GLOBALS_H
#include "../const.h"
#include "enums.h"
#include "../define.h"

static long
                  Nperiods,                 // Number of reporting periods
                  TotalStepCount,           // Total routing steps used 
                  ReportStepCount,          // Reporting routing steps used
                  NonConvergeCount;         // Number of non-converging steps

static int
                  Nobjects[MAX_OBJ_TYPES],  // Number of each object type
                  Nnodes[MAX_NODE_TYPES],   // Number of each node sub-type
                  Nlinks[MAX_LINK_TYPES],   // Number of each link sub-type
                  UnitSystem,               // Unit system
                  FlowUnits,                // Flow units
                  InfilModel,               // Infiltration method
                  RouteModel,               // Flow routing method
                  ForceMainEqn,             // Flow equation for force mains
                  LinkOffsets,              // Link offset convention
                  SurchargeMethod,          // EXTRAN or SLOT method 
                  AllowPonding = TRUE,             // Allow water to pond at nodes
                  InertDamping,             // Degree of inertial damping
                  NormalFlowLtd = BOTH,     // Normal flow limited
                  SlopeWeighting,           // Use slope weighting
                  Compatibility,            // SWMM 5/3/4 compatibility
                  SkipSteadyState,          // Skip over steady state periods
                  IgnoreRainfall,           // Ignore rainfall/runoff
                  IgnoreRDII,               // Ignore RDII
                  IgnoreSnowmelt,           // Ignore snowmelt
                  IgnoreGwater,             // Ignore groundwater
                  IgnoreRouting,            // Ignore flow routing
                  IgnoreQuality,            // Ignore water quality
                  ErrorCode,                // Error code number
                  Warnings,                 // Number of warning messages
                  WetStep,                  // Runoff wet time step (sec)
                  DryStep,                  // Runoff dry time step (sec)
                  ReportStep,               // Reporting time step (sec)
                  RuleStep,                 // Rule evaluation time step (sec)
                  SweepStart,               // Day of year when sweeping starts
                  SweepEnd,                 // Day of year when sweeping ends
                  MaxTrials = 16,                // Max. trials for DW routing
                  NumThreads,               // Number of parallel threads used
                  NumEvents;                // Number of detailed events

static double
                  RouteStep,                // Routing time step (sec)
                  MinRouteStep = 0.1,             // Minimum variable time step (sec)
                  InitRouteStep = 10,
                  LengtheningStep,          // Time step for lengthening (sec)
                  StartDryDays,             // Antecedent dry days
                  CourantFactor = COURANT_FACTOR,    // Courant time step factor
                  MinSurfArea = 1,              // Minimum nodal surface area
                  MinSlope,                 // Minimum conduit slope
                  RunoffError,              // Runoff continuity error
                  GwaterError,              // Groundwater continuity error
                  FlowError,                // Flow routing error
                  QualError,                // Quality routing error
                  HeadTol = HEAD_TOL,                 // DW routing head tolerance (m); host alias
                  SysFlowTol,               // Tolerance for steady system flow
                  LatFlowTol,              // Tolerance for steady nodal inflow
                  VariableStep = 0.0;

#endif
