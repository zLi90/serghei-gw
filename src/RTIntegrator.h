/* -*- mode: c++; c-default-style: "linux" -*- */
/**
 * @file RTIntegrator.h
 * @brief Time integration and global mass-balance tracking for groundwater
 *        reactive transport.
 *
 * This header defines the RTIntegrator class, which computes global (domain-
 * wide) mass-balance diagnostics for each solute species in the groundwater
 * reactive transport system.  It is called once per transport time step to:
 *
 *   1. Accumulate the total dissolved (liquid-phase) mass of each species.
 *   2. Accumulate the total adsorbed (solid-phase) mass of each species.
 *   3. Compute surface-subsurface mass exchange fluxes (when coupled to
 *      the surface transport module).
 *   4. Sum boundary mass inflows and outflows.
 *   5. Sum source/sink contributions (injection wells, drainage, etc.).
 *   6. Compute reaction-induced mass change rates (before/after snapshot).
 *   7. Reduce all quantities across MPI ranks via a single MPI_Allreduce
 *      call for optimal communication efficiency.
 *
 * Physical units (consistent throughout):
 *   - Mass                : mg  (milligrams)
 *   - Mass flux / rate    : mg/s (milligrams per second)
 *   - Concentration c     : mg/L
 *   - Water content wc    : dimensionless (volumetric fraction)
 *   - Cell volume         : m^3
 *   - Conversion factor   : 1000.0  (mg/L * m^3 = 1000 * mg)
 *
 * @see RTIntegratorSW.h  Surface-water counterpart
 * @see RTStateGW.h       Groundwater reactive-transport state
 * @see GwDomain.h        Groundwater domain geometry
 */
#ifndef _RT_INTEGRATOR_H_
#define _RT_INTEGRATOR_H_

#if SERGHEI_SUBSURFACE_TRANSPORT

#include "RTStateGW.h"
#include "SourceSink.h"
#include "GwDomain.h"
#include "RTBCGW.h"
#include "Indexing.h"
#include <vector>
#include <algorithm>
#include <cmath>

/**
 * @class RTIntegrator
 * @brief Per-species mass-balance integrator for groundwater reactive transport.
 *
 * For each solute species, this class tracks the following global (domain-wide)
 * quantities at every transport time step:
 *
 *   - LiquidMassTot  : total dissolved mass in the groundwater domain [mg]
 *   - SolidMassTot   : total adsorbed mass on the solid phase [mg]
 *   - MassExch       : net surface-subsurface exchange flux [mg/s]
 *   - QMassInBC      : total mass inflow through boundaries [mg/s]
 *   - QMassOutBC     : total mass outflow through boundaries [mg/s]
 *   - QMassInSS      : total mass inflow from sources [mg/s]
 *   - QMassOutSS     : total mass outflow from sinks (transpiration, drainage) [mg/s]
 *   - QMassReaction  : net mass change rate due to chemical reactions [mg/s]
 *
 * All quantities are computed locally on each MPI rank and then reduced
 * globally using a single MPI_Allreduce call with a pre-allocated buffer
 * to minimize communication overhead.
 */
class RTIntegrator
{
    Kokkos::Timer timer;

public:
    int n_mass; ///< Number of solute species (constituents) being transported

    /* ------------------------------------------------------------------ */
    /*  Per-species integral quantities (local to this MPI rank)           */
    /* ------------------------------------------------------------------ */

    // --- Total mass inventories [mg] ---
    realArr LiquidMassTot_spec;      ///< Total dissolved mass per species (local rank) [mg]
    realArr LiquidMassTot_spec_glob; ///< Total dissolved mass per species (global, MPI-reduced) [mg]
    realArr SolidMassTot_spec;       ///< Total adsorbed mass per species (local rank) [mg]
    realArr SolidMassTot_spec_glob;  ///< Total adsorbed mass per species (global, MPI-reduced) [mg]
    realArr TotalMassTot_spec_glob;  ///< Liquid + solid mass per species (global) [mg]

    // --- Surface-subsurface exchange flux [mg/s] ---
    realArr MassExch_spec;            ///< Net exchange flux per species (local rank) [mg/s]
    realArr MassExch_spec_glob;       ///< Net exchange flux per species (global, MPI-reduced) [mg/s]
    realArr ExchangeIntoGW_spec_glob; ///< Exchange flux entering groundwater per species [mg/s]

    // --- Boundary fluxes [mg/s] ---
    realArr QMassInBC_spec;                                 ///< Boundary mass inflow per species (local rank) [mg/s]
    realArr QMassOutBC_spec;                                ///< Boundary mass outflow per species (local rank) [mg/s]
    realArr QMassInBC_spec_glob;                            ///< Boundary mass inflow per species (global) [mg/s]
    realArr QMassOutBC_spec_glob;                           ///< Boundary mass outflow per species (global) [mg/s]
    std::vector<std::vector<real>> BoundaryInBC_spec_glob;  ///< Per-RTBC inflow rate [species][boundary] [mg/s]
    std::vector<std::vector<real>> BoundaryOutBC_spec_glob; ///< Per-RTBC outflow rate [species][boundary] [mg/s]

    // --- Source/sink fluxes [mg/s] ---
    realArr QMassInSS_spec;       ///< Source mass inflow per species (local rank) [mg/s]
    realArr QMassOutSS_spec;      ///< Sink mass outflow per species (local rank) [mg/s]
    realArr QMassInSS_spec_glob;  ///< Source mass inflow per species (global) [mg/s]
    realArr QMassOutSS_spec_glob; ///< Sink mass outflow per species (global) [mg/s]

    // --- Reaction flux [mg/s] ---
    realArr QMassReaction_spec;      ///< Reaction-induced mass change rate per species (local) [mg/s]
    realArr QMassReaction_spec_glob; ///< Reaction-induced mass change rate per species (global) [mg/s]

    // --- Accumulated process masses between time-series outputs [mg] ---
    realArr AccumMassExch_spec;                                    ///< Accumulated SW-GW exchange mass, GW sign convention before negation [mg]
    realArr AccumMassInBC_spec;                                    ///< Accumulated boundary inflow mass [mg]
    realArr AccumMassOutBC_spec;                                   ///< Accumulated boundary outflow mass [mg]
    std::vector<std::vector<real>> AccumMassInBC_byBoundary_spec;  ///< Accumulated per-RTBC inflow mass [species][boundary] [mg]
    std::vector<std::vector<real>> AccumMassOutBC_byBoundary_spec; ///< Accumulated per-RTBC outflow mass [species][boundary] [mg]
    realArr AccumMassInSS_spec;                                    ///< Accumulated source inflow mass [mg]
    realArr AccumMassOutSS_spec;                                   ///< Accumulated sink outflow mass [mg]
    realArr AccumMassReaction_spec;                                ///< Accumulated reaction mass change [mg]
    realArr AccumTransportHalf1ResidualMass_spec;                  ///< Accumulated Step-1 transport closure error [mg]
    realArr AccumTransportHalf3ResidualMass_spec;                  ///< Accumulated Step-3 transport closure error [mg]
    realArr AccumNegativeClippedMass_spec;                         ///< Mass added by clipping negative concentrations to zero [mg]

    // --- Direct mass-balance diagnostics [mg/s] ---
    realArr QMassSorption_spec_glob;                ///< Net liquid-to-solid sorption rate [mg/s]
    realArr StorageChange_spec_glob;                ///< d(LiquidMass)/dt between output samples [mg/s]
    realArr BalanceRHS_spec_glob;                   ///< Liquid-phase in - out + exchange + reaction - sorption [mg/s]
    realArr BalanceResidual_spec_glob;              ///< StorageChange - BalanceRHS [mg/s]
    realArr BalanceResidualRel_spec_glob;           ///< Net-balance relative residual [-]
    realArr BalanceResidualRelThroughput_spec_glob; ///< Residual relative to total process throughput [-]
    realArr TransportHalf1Residual_spec_glob;       ///< Step-1 transport closure residual [mg/s]
    realArr TransportHalf3Residual_spec_glob;       ///< Step-3 transport closure residual [mg/s]
    realArr NegativeClippedMassRate_spec_glob;      ///< Negative-concentration clipping mass rate [mg/s]
    realArr PrevLiquidMass_spec_glob;               ///< Previous output sample liquid mass [mg]
    realArr PrevSolidMass_spec_glob;                ///< Previous output sample solid mass [mg]
    real lastBalanceTime;                           ///< Previous diagnostic time [s]
    bool hasBalanceHistory;                         ///< False until the first integrate() call initializes history

    // --- Temporary storage for reaction flux computation ---
    realArr tempLiquidMassBefore; ///< Liquid mass snapshot before reaction step [mg]
    realArr tempSolidMassBefore;  ///< Solid mass snapshot before reaction step [mg]

    /* ------------------------------------------------------------------ */
    /*  Pointers to external objects                                       */
    /* ------------------------------------------------------------------ */
    SourceSink *ss;            ///< Pointer to the source/sink manager
    RTStateGW *rt;             ///< Pointer to the groundwater reactive-transport state
    GwDomain *gdom;            ///< Pointer to the groundwater domain geometry
    std::vector<RTBCGW> *rtbc; ///< Pointer to the vector of boundary condition objects

    /* ------------------------------------------------------------------ */
    /*  Cell count statistics                                              */
    /* ------------------------------------------------------------------ */
    int ncellsBC;      ///< Number of boundary condition cells (local rank)
    int ncellsBC_glob; ///< Number of boundary condition cells (global, MPI-reduced)
    int ncellsIT;      ///< Number of interior (integration) cells (local rank)
    int ncellsIT_glob; ///< Number of interior cells (global, MPI-reduced)

    /* ------------------------------------------------------------------ */
    /*  Pre-allocated MPI communication buffers                            */
    /*  Promoted to class members to avoid dynamic allocation every        */
    /*  time step.  Each species contributes 8 double values; one extra    */
    /*  integer is appended for the boundary cell count.                   */
    /* ------------------------------------------------------------------ */
    std::vector<real> send_buf; ///< Send buffer for MPI_Allreduce [total_mpi_vars]
    std::vector<real> recv_buf; ///< Receive buffer for MPI_Allreduce [total_mpi_vars]
    int total_mpi_vars;         ///< Total number of variables in the MPI buffer

    bool hasAdsorption(int iSpec) const
    {
        if (rt == nullptr || iSpec < 0 || iSpec >= n_mass)
            return false;
        return rt->Sorption_Type(iSpec) != 0.0;
    }

    bool hasReaction() const
    {
        if (rt == nullptr)
            return false;
        return rt->ReactionModule != 0.0 || rt->Nitrogen_Cycle_Simulation != 0.0;
    }

    /* ================================================================== */
    /*  Memory allocation                                                  */
    /* ================================================================== */

    /**
     * @brief Allocate all Kokkos views and MPI buffers for the given number
     *        of solute species.
     *
     * All views are initialised to zero.  The MPI send/receive buffers are
     * sized to hold 8 doubles per species plus 1 extra for the boundary
     * cell count.
     *
     * @param n_mass_  Number of solute species to track.
     */
    void allocate(int n_mass_)
    {
        n_mass = n_mass_;

        // Total mass inventories [mg]
        LiquidMassTot_spec = realArr("LiquidMassTot_spec", n_mass);
        LiquidMassTot_spec_glob = realArr("LiquidMassTot_spec_glob", n_mass);
        SolidMassTot_spec = realArr("SolidMassTot_spec", n_mass);
        SolidMassTot_spec_glob = realArr("SolidMassTot_spec_glob", n_mass);
        TotalMassTot_spec_glob = realArr("TotalMassTot_spec_glob", n_mass);
        MassExch_spec = realArr("MassExch_spec", n_mass);
        MassExch_spec_glob = realArr("MassExch_spec_glob", n_mass);
        ExchangeIntoGW_spec_glob = realArr("ExchangeIntoGW_spec_glob", n_mass);

        // Boundary fluxes [mg/s]
        QMassInBC_spec = realArr("QMassInBC_spec", n_mass);
        QMassInBC_spec_glob = realArr("QMassInBC_spec_glob", n_mass);
        QMassOutBC_spec = realArr("QMassOutBC_spec", n_mass);
        QMassOutBC_spec_glob = realArr("QMassOutBC_spec_glob", n_mass);

        // Source/sink fluxes [mg/s]
        QMassInSS_spec = realArr("QMassInSS_spec", n_mass);
        QMassInSS_spec_glob = realArr("QMassInSS_spec_glob", n_mass);
        QMassOutSS_spec = realArr("QMassOutSS_spec", n_mass);
        QMassOutSS_spec_glob = realArr("QMassOutSS_spec_glob", n_mass);

        // Reaction flux [mg/s]
        QMassReaction_spec = realArr("QMassReaction_spec", n_mass);
        QMassReaction_spec_glob = realArr("QMassReaction_spec_glob", n_mass);
        AccumMassExch_spec = realArr("AccumMassExch_spec", n_mass);
        AccumMassInBC_spec = realArr("AccumMassInBC_spec", n_mass);
        AccumMassOutBC_spec = realArr("AccumMassOutBC_spec", n_mass);
        AccumMassInSS_spec = realArr("AccumMassInSS_spec", n_mass);
        AccumMassOutSS_spec = realArr("AccumMassOutSS_spec", n_mass);
        AccumMassReaction_spec = realArr("AccumMassReaction_spec", n_mass);
        AccumTransportHalf1ResidualMass_spec = realArr("AccumTransportHalf1ResidualMass_spec", n_mass);
        AccumTransportHalf3ResidualMass_spec = realArr("AccumTransportHalf3ResidualMass_spec", n_mass);
        AccumNegativeClippedMass_spec = realArr("AccumNegativeClippedMass_spec", n_mass);
        QMassSorption_spec_glob = realArr("QMassSorption_spec_glob", n_mass);
        StorageChange_spec_glob = realArr("StorageChange_spec_glob", n_mass);
        BalanceRHS_spec_glob = realArr("BalanceRHS_spec_glob", n_mass);
        BalanceResidual_spec_glob = realArr("BalanceResidual_spec_glob", n_mass);
        BalanceResidualRel_spec_glob = realArr("BalanceResidualRel_spec_glob", n_mass);
        BalanceResidualRelThroughput_spec_glob = realArr("BalanceResidualRelThroughput_spec_glob", n_mass);
        TransportHalf1Residual_spec_glob = realArr("TransportHalf1Residual_spec_glob", n_mass);
        TransportHalf3Residual_spec_glob = realArr("TransportHalf3Residual_spec_glob", n_mass);
        NegativeClippedMassRate_spec_glob = realArr("NegativeClippedMassRate_spec_glob", n_mass);
        PrevLiquidMass_spec_glob = realArr("PrevLiquidMass_spec_glob", n_mass);
        PrevSolidMass_spec_glob = realArr("PrevSolidMass_spec_glob", n_mass);

        // Temporary storage for reaction before/after comparison [mg]
        tempLiquidMassBefore = realArr("tempLiquidMassBefore", n_mass);
        tempSolidMassBefore = realArr("tempSolidMassBefore", n_mass);

        // Zero-initialise all Kokkos views
        Kokkos::deep_copy(LiquidMassTot_spec, 0.0);
        Kokkos::deep_copy(LiquidMassTot_spec_glob, 0.0);
        Kokkos::deep_copy(SolidMassTot_spec, 0.0);
        Kokkos::deep_copy(SolidMassTot_spec_glob, 0.0);
        Kokkos::deep_copy(TotalMassTot_spec_glob, 0.0);
        Kokkos::deep_copy(MassExch_spec, 0.0);
        Kokkos::deep_copy(MassExch_spec_glob, 0.0);
        Kokkos::deep_copy(ExchangeIntoGW_spec_glob, 0.0);
        Kokkos::deep_copy(QMassInBC_spec, 0.0);
        Kokkos::deep_copy(QMassInBC_spec_glob, 0.0);
        Kokkos::deep_copy(QMassOutBC_spec, 0.0);
        Kokkos::deep_copy(QMassOutBC_spec_glob, 0.0);
        Kokkos::deep_copy(QMassInSS_spec, 0.0);
        Kokkos::deep_copy(QMassInSS_spec_glob, 0.0);
        Kokkos::deep_copy(QMassOutSS_spec, 0.0);
        Kokkos::deep_copy(QMassOutSS_spec_glob, 0.0);
        Kokkos::deep_copy(QMassReaction_spec, 0.0);
        Kokkos::deep_copy(QMassReaction_spec_glob, 0.0);
        Kokkos::deep_copy(AccumMassExch_spec, 0.0);
        Kokkos::deep_copy(AccumMassInBC_spec, 0.0);
        Kokkos::deep_copy(AccumMassOutBC_spec, 0.0);
        Kokkos::deep_copy(AccumMassInSS_spec, 0.0);
        Kokkos::deep_copy(AccumMassOutSS_spec, 0.0);
        Kokkos::deep_copy(AccumMassReaction_spec, 0.0);
        Kokkos::deep_copy(AccumTransportHalf1ResidualMass_spec, 0.0);
        Kokkos::deep_copy(AccumTransportHalf3ResidualMass_spec, 0.0);
        Kokkos::deep_copy(AccumNegativeClippedMass_spec, 0.0);
        Kokkos::deep_copy(QMassSorption_spec_glob, 0.0);
        Kokkos::deep_copy(StorageChange_spec_glob, 0.0);
        Kokkos::deep_copy(BalanceRHS_spec_glob, 0.0);
        Kokkos::deep_copy(BalanceResidual_spec_glob, 0.0);
        Kokkos::deep_copy(BalanceResidualRel_spec_glob, 0.0);
        Kokkos::deep_copy(BalanceResidualRelThroughput_spec_glob, 0.0);
        Kokkos::deep_copy(TransportHalf1Residual_spec_glob, 0.0);
        Kokkos::deep_copy(TransportHalf3Residual_spec_glob, 0.0);
        Kokkos::deep_copy(NegativeClippedMassRate_spec_glob, 0.0);
        Kokkos::deep_copy(PrevLiquidMass_spec_glob, 0.0);
        Kokkos::deep_copy(PrevSolidMass_spec_glob, 0.0);
        Kokkos::deep_copy(tempLiquidMassBefore, 0.0);
        Kokkos::deep_copy(tempSolidMassBefore, 0.0);
        lastBalanceTime = 0.0;
        hasBalanceHistory = false;

        // Pre-allocate MPI buffers: 8 double values per species + 1 cell count
        total_mpi_vars = n_mass * 8 + 1;
        send_buf.resize(total_mpi_vars, 0.0);
        recv_buf.resize(total_mpi_vars, 0.0);
        BoundaryInBC_spec_glob.assign(n_mass, std::vector<real>());
        BoundaryOutBC_spec_glob.assign(n_mass, std::vector<real>());
        AccumMassInBC_byBoundary_spec.assign(n_mass, std::vector<real>());
        AccumMassOutBC_byBoundary_spec.assign(n_mass, std::vector<real>());
    }

    /* ================================================================== */
    /*  Initialization                                                     */
    /* ================================================================== */

    /**
     * @brief Initialize the integrator by binding external objects and
     *        allocating storage.
     *
     * @param rt_    Groundwater reactive-transport state.
     * @param gdom_  Groundwater domain geometry.
     * @param rtbc_  Vector of reactive-transport boundary conditions.
     * @param ss_    Source/sink manager.
     */
    void initialize(RTStateGW &rt_, GwDomain &gdom_, std::vector<RTBCGW> &rtbc_, SourceSink &ss_)
    {
        rt = &rt_;
        gdom = &gdom_;
        rtbc = &rtbc_;
        ss = &ss_;
        allocate(rt_.n_mass);
        BoundaryInBC_spec_glob.assign(n_mass, std::vector<real>(rtbc_.size(), 0.0));
        BoundaryOutBC_spec_glob.assign(n_mass, std::vector<real>(rtbc_.size(), 0.0));
        AccumMassInBC_byBoundary_spec.assign(n_mass, std::vector<real>(rtbc_.size(), 0.0));
        AccumMassOutBC_byBoundary_spec.assign(n_mass, std::vector<real>(rtbc_.size(), 0.0));
    }

    void clearAccumulatedMasses()
    {
        Kokkos::deep_copy(AccumMassExch_spec, 0.0);
        Kokkos::deep_copy(AccumMassInBC_spec, 0.0);
        Kokkos::deep_copy(AccumMassOutBC_spec, 0.0);
        Kokkos::deep_copy(AccumMassInSS_spec, 0.0);
        Kokkos::deep_copy(AccumMassOutSS_spec, 0.0);
        Kokkos::deep_copy(AccumMassReaction_spec, 0.0);
        Kokkos::deep_copy(AccumTransportHalf1ResidualMass_spec, 0.0);
        Kokkos::deep_copy(AccumTransportHalf3ResidualMass_spec, 0.0);
        Kokkos::deep_copy(AccumNegativeClippedMass_spec, 0.0);
        for (int iSpec = 0; iSpec < (int)AccumMassInBC_byBoundary_spec.size(); ++iSpec)
        {
            std::fill(AccumMassInBC_byBoundary_spec[iSpec].begin(), AccumMassInBC_byBoundary_spec[iSpec].end(), 0.0);
            std::fill(AccumMassOutBC_byBoundary_spec[iSpec].begin(), AccumMassOutBC_byBoundary_spec[iSpec].end(), 0.0);
        }
    }

    /* ================================================================== */
    /*  Reaction flux computation (two-phase snapshot approach)            */
    /* ================================================================== */

    /**
     * @brief Record liquid and solid inventories before the reaction step.
     *
     * The reaction diagnostic is a total chemical-reaction mass change.  A
     * liquid-only snapshot would incorrectly classify equilibrium sorption
     * or desorption occurring during the reaction solve as chemical reaction.
     */
    void recordLiquidMassBefore(RTStateGW const &rt, GwState const &gw, GwDomain const &gdom)
    {
        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            real liquid_local = 0.0;
            Kokkos::parallel_reduce("RTGW_MassBeforeReact", gdom.nCell, KOKKOS_LAMBDA(int idx, real &tmp) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idx, kk, jj, ii);
                iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
                real Vol = gdom.dx * gdom.dy * gdom.dz(iGlob);
                // Use the same water-content time layer as the RT storage term.
                real wc_store = (gw.wc_new(iGlob) > 0.0) ? gw.wc_new(iGlob) : gw.wc(iGlob, 1);
                tmp += wc_store * rt.c(iSpec, iGlob, 1) * Vol * 1000.0; }, Kokkos::Sum<real>(liquid_local));
            tempLiquidMassBefore(iSpec) = liquid_local;

            real solid_local = 0.0;
            if (rt.Sorption_Type(iSpec) != 0.0)
            {
                Kokkos::parallel_reduce("RTGW_SolidMassBeforeReact", gdom.nCell, KOKKOS_LAMBDA(int idx, real &tmp) {
                    int ii, jj, kk, iGlob;
                    gdom.unpackIndices(idx, kk, jj, ii);
                    iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
                    real Vol = gdom.dx * gdom.dy * gdom.dz(iGlob);
                    tmp += rt.rho_b * rt.c_solid(iSpec, iGlob, 1) * Vol * 1000.0;
                }, Kokkos::Sum<real>(solid_local));
            }
            tempSolidMassBefore(iSpec) = solid_local;
        }
        Kokkos::fence();
    }

    /**
     * @brief Compute the chemical reaction mass change after the reaction step.
     *
     * Reaction is measured on liquid plus solid inventory.  Sorption remains
     * a separate interval-level liquid/solid transfer diagnosed from the
     * change in solid inventory between time-series outputs.
     */
    void finalizeReactionFlux(RTStateGW const &rt, GwState const &gw, GwDomain const &gdom, real dt_react)
    {
        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            real liquid_after = 0.0;
            Kokkos::parallel_reduce("RTGW_MassAfterReact", gdom.nCell, KOKKOS_LAMBDA(int idx, real &tmp) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idx, kk, jj, ii);
                iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
                real Vol = gdom.dx * gdom.dy * gdom.dz(iGlob);
                real wc_store = (gw.wc_new(iGlob) > 0.0) ? gw.wc_new(iGlob) : gw.wc(iGlob, 1);
                tmp += wc_store * rt.c(iSpec, iGlob, 1) * Vol * 1000.0; }, Kokkos::Sum<real>(liquid_after));

            real solid_after = 0.0;
            if (rt.Sorption_Type(iSpec) != 0.0)
            {
                Kokkos::parallel_reduce("RTGW_SolidMassAfterReact", gdom.nCell, KOKKOS_LAMBDA(int idx, real &tmp) {
                    int ii, jj, kk, iGlob;
                    gdom.unpackIndices(idx, kk, jj, ii);
                    iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
                    real Vol = gdom.dx * gdom.dy * gdom.dz(iGlob);
                    tmp += rt.rho_b * rt.c_solid(iSpec, iGlob, 1) * Vol * 1000.0;
                }, Kokkos::Sum<real>(solid_after));
            }

            const real total_before = tempLiquidMassBefore(iSpec) + tempSolidMassBefore(iSpec);
            const real total_after = liquid_after + solid_after;
            const real mass_change = total_after - total_before;

            // Positive means net chemical production; negative means consumption.
            QMassReaction_spec(iSpec) = mass_change / dt_react;
            AccumMassReaction_spec(iSpec) += mass_change;
        }
        Kokkos::fence();
    }

    std::vector<real> computeLocalTotalMasses(RTStateGW const &rt, GwDomain const &gdom,
                                              realArr const &water_content) const
    {
        std::vector<real> masses(n_mass, 0.0);
        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            real local_mass = 0.0;
            const bool include_solid = (rt.Sorption_Type(iSpec) != 0.0);
            Kokkos::parallel_reduce("RTGW_TransportAuditMass", gdom.nCell, KOKKOS_LAMBDA(int idom, real &tmp) {
                    int ii, jj, kk;
                    gdom.unpackIndices(idom, kk, jj, ii);
                    const int iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc
                                    + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;
                    const real volume = gdom.dx * gdom.dy * gdom.dz(iGlob);
                    const real liquid = water_content(iGlob) * rt.c(iSpec, iGlob, 1);
                    const real solid = include_solid ? rt.rho_b * rt.c_solid(iSpec, iGlob, 1) : 0.0;
                    tmp += (liquid + solid) * volume * 1000.0; }, Kokkos::Sum<real>(local_mass));
            masses[iSpec] = local_mass;
        }
        Kokkos::fence();
        return masses;
    }

    void accumulateTransportFluxes(GwState const &gw, RTStateGW const &rt, GwDomain const &gdom,
                                   std::vector<RTBCGW> &rtbc, real dt_part,
                                   std::vector<real> const *mass_before = nullptr,
                                   std::vector<real> const *mass_after = nullptr,
                                   int transport_half = 0)
    {
        if (dt_part <= 0.0)
            return;

        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            real exch_tmp = 0.0;
#if SERGHEI_SURFACE_TRANSPORT
            Kokkos::parallel_reduce("RTGW_AccumExch", gdom.ny * gdom.nx, KOKKOS_LAMBDA(int idx, real &tmp) {
                int ii, jj, iGlob;
                unpackIndicesUniformGrid(idx, gdom.ny, gdom.nx, jj, ii);
                iGlob = jj * gdom.nx + ii;
                tmp += rt.ConQss(iSpec, iGlob) * gdom.dx * gdom.dy * 1000.0; }, Kokkos::Sum<real>(exch_tmp));
#endif
            AccumMassExch_spec(iSpec) += exch_tmp * dt_part;

            real qin_bc = 0.0;
            real qout_bc = 0.0;
            if (iSpec < (int)AccumMassInBC_byBoundary_spec.size() &&
                AccumMassInBC_byBoundary_spec[iSpec].size() != rtbc.size())
            {
                AccumMassInBC_byBoundary_spec[iSpec].assign(rtbc.size(), 0.0);
                AccumMassOutBC_byBoundary_spec[iSpec].assign(rtbc.size(), 0.0);
                BoundaryInBC_spec_glob[iSpec].assign(rtbc.size(), 0.0);
                BoundaryOutBC_spec_glob[iSpec].assign(rtbc.size(), 0.0);
            }
            for (int k = 0; k < rtbc.size(); k++)
            {
                if (rtbc[k].ncellsBC > 0 && iSpec < rtbc[k].spec_mass_in_rate.size())
                {
                    real bc_in = rtbc[k].spec_mass_in_rate[iSpec];
                    real bc_out = rtbc[k].spec_mass_out_rate[iSpec];
                    qin_bc += bc_in;
                    qout_bc += bc_out;
                    AccumMassInBC_byBoundary_spec[iSpec][k] += bc_in * dt_part;
                    AccumMassOutBC_byBoundary_spec[iSpec][k] += bc_out * dt_part;
                }
            }
            AccumMassInBC_spec(iSpec) += qin_bc * dt_part;
            AccumMassOutBC_spec(iSpec) += qout_bc * dt_part;

            real ss_in_tmp = 0.0;
            real ss_out_tmp = 0.0;
            real tscf = rt.tscf(iSpec);

            real ss_out_local = 0.0;
            Kokkos::parallel_reduce("RTGW_AccumSSOut", gdom.nCell, KOKKOS_LAMBDA(int idx, real &tmp_out) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idx, kk, jj, ii);
                iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;

                real transp_rate = gw.transpGW(iGlob);
                real evap_rate = gw.evapGW(iGlob);
                real area = gdom.dx * gdom.dy;
                real C_curr = rt.c(iSpec, iGlob, 0);

                if (transp_rate < 0.0 && tscf > 0.0)
                    tmp_out += (-transp_rate) * area * C_curr * 1000.0 * tscf;
                if (evap_rate < 0.0)
                    tmp_out += (-evap_rate) * area * C_curr * 1000.0 * 0.0; }, Kokkos::Sum<real>(ss_out_local));
            ss_out_tmp += ss_out_local;

            if (ss != nullptr)
            {
                for (size_t k = 0; k < ss->gwss.size(); k++)
                {
                    const GwSS &ss_obj = ss->gwss[k];
                    if (ss_obj.ncellsIT <= 0 || ss_obj.sstype != 3)
                        continue;

                    auto icells_local = ss_obj.icells;
                    real drain_Cd_local = ss_obj.drain_Cd;
                    real local_out = 0.0;

                    Kokkos::parallel_reduce("RTGW_AccumDrain", ss_obj.ncellsIT, KOKKOS_LAMBDA(int idx, real &tmp_drain) {
                        int iGlob = icells_local(idx);
                        real h_val = gw.h(iGlob, 1);
                        if (h_val > 0.0)
                        {
                            int ivg = gw.soilID(iGlob) * NVG;
                            real ks = gw.vgTable(ivg);
                            real Gamma = 4.0 * drain_Cd_local * ks / (gdom.dx * gdom.dy);
                            real water_flux_out = Gamma * h_val * gdom.dx * gdom.dy * gdom.dz(iGlob);
                            tmp_drain += water_flux_out * rt.c(iSpec, iGlob, 0) * 1000.0;
                        } }, Kokkos::Sum<real>(local_out));
                    ss_out_tmp += local_out;
                }

                for (size_t k = 0; k < ss->rtgwss.size(); k++)
                {
                    const RTGwSS &rt_ss = ss->rtgwss[k];
                    if (rt_ss.ncellsIT <= 0)
                        continue;

                    real C_inject = 0.0;
                    if (rt_ss.spec_sstype[iSpec] == 0)
                    {
                        C_inject = rt_ss.spec_ssval_const[iSpec];
                    }
                    else if (rt_ss.spec_sstype[iSpec] == 1 || rt_ss.spec_sstype[iSpec] == 2)
                    {
                        const auto &ts = rt_ss.spec_ts[iSpec];
                        real t = gdom.etime;
                        int ii_ts = 0;
                        while (ii_ts < ts.np - 1 && t >= ts.time(ii_ts + 1))
                            ii_ts++;

                        if (rt_ss.spec_sstype[iSpec] == 2)
                        {
                            C_inject = ts.value(ii_ts);
                        }
                        else
                        {
                            int jj_ts = (ii_ts == ts.np - 1) ? ii_ts : ii_ts + 1;
                            if (ii_ts == jj_ts)
                                C_inject = ts.value(ii_ts);
                            else
                                C_inject = ts.value(ii_ts) + (ts.value(jj_ts) - ts.value(ii_ts)) / (ts.time(jj_ts) - ts.time(ii_ts)) * (t - ts.time(ii_ts));
                        }
                    }

                    if (C_inject > 0.0)
                    {
                        real local_in = 0.0;
                        int rtsstype = rt_ss.rtsstype;
                        auto icells_local = rt_ss.icells;
                        auto col_depth_local = rt_ss.col_depth;

                        Kokkos::parallel_reduce("RTGW_AccumSSIn", rt_ss.ncellsIT, KOKKOS_LAMBDA(int idx, real &val) {
                            int iGlob = icells_local(idx);
                            if (rtsstype == 1)
                            {
                                val += C_inject * gdom.dx * gdom.dy * gdom.dz(iGlob);
                            }
                            else if (rtsstype == 3)
                            {
                                int nx_h = gdom.nxhc;
                                int ny_h = gdom.nyhc;
                                int rem = iGlob % (nx_h * ny_h);
                                int jj_h = rem / nx_h;
                                int ii_h = rem % nx_h;
                                int iGlobSW = jj_h * nx_h + ii_h;
                                real H_col = col_depth_local(iGlobSW);
                                if (H_col > 1e-6)
                                    val += C_inject * gdom.dx * gdom.dy * (gdom.dz(iGlob) / H_col);
                            } }, Kokkos::Sum<real>(local_in));
                        ss_in_tmp += local_in;
                    }
                }
            }

            AccumMassInSS_spec(iSpec) += ss_in_tmp * dt_part;
            AccumMassOutSS_spec(iSpec) += ss_out_tmp * dt_part;

            if (mass_before != nullptr && mass_after != nullptr &&
                iSpec < (int)mass_before->size() && iSpec < (int)mass_after->size())
            {
                const real net_transport_rate = qin_bc - qout_bc + ss_in_tmp - ss_out_tmp - exch_tmp;
                const real residual_mass = ((*mass_after)[iSpec] - (*mass_before)[iSpec]) - net_transport_rate * dt_part;
                if (transport_half == 1)
                    AccumTransportHalf1ResidualMass_spec(iSpec) += residual_mass;
                else if (transport_half == 3)
                    AccumTransportHalf3ResidualMass_spec(iSpec) += residual_mass;
            }
        }
    }

    bool hasStrongDirichletBoundary(std::vector<RTBCGW> const &rtbc, int iSpec) const
    {
        for (size_t k = 0; k < rtbc.size(); ++k)
        {
            if (rtbc[k].ncellsBC <= 0)
                continue;

            int current_bctype = rtbc[k].rtbctype;
            if (iSpec >= 0 && iSpec < (int)rtbc[k].spec_bctype.size())
                current_bctype = rtbc[k].spec_bctype[iSpec];

            if (current_bctype == SUB_RT_BC_Dirichlet_CONST ||
                current_bctype == SUB_RT_BC_Dirichlet_T)
                return true;
        }
        return false;
    }

    bool hasImplicitOpenBoundary(std::vector<RTBCGW> const &rtbc, int iSpec) const
    {
        for (size_t k = 0; k < rtbc.size(); ++k)
        {
            if (rtbc[k].ncellsBC <= 0)
                continue;

            int current_bctype = rtbc[k].rtbctype;
            if (iSpec >= 0 && iSpec < (int)rtbc[k].spec_bctype.size())
                current_bctype = rtbc[k].spec_bctype[iSpec];

            if (current_bctype == SUB_RT_BC_FD)
                return true;
        }
        return false;
    }

    void alignExchangeWithSurface(realArr const &exchange_into_sw_glob)
    {
        /*
         * Strict diagnostic mode: do not overwrite the groundwater exchange
         * with the surface-water diagnostic exchange, and do not distribute any
         * SW-GW mismatch into ordinary groundwater boundary fluxes.  The GW
         * time series reports its raw ConQss integral as ExchangeIntoGW.
         */
        (void)exchange_into_sw_glob;
    }

    /* ================================================================== */
    /*  Main integration (mass-balance computation per time step)          */
    /* ================================================================== */

    /**
     * @brief Compute all mass-balance quantities for the current time step
     *        and perform a single MPI_Allreduce to obtain global values.
     *
     * The integration proceeds as follows for each species:
     *
     *   1. Total liquid mass:  sum over all cells of wc * c * Volume * 1000 [mg]
     *   2. Total solid mass:   sum over all cells of c_solid * Volume * 1000 [mg]
     *   3. Surface-subsurface exchange: sum of ConQss * Area * 1000 [mg/s]
     *      (only if SERGHEI_SURFACE_TRANSPORT is enabled)
     *   4. Boundary inflows/outflows: accumulated from boundary condition objects [mg/s]
     *   5. Source/sink inflows: injection wells, column sources [mg/s]
     *      Source/sink outflows: transpiration, drainage [mg/s]
     *   6. Reaction flux: pre-computed by finalizeReactionFlux() [mg/s]
     *
     * All local values are packed into a contiguous send buffer (8 doubles
     * per species + 1 cell count) and reduced globally via a single
     * MPI_Allreduce call.
     *
     * @param gw    Groundwater state (water content, head, soil properties).
     * @param rt    Groundwater reactive-transport state (concentrations).
     * @param gdom  Groundwater domain (cell geometry).
     * @param rtbc  Vector of reactive-transport boundary conditions.
     */
    void integrate(GwState const &gw, RTStateGW const &rt, GwDomain const &gdom,
                   std::vector<RTBCGW> &rtbc, real balance_time = -1.0)
    {
        int ierr = 0;

        // Zero out the send buffer for this time step
        std::fill(send_buf.begin(), send_buf.end(), 0.0);

        // Count the total number of boundary condition cells on this rank
        ncellsBC = 0;
        for (int k = 0; k < rtbc.size(); k++)
            ncellsBC += rtbc[k].ncellsBC;

        // Pack the boundary cell count at the end of the buffer
        send_buf[n_mass * 8] = (real)ncellsBC;

        // Iterate over each solute species
        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            // ==================================================================
            // 1. Total liquid-phase mass [mg]
            //    m_liquid = sum( wc_new * c * dx * dy * dz * 1000 )
            //    where wc_new is the same water-content time layer used by
            //    the RT storage term, c = concentration [mg/L]
            // ==================================================================
            real liquid_tmp = 0.0;
            Kokkos::parallel_reduce("RTGW_Liquid", gdom.nCell, KOKKOS_LAMBDA(int idx, real &tmp) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idx, kk, jj, ii);
                iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
                real Vol = gdom.dx * gdom.dy * gdom.dz(iGlob);
                real wc_store = (gw.wc_new(iGlob) > 0.0) ? gw.wc_new(iGlob) : gw.wc(iGlob, 1);
                tmp += wc_store * rt.c(iSpec, iGlob, 1) * Vol * 1000.0; }, Kokkos::Sum<real>(liquid_tmp));
            LiquidMassTot_spec(iSpec) = liquid_tmp;
            send_buf[iSpec * 8 + 0] = liquid_tmp;

            // ==================================================================
            // 2. Total solid-phase (adsorbed) mass [mg]
            //    m_solid = sum( c_solid * dx * dy * dz * 1000 )
            // ==================================================================
            real solid_tmp = 0.0;
            if (rt.Sorption_Type(iSpec) != 0.0)
            {
                Kokkos::parallel_reduce("RTGW_Solid", gdom.nCell, KOKKOS_LAMBDA(int idx, real &tmp) {
                    int ii, jj, kk, iGlob;
                    gdom.unpackIndices(idx, kk, jj, ii);
                    iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;
                    real Vol = gdom.dx * gdom.dy * gdom.dz(iGlob);
                    // Multiply by the bulk density rt.rho_b, as c_solid is in mg/kg
                    tmp += rt.rho_b * rt.c_solid(iSpec, iGlob, 1) * Vol * 1000.0; }, Kokkos::Sum<real>(solid_tmp));
            }
            SolidMassTot_spec(iSpec) = solid_tmp;
            send_buf[iSpec * 8 + 1] = solid_tmp;

            // ==================================================================
            // 3. Surface-subsurface mass exchange flux [mg/s]
            //    ConQss is the exchange concentration flux [mg/L/s]
            //    multiplied by cell area and unit conversion factor.
            //    Only computed when coupled to the surface transport module.
            // ==================================================================
            real exch_tmp = 0.0;
#if SERGHEI_SURFACE_TRANSPORT
            Kokkos::parallel_reduce("RTGW_Exch", gdom.ny * gdom.nx, KOKKOS_LAMBDA(int idx, real &tmp) {
                int ii, jj, iGlob;
                unpackIndicesUniformGrid(idx, gdom.ny, gdom.nx, jj, ii);
                iGlob = jj*gdom.nx + ii;
                tmp += rt.ConQss(iSpec, iGlob) * gdom.dx * gdom.dy * 1000.0; }, Kokkos::Sum<real>(exch_tmp));
#endif
            MassExch_spec(iSpec) = exch_tmp;
            send_buf[iSpec * 8 + 2] = AccumMassExch_spec(iSpec);

            // ==================================================================
            // 4. Boundary mass inflow and outflow [mg/s]
            //    Accumulated from all boundary condition objects for this
            //    species.  Only boundaries that have active cells on this
            //    rank contribute.
            // ==================================================================
            real qin_bc = 0.0;
            real qout_bc = 0.0;
            if (iSpec < (int)BoundaryInBC_spec_glob.size() &&
                BoundaryInBC_spec_glob[iSpec].size() != rtbc.size())
            {
                BoundaryInBC_spec_glob[iSpec].assign(rtbc.size(), 0.0);
                BoundaryOutBC_spec_glob[iSpec].assign(rtbc.size(), 0.0);
            }
            for (int k = 0; k < rtbc.size(); k++)
            {
                real bc_in = 0.0;
                real bc_out = 0.0;
                // Only include contributions from boundaries that have cells
                // on this MPI rank
                if (rtbc[k].ncellsBC > 0)
                {
                    if (iSpec < rtbc[k].spec_mass_in_rate.size())
                    {
                        bc_in = rtbc[k].spec_mass_in_rate[iSpec];
                        bc_out = rtbc[k].spec_mass_out_rate[iSpec];
                        qin_bc += bc_in;
                        qout_bc += bc_out;
                    }
                }
                if (iSpec < (int)BoundaryInBC_spec_glob.size() &&
                    k < (int)BoundaryInBC_spec_glob[iSpec].size())
                {
                    BoundaryInBC_spec_glob[iSpec][k] = bc_in;
                    BoundaryOutBC_spec_glob[iSpec][k] = bc_out;
                }
            }
            QMassInBC_spec(iSpec) = qin_bc;
            QMassOutBC_spec(iSpec) = qout_bc;
            send_buf[iSpec * 8 + 3] = AccumMassInBC_spec(iSpec);
            send_buf[iSpec * 8 + 4] = AccumMassOutBC_spec(iSpec);

            // ==================================================================
            // 5. Source/sink mass fluxes [mg/s]
            //    Outflows: transpiration and evaporation (upward water flux
            //              carrying solute out of the groundwater domain).
            //    Inflows:  injection wells and reactive-transport-specific
            //              source terms.
            // ==================================================================
            real ss_in_tmp = 0.0;
            real ss_out_tmp = 0.0;

            // Default Transpiration Stream Concentration Factor (TSCF) [-]
            // Fraction of groundwater concentration taken up by plant roots
            real tscf = rt.tscf(iSpec);

            // ------------------------------------------------------------------
            // 5a. Transpiration and evaporation outflows [mg/s]
            //     Transpiration removes solute proportional to TSCF * c.
            //     Evaporation does NOT transport solute (factor = 0.0).
            // ------------------------------------------------------------------
            real ss_out_local = 0.0;

            Kokkos::parallel_reduce("RTGW_SSOut", gdom.nCell, KOKKOS_LAMBDA(int idx, real &tmp_out) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idx, kk, jj, ii);
                iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;

                real transp_rate = gw.transpGW(iGlob);   // Transpiration water flux [m/s] (negative = upward)
                real evap_rate = gw.evapGW(iGlob);       // Evaporation water flux [m/s] (negative = upward)
                real area = gdom.dx * gdom.dy;            // Cell area [m^2]
                real C_curr = rt.c(iSpec, iGlob, 0);     // Current concentration [mg/L]

                // Transpiration removes solute: flux = |rate| * area * c * TSCF * 1000
                if (transp_rate < 0.0 && tscf > 0.0) {
                    tmp_out += (-transp_rate) * area * C_curr * 1000.0 * tscf;
                }
                // Evaporation does NOT transport solute (pure water loss)
                if (evap_rate < 0.0) {
                    tmp_out += (-evap_rate) * area * C_curr * 1000.0 * 0.0;
                } }, Kokkos::Sum<real>(ss_out_local));
            ss_out_tmp += ss_out_local;

            // ------------------------------------------------------------------
            // 5b. Drainage outflows from general groundwater source/sink objects
            //     (type 3 = drainage).  Solute is carried out with the
            //     drainage water flux.
            // ------------------------------------------------------------------
            if (ss != nullptr)
            {
                for (size_t k = 0; k < ss->gwss.size(); k++)
                {
                    const GwSS &ss_obj = ss->gwss[k];
                    if (ss_obj.ncellsIT <= 0)
                        continue;

                    if (ss_obj.sstype == 3)
                    {
                        auto icells_local = ss_obj.icells;
                        real drain_Cd_local = ss_obj.drain_Cd;
                        real local_out = 0.0;

                        Kokkos::parallel_reduce("RTGW_Drain", ss_obj.ncellsIT, KOKKOS_LAMBDA(int idx, real &tmp_drain) {
                            int iGlob = icells_local(idx);
                            real h_val = gw.h(iGlob, 1);     // Hydraulic head [m]
                            if (h_val > 0.0) {
                                int ivg = gw.soilID(iGlob) * NVG;
                                real ks = gw.vgTable(ivg);     // Saturated hydraulic conductivity [m/s]
                                // Drainage conductance: Gamma = 4 * Cd * Ks / (dx * dy)
                                real Gamma = 4.0 * drain_Cd_local * ks / (gdom.dx * gdom.dy);
                                // Drainage water flux = Gamma * h * dx * dy * dz [m^3/s]
                                real water_flux_out = Gamma * h_val * gdom.dx * gdom.dy * gdom.dz(iGlob);
                                // Solute mass flux = water_flux * concentration * 1000 [mg/s]
                                tmp_drain += water_flux_out * rt.c(iSpec, iGlob, 0) * 1000.0;
                            } }, Kokkos::Sum<real>(local_out));
                        ss_out_tmp += local_out;
                    }
                }

                // ------------------------------------------------------------------
                // 5c. Reactive-transport specific source inflows
                //     These inject solute at prescribed concentrations, either
                //     constant or time-varying (linearly interpolated from a
                //     time series).
                // ------------------------------------------------------------------
                for (size_t k = 0; k < ss->rtgwss.size(); k++)
                {
                    const RTGwSS &rt_ss = ss->rtgwss[k];
                    if (rt_ss.ncellsIT <= 0)
                        continue;

                    // Determine the injection concentration for this species
                    real C_inject = 0.0;
                    if (rt_ss.spec_sstype[iSpec] == 0)
                    {
                        // Type 0: constant concentration
                        C_inject = rt_ss.spec_ssval_const[iSpec];
                    }
                    else if (rt_ss.spec_sstype[iSpec] == 1 || rt_ss.spec_sstype[iSpec] == 2)
                    {
                        // Type 1: linearly interpolated time series
                        // Type 2: step-wise constant time series
                        // Use const reference to avoid deep-copying the TimeSeries object
                        const auto &ts = rt_ss.spec_ts[iSpec];
                        real t = gdom.etime;
                        int ii_ts = 0;
                        while (ii_ts < ts.np - 1 && t >= ts.time(ii_ts + 1))
                            ii_ts++;

                        if (rt_ss.spec_sstype[iSpec] == 2)
                        {
                            // Step-wise constant: use the value at the current interval
                            C_inject = ts.value(ii_ts);
                        }
                        else
                        {
                            // Linear interpolation between time series points
                            int jj_ts = (ii_ts == ts.np - 1) ? ii_ts : ii_ts + 1;
                            if (ii_ts == jj_ts)
                                C_inject = ts.value(ii_ts);
                            else
                                C_inject = ts.value(ii_ts) + (ts.value(jj_ts) - ts.value(ii_ts)) / (ts.time(jj_ts) - ts.time(ii_ts)) * (t - ts.time(ii_ts));
                        }
                    }

                    // Compute the injection mass flux if concentration is positive
                    if (C_inject > 0.0)
                    {
                        real local_in = 0.0;
                        int rtsstype = rt_ss.rtsstype;
                        auto icells_local = rt_ss.icells;
                        auto col_depth_local = rt_ss.col_depth;

                        Kokkos::parallel_reduce("RTGW_SSIn", rt_ss.ncellsIT, KOKKOS_LAMBDA(int idx, real &val) {
                            int iGlob = icells_local(idx);
                            if (rtsstype == 1) {
                                // Type 1: volumetric injection rate [mg/s]
                                // = concentration * cell_volume
                                val += C_inject * gdom.dx * gdom.dy * gdom.dz(iGlob);
                            } else if (rtsstype == 3) {
                                // Type 3: column-based injection
                                // Distribute mass flux across the column based on
                                // the ratio of cell thickness to total column depth
                                int nx_h = gdom.nxhc; int ny_h = gdom.nyhc;
                                int rem = iGlob % (nx_h * ny_h);
                                int jj_h = rem / nx_h;  int ii_h = rem % nx_h;
                                int iGlobSW = jj_h * nx_h + ii_h;
                                real H_col = col_depth_local(iGlobSW);  // Column depth [m]
                                if (H_col > 1e-6) {
                                    val += C_inject * gdom.dx * gdom.dy * (gdom.dz(iGlob) / H_col);
                                }
                            } }, Kokkos::Sum<real>(local_in));
                        ss_in_tmp += local_in;
                    }
                }
            }
            QMassInSS_spec(iSpec) = ss_in_tmp;
            QMassOutSS_spec(iSpec) = ss_out_tmp;
            send_buf[iSpec * 8 + 5] = AccumMassInSS_spec(iSpec);
            send_buf[iSpec * 8 + 6] = AccumMassOutSS_spec(iSpec);

            // ==================================================================
            // 6. Reaction-induced mass change rate [mg/s]
            //    Pre-computed by finalizeReactionFlux().
            // ==================================================================
            send_buf[iSpec * 8 + 7] = hasReaction() ? AccumMassReaction_spec(iSpec) : 0.0;
        }

        // ==================================================================
        // MPI global reduction: single MPI_Allreduce replaces the previous
        // per-variable aggregation, dramatically reducing communication
        // overhead for large numbers of species.
        // ==================================================================
        SERGHEI_MPI_CHECK(MPI_Allreduce(send_buf.data(), recv_buf.data(), total_mpi_vars, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD));

        // Unpack the global results from the receive buffer
        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            LiquidMassTot_spec_glob(iSpec) = recv_buf[iSpec * 8 + 0]; // Total liquid mass [mg]
            SolidMassTot_spec_glob(iSpec) = recv_buf[iSpec * 8 + 1];  // Total solid mass [mg]
            MassExch_spec_glob(iSpec) = recv_buf[iSpec * 8 + 2];      // Accumulated exchange mass [mg]
            QMassInBC_spec_glob(iSpec) = recv_buf[iSpec * 8 + 3];     // Accumulated boundary inflow [mg]
            QMassOutBC_spec_glob(iSpec) = recv_buf[iSpec * 8 + 4];    // Accumulated boundary outflow [mg]
            QMassInSS_spec_glob(iSpec) = recv_buf[iSpec * 8 + 5];     // Accumulated source inflow [mg]
            QMassOutSS_spec_glob(iSpec) = recv_buf[iSpec * 8 + 6];    // Accumulated sink outflow [mg]
            QMassReaction_spec_glob(iSpec) = recv_buf[iSpec * 8 + 7]; // Accumulated reaction change [mg]
        }
        ncellsBC_glob = (int)recv_buf[n_mass * 8];

        int n_rtbc = (int)rtbc.size();
        std::vector<real> bc_in_send(n_mass * n_rtbc, 0.0);
        std::vector<real> bc_out_send(n_mass * n_rtbc, 0.0);
        std::vector<real> bc_in_recv(n_mass * n_rtbc, 0.0);
        std::vector<real> bc_out_recv(n_mass * n_rtbc, 0.0);
        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            for (int k = 0; k < n_rtbc; ++k)
            {
                if (iSpec < (int)AccumMassInBC_byBoundary_spec.size() &&
                    k < (int)AccumMassInBC_byBoundary_spec[iSpec].size())
                {
                    bc_in_send[iSpec * n_rtbc + k] = AccumMassInBC_byBoundary_spec[iSpec][k];
                    bc_out_send[iSpec * n_rtbc + k] = AccumMassOutBC_byBoundary_spec[iSpec][k];
                }
            }
        }
        if (n_rtbc > 0)
        {
            SERGHEI_MPI_CHECK(MPI_Allreduce(bc_in_send.data(), bc_in_recv.data(), n_mass * n_rtbc, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD));
            SERGHEI_MPI_CHECK(MPI_Allreduce(bc_out_send.data(), bc_out_recv.data(), n_mass * n_rtbc, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD));
        }

        std::vector<real> half_residual_recv(3 * n_mass, 0.0);
#if SERGHEI_OUTPUT_MASS_CONSERVATION_METRICS
        std::vector<real> half_residual_send(3 * n_mass, 0.0);
        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            half_residual_send[iSpec] = AccumTransportHalf1ResidualMass_spec(iSpec);
            half_residual_send[n_mass + iSpec] = AccumTransportHalf3ResidualMass_spec(iSpec);
            half_residual_send[2 * n_mass + iSpec] = AccumNegativeClippedMass_spec(iSpec);
        }
        if (n_mass > 0)
        {
            SERGHEI_MPI_CHECK(MPI_Allreduce(half_residual_send.data(), half_residual_recv.data(), 3 * n_mass,
                                            SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD));
        }
#endif

        real current_time = (balance_time >= 0.0) ? balance_time : gdom.etime;
        real dt_balance = current_time - lastBalanceTime;
        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            const real liquid_mass = LiquidMassTot_spec_glob(iSpec);
            const real solid_mass = (rt.Sorption_Type(iSpec) != 0.0)
                                      ? SolidMassTot_spec_glob(iSpec) : 0.0;
            const real total_mass = liquid_mass + solid_mass;
            real storage_change = 0.0;
            bool has_interval = hasBalanceHistory && dt_balance > 0.0;
            real inv_dt = has_interval ? (1.0 / dt_balance) : 0.0;
            real exch_rate = MassExch_spec_glob(iSpec) * inv_dt;
            real qin_bc_rate = QMassInBC_spec_glob(iSpec) * inv_dt;
            real qout_bc_rate = QMassOutBC_spec_glob(iSpec) * inv_dt;
            real qin_ss_rate = QMassInSS_spec_glob(iSpec) * inv_dt;
            real qout_ss_rate = QMassOutSS_spec_glob(iSpec) * inv_dt;
            real reaction_rate = QMassReaction_spec_glob(iSpec) * inv_dt;
            real exchange_into_gw = -exch_rate;
            real sorption_rate = 0.0;

            if (has_interval)
            {
                storage_change = (liquid_mass - PrevLiquidMass_spec_glob(iSpec)) / dt_balance;
                if (rt.Sorption_Type(iSpec) != 0.0)
                    sorption_rate = (solid_mass - PrevSolidMass_spec_glob(iSpec)) / dt_balance;
            }

            if (iSpec < (int)BoundaryInBC_spec_glob.size())
            {
                BoundaryInBC_spec_glob[iSpec].assign(n_rtbc, 0.0);
                BoundaryOutBC_spec_glob[iSpec].assign(n_rtbc, 0.0);
                for (int k = 0; k < n_rtbc; ++k)
                {
                    BoundaryInBC_spec_glob[iSpec][k] = bc_in_recv[iSpec * n_rtbc + k] * inv_dt;
                    BoundaryOutBC_spec_glob[iSpec][k] = bc_out_recv[iSpec * n_rtbc + k] * inv_dt;
                }
            }

            real rhs = qin_bc_rate - qout_bc_rate + qin_ss_rate - qout_ss_rate
                     + exchange_into_gw + reaction_rate - sorption_rate;
            real residual = has_interval ? storage_change - rhs : 0.0;
            real scale = std::max(std::fabs(storage_change), std::fabs(rhs));
            real throughput_scale = std::fabs(qin_bc_rate) + std::fabs(qout_bc_rate)
                                  + std::fabs(qin_ss_rate) + std::fabs(qout_ss_rate)
                                  + std::fabs(exchange_into_gw) + std::fabs(reaction_rate)
                                  + std::fabs(sorption_rate);
            // Avoid magnifying machine-level residuals when the groundwater
            // control volume has no physically meaningful solute throughput.
            constexpr real relative_rate_floor = 1.0e-12;

            TotalMassTot_spec_glob(iSpec) = total_mass;
            SolidMassTot_spec_glob(iSpec) = solid_mass;
            MassExch_spec_glob(iSpec) = exch_rate;
            QMassInBC_spec_glob(iSpec) = qin_bc_rate;
            QMassOutBC_spec_glob(iSpec) = qout_bc_rate;
            QMassInSS_spec_glob(iSpec) = qin_ss_rate;
            QMassOutSS_spec_glob(iSpec) = qout_ss_rate;
            QMassReaction_spec_glob(iSpec) = reaction_rate;
            QMassSorption_spec_glob(iSpec) = sorption_rate;
            ExchangeIntoGW_spec_glob(iSpec) = exchange_into_gw;
            StorageChange_spec_glob(iSpec) = storage_change;
            BalanceRHS_spec_glob(iSpec) =
#if SERGHEI_OUTPUT_MASS_CONSERVATION_METRICS
                rhs;
#else
                0.0;
#endif
            BalanceResidual_spec_glob(iSpec) = residual;
            BalanceResidualRel_spec_glob(iSpec) =
#if SERGHEI_OUTPUT_MASS_CONSERVATION_METRICS
                (has_interval && scale > relative_rate_floor) ? residual / scale : 0.0;
#else
                0.0;
#endif
            BalanceResidualRelThroughput_spec_glob(iSpec) =
                (has_interval && throughput_scale > relative_rate_floor) ? residual / throughput_scale : 0.0;
            TransportHalf1Residual_spec_glob(iSpec) =
#if SERGHEI_OUTPUT_MASS_CONSERVATION_METRICS
                has_interval ? half_residual_recv[iSpec] * inv_dt : 0.0;
#else
                0.0;
#endif
            TransportHalf3Residual_spec_glob(iSpec) =
#if SERGHEI_OUTPUT_MASS_CONSERVATION_METRICS
                has_interval ? half_residual_recv[n_mass + iSpec] * inv_dt : 0.0;
#else
                0.0;
#endif
            NegativeClippedMassRate_spec_glob(iSpec) =
#if SERGHEI_OUTPUT_MASS_CONSERVATION_METRICS
                has_interval ? half_residual_recv[2 * n_mass + iSpec] * inv_dt : 0.0;
#else
                0.0;
#endif
            PrevLiquidMass_spec_glob(iSpec) = liquid_mass;
            PrevSolidMass_spec_glob(iSpec) = solid_mass;
        }
        lastBalanceTime = current_time;
        hasBalanceHistory = true;
        clearAccumulatedMasses();

        // Note: MPI_Barrier removed to avoid unnecessary synchronisation stalls
    }
};
#endif
#endif
