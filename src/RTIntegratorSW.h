/* -*- mode: c++; c-default-style: "linux" -*- */
/**
 * @file RTIntegratorSW.h
 * @brief Time integration and global mass-balance tracking for surface water
 *        reactive transport.
 *
 * This header defines the RTIntegratorSW class, which computes global
 * (domain-wide) mass-balance diagnostics for each solute species in the
 * surface water reactive transport system.  It is called once per transport
 * time step to:
 *
 *   1. Accumulate the total dissolved mass of each species in the surface
 *      water layer.
 *   2. Compute surface-subsurface mass exchange fluxes (when coupled to
 *      the subsurface transport module).
 *   3. Compute reaction-induced mass change rates (before/after snapshot).
 *   4. Compute rainfall solute input fluxes.
 *   5. Sum boundary mass inflows and outflows.
 *   6. Sum source/sink contributions (injection points, etc.).
 *   7. Reduce all quantities across MPI ranks via a single MPI_Allreduce
 *      call for optimal communication efficiency.
 *
 * Physical units (consistent throughout):
 *   - Mass                : mg  (milligrams)
 *   - Mass flux / rate    : mg/s (milligrams per second)
 *   - Concentration c     : mg/L
 *   - Water depth h       : m   (meters)
 *   - Cell area           : m^2
 *   - Conversion factor   : 1000.0  (mg/L * m^3 = 1000 * mg)
 *
 * @see RTIntegrator.h    Groundwater counterpart
 * @see RTStateSW.h       Surface water reactive-transport state
 * @see Domain.h          Surface water domain geometry
 */
#ifndef _RT_INTEGRATOR_SW_H_
#define _RT_INTEGRATOR_SW_H_

#if SERGHEI_SURFACE_TRANSPORT

#if SERGHEI_RE_MODEL // Renamed from SERGHEI_SUBSURFACE_MODEL to SERGHEI_RE_MODEL
#include "RTStateGW.h"
#endif
#include "RTStateSW.h"
#include "SourceSink.h"
#include "Domain.h"
#include "State.h"
#include "Indexing.h"
#include "BC.h"
#include <vector>
#include <math.h>
#include <cmath>
#include <iostream>
#include <algorithm> // for std::fill

/**
 * @class RTIntegratorSW
 * @brief Per-species mass-balance integrator for surface water reactive
 *        transport.
 *
 * For each solute species, this class tracks the following global
 * (domain-wide) quantities at every transport time step:
 *
 *   - MassTot        : total dissolved mass in the surface water layer [mg]
 *   - MassExch       : net surface-subsurface exchange flux [mg/s]
 *   - QMassInBC      : total mass inflow through boundaries [mg/s]
 *   - QMassOutBC     : total mass outflow through boundaries [mg/s]
 *   - QMassInSS      : total mass inflow from sources [mg/s]
 *   - QMassOutSS     : total mass outflow from sinks [mg/s]
 *   - QMassReaction  : net mass change rate due to chemical reactions [mg/s]
 *   - QMassRain      : mass input from rainfall [mg/s]
 *
 * All quantities are computed locally on each MPI rank and then reduced
 * globally using a single MPI_Allreduce call with a pre-allocated buffer
 * to minimize communication overhead.
 */
class RTIntegratorSW
{
    Kokkos::Timer timer;

public:
    /**
     * Dynamic-array Kokkos reduction used by the surface-water mass audits.
     * The largest result layout (timeseries audit) is:
     *   [0, n_mass)            total/storage mass
     *   [n_mass, 2*n_mass)     exchange mass flux
     *   [2*n_mass, 3*n_mass)   rainfall mass flux
     *
     * Using a dynamic array keeps the implementation valid for any configured
     * number of species while reducing all species in one domain traversal.
     */
    struct FusedReactionMassFunctor
    {
        using value_type = real[];

        realArr2 h4rtsw;
        realArr3 concentration;
        boolArr nodata;
        int nx;
        int hc;
        int rt_nxhc;
        int rt_hc;
        int nSpecies;
        int value_count;
        real areaConversion;

        FusedReactionMassFunctor(State const &state, RTStateSW const &rtsw,
                                 Domain const &dom)
            : h4rtsw(state.h4rtsw), concentration(rtsw.c), nodata(state.isnodata),
              nx(dom.nx), hc(dom.hc), rt_nxhc(rtsw.rt_nxhc), rt_hc(RTStateSW::rt_hc), nSpecies(rtsw.n_mass),
              value_count(rtsw.n_mass),
              areaConversion(dom.cellArea() * 1000.0)
        {
        }

        KOKKOS_INLINE_FUNCTION
        int cellIndex(const int iGlob) const
        {
            const int i = iGlob % nx;
            const int j = iGlob / nx;
            return (j + rt_hc) * rt_nxhc + rt_hc + i;
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int iGlob, value_type update) const
        {
            const int ii = cellIndex(iGlob);
            const int hi = (hc + iGlob / nx) * (nx + 2 * hc) + hc + iGlob % nx;
            if (!nodata(hi))
            {
                const real hArea = h4rtsw(hi, 0) * areaConversion;
                for (int iSpec = 0; iSpec < nSpecies; ++iSpec)
                    update[iSpec] += hArea * concentration(iSpec, ii, 1);
            }
        }

        KOKKOS_INLINE_FUNCTION
        void init(value_type update) const
        {
            for (int i = 0; i < value_count; ++i)
                update[i] = 0.0;
        }

        KOKKOS_INLINE_FUNCTION
        void join(value_type update, const value_type source) const
        {
            for (int i = 0; i < value_count; ++i)
                update[i] += source[i];
        }
    };

    /**
     * One traversal for all transport-side mass increments.  The reduction
     * layout is [exchange species][dry-loss species][rain species].
     */
    struct FusedAppliedTransportMassFunctor
    {
        using value_type = real[];

        realArr2 exchangeMass;
        realArr2 dryLossMass;
        realArr2 rainMass;
        int nSpecies;
        int value_count;
        int rt_nx;
        int rt_nxhc;
        int rt_hc;

        FusedAppliedTransportMassFunctor(RTStateSW const &rtsw)
            : exchangeMass(rtsw.exchange_mass_step),
              dryLossMass(rtsw.dry_loss_mass_step),
              rainMass(rtsw.rain_mass_step),
              nSpecies(rtsw.n_mass), value_count(3 * rtsw.n_mass),
              rt_nx(rtsw.rt_nx), rt_nxhc(rtsw.rt_nxhc), rt_hc(RTStateSW::rt_hc)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int idom, value_type update) const
        {
            const int i = idom % rt_nx;
            const int j = idom / rt_nx;
            const int rt_idx = (j + rt_hc) * rt_nxhc + i + rt_hc;
            for (int iSpec = 0; iSpec < nSpecies; ++iSpec)
            {
                update[iSpec] += exchangeMass(iSpec, rt_idx);
                update[nSpecies + iSpec] += dryLossMass(iSpec, rt_idx);
                update[2 * nSpecies + iSpec] += rainMass(iSpec, rt_idx);
            }
        }

        KOKKOS_INLINE_FUNCTION
        void init(value_type update) const
        {
            for (int i = 0; i < value_count; ++i)
                update[i] = 0.0;
        }

        KOKKOS_INLINE_FUNCTION
        void join(value_type update, const value_type source) const
        {
            for (int i = 0; i < value_count; ++i)
                update[i] += source[i];
        }
    };

    struct FusedDissolvedMassFunctor
    {
        using value_type = real[];
        realArr dissolvedMass;
        int nSpecies;
        int sourceSpecies;
        int value_count;

        FusedDissolvedMassFunctor(realArr dissolvedMass_, int nSpecies_, int sourceSpecies_)
            : dissolvedMass(dissolvedMass_), nSpecies(nSpecies_), sourceSpecies(sourceSpecies_), value_count(nSpecies_)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int iSpec, value_type update) const
        {
            if (iSpec < sourceSpecies)
                update[iSpec] += dissolvedMass(iSpec);
        }

        KOKKOS_INLINE_FUNCTION
        void init(value_type update) const
        {
            for (int i = 0; i < value_count; ++i)
                update[i] = 0.0;
        }

        KOKKOS_INLINE_FUNCTION
        void join(value_type update, const value_type source) const
        {
            for (int i = 0; i < value_count; ++i)
                update[i] += source[i];
        }
    };

    struct FusedTransportFluxFunctor
    {
        using value_type = real[];

#if SERGHEI_RE_MODEL
        realArr2 exchangeConcentration;
#endif
        realArr rainConcentration;
        realArr rainRate;
        boolArr nodata;
        int nx;
        int hc;
        int rt_nxhc;
        int rt_hc;
        int nSpecies;
        int value_count;
        bool rainActive;
        real areaConversion;

        FusedTransportFluxFunctor(RTStateSW const &rtsw, State const &state,
                                  Domain const &dom, SourceSinkData const &ss)
#if SERGHEI_RE_MODEL
            : exchangeConcentration(rtsw.ConQss),
#else
            :
#endif
              rainConcentration(rtsw.RainCon_arr), rainRate(ss.rainRate),
              nodata(state.isnodata), nx(dom.nx), hc(dom.hc), rt_nxhc(rtsw.rt_nxhc), rt_hc(RTStateSW::rt_hc),
              nSpecies(rtsw.n_mass), value_count(2 * rtsw.n_mass),
              rainActive(dom.isRain), areaConversion(dom.cellArea() * 1000.0)
        {
        }

        KOKKOS_INLINE_FUNCTION
        int cellIndex(const int iGlob) const
        {
            const int i = iGlob % nx;
            const int j = iGlob / nx;
            return (j + rt_hc) * rt_nxhc + rt_hc + i;
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int iGlob, value_type update) const
        {
            const int ii = cellIndex(iGlob);
            const int hi = (hc + iGlob / nx) * (nx + 2 * hc) + hc + iGlob % nx;
            if (nodata(hi))
                return;

#if SERGHEI_RE_MODEL
            for (int iSpec = 0; iSpec < nSpecies; ++iSpec)
                update[iSpec] += exchangeConcentration(iSpec, iGlob) * areaConversion;
#endif

            if (rainActive)
            {
                const real rateArea = rainRate(iGlob) * areaConversion;
                for (int iSpec = 0; iSpec < nSpecies; ++iSpec)
                {
                    const real rainCon = rainConcentration(iSpec);
                    if (rainCon > 0.0)
                        update[nSpecies + iSpec] += rateArea * rainCon;
                }
            }
        }

        KOKKOS_INLINE_FUNCTION
        void init(value_type update) const
        {
            for (int i = 0; i < value_count; ++i)
                update[i] = 0.0;
        }

        KOKKOS_INLINE_FUNCTION
        void join(value_type update, const value_type source) const
        {
            for (int i = 0; i < value_count; ++i)
                update[i] += source[i];
        }
    };

    struct FusedTimeseriesMassFunctor
    {
        using value_type = real[];

        realArr waterDepth;
        realArr3 concentration;
#if SERGHEI_RE_MODEL
        realArr2 exchangeConcentration;
#endif
        realArr rainConcentration;
        realArr rainRate;
        boolArr nodata;
        int nx;
        int hc;
        int rt_nxhc;
        int rt_hc;
        int nSpecies;
        int value_count;
        bool rainActive;
        real areaConversion;

        FusedTimeseriesMassFunctor(RTStateSW const &rtsw, State const &state,
                                   Domain const &dom, SourceSinkData const &ss,
                                   const bool rainActive_)
            : waterDepth(state.h), concentration(rtsw.c),
#if SERGHEI_RE_MODEL
              exchangeConcentration(rtsw.ConQss),
#endif
              rainConcentration(rtsw.RainCon_arr), rainRate(ss.rainRate),
              nodata(state.isnodata), nx(dom.nx), hc(dom.hc), rt_nxhc(rtsw.rt_nxhc), rt_hc(RTStateSW::rt_hc),
              nSpecies(rtsw.n_mass), value_count(3 * rtsw.n_mass),
              rainActive(rainActive_), areaConversion(dom.cellArea() * 1000.0)
        {
        }

        KOKKOS_INLINE_FUNCTION
        int cellIndex(const int iGlob) const
        {
            const int i = iGlob % nx;
            const int j = iGlob / nx;
            return (j + rt_hc) * rt_nxhc + rt_hc + i;
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(const int iGlob, value_type update) const
        {
            const int ii = cellIndex(iGlob);
            const int hi = (hc + iGlob / nx) * (nx + 2 * hc) + hc + iGlob % nx;
            if (nodata(hi))
                return;

            const real area = areaConversion;
            const real hArea = waterDepth(hi) * area;
            for (int iSpec = 0; iSpec < nSpecies; ++iSpec)
            {
                update[iSpec] += hArea * concentration(iSpec, ii, 1);
#if SERGHEI_RE_MODEL
                update[nSpecies + iSpec] += exchangeConcentration(iSpec, iGlob) * area;
#endif
            }

            if (rainActive)
            {
                const real rateArea = rainRate(iGlob) * area;
                for (int iSpec = 0; iSpec < nSpecies; ++iSpec)
                {
                    const real rainCon = rainConcentration(iSpec);
                    if (rainCon > 0.0)
                        update[2 * nSpecies + iSpec] += rateArea * rainCon;
                }
            }
        }

        KOKKOS_INLINE_FUNCTION
        void init(value_type update) const
        {
            for (int i = 0; i < value_count; ++i)
                update[i] = 0.0;
        }

        KOKKOS_INLINE_FUNCTION
        void join(value_type update, const value_type source) const
        {
            for (int i = 0; i < value_count; ++i)
                update[i] += source[i];
        }
    };

public:
    int n_mass; ///< Number of solute species (constituents) being transported
    realArr dissolvedMassScratch; ///< Device-side source/species reduction scratch
    Kokkos::View<real *, Kokkos::HostSpace> dissolvedMassHost; ///< One host snapshot for all sources
    bool dissolvedMassScratchReady = false;

    /* ------------------------------------------------------------------ */
    /*  Per-species integral quantities (local to this MPI rank)           */
    /* ------------------------------------------------------------------ */

    // --- Total mass inventory [mg] ---
    realArr MassTot_spec;      ///< Total dissolved mass per species (local rank) [mg]
    realArr MassTot_spec_glob; ///< Total dissolved mass per species (global, MPI-reduced) [mg]

    // --- Surface-subsurface exchange flux [mg/s] ---
    realArr MassExch_spec;            ///< Net exchange flux per species (local rank) [mg/s]
    realArr MassExch_spec_glob;       ///< Net exchange flux per species (global, MPI-reduced) [mg/s]
    realArr ExchangeIntoSW_spec_glob; ///< Exchange flux entering surface water per species [mg/s]

    // --- Boundary fluxes [mg/s] ---
    realArr QMassInBC_spec;       ///< Boundary mass inflow per species (local rank) [mg/s]
    realArr QMassInBC_spec_glob;  ///< Boundary mass inflow per species (global) [mg/s]
    realArr QMassOutBC_spec;      ///< Boundary mass outflow per species (local rank) [mg/s]
    realArr QMassOutBC_spec_glob; ///< Boundary mass outflow per species (global) [mg/s]

    // --- Source/sink fluxes [mg/s] ---
    realArr QMassInSS_spec;       ///< Source mass inflow per species (local rank) [mg/s]
    realArr QMassInSS_spec_glob;  ///< Source mass inflow per species (global) [mg/s]
    realArr QMassOutSS_spec;      ///< Sink mass outflow per species (local rank) [mg/s]
    realArr QMassOutSS_spec_glob; ///< Sink mass outflow per species (global) [mg/s]

    // --- Reaction flux [mg/s] ---
    realArr QMassReaction_spec;      ///< Reaction-induced mass change rate per species (local) [mg/s]
    realArr QMassReaction_spec_glob; ///< Reaction-induced mass change rate per species (global) [mg/s]

    // --- Rainfall solute input [mg/s] ---
    realArr QMassRain_spec;      ///< Rainfall solute mass input per species (local rank) [mg/s]
    realArr QMassRain_spec_glob; ///< Rainfall solute mass input per species (global) [mg/s]

    // --- Accumulated process masses between time-series outputs [mg] ---
    realArr AccumMassExch_spec;     ///< Accumulated SW-GW exchange mass [mg]
    realArr AccumMassReaction_spec; ///< Accumulated reaction mass change [mg]
    realArr AccumMassRain_spec;     ///< Accumulated rainfall input mass [mg]
    realArr AccumMassInBC_spec;     ///< Accumulated boundary inflow mass [mg]
    realArr AccumMassOutBC_spec;    ///< Accumulated boundary outflow mass [mg]
    realArr AccumBalanceMassInBC_spec;  ///< Boundary inflow actually integrated by transport half-steps [mg]
    realArr AccumBalanceMassOutBC_spec; ///< Boundary outflow actually integrated by transport half-steps [mg]
    realArr AccumMassInSS_spec;     ///< Accumulated source inflow mass [mg]
    realArr AccumMassOutSS_spec;    ///< Accumulated sink outflow mass [mg]
    realArr AccumMassDryLoss_spec;  ///< Liquid mass removed when cells dry [mg]
    realArr QMassDryLoss_spec_glob; ///< Dry-cell liquid mass loss rate [mg/s]

    // --- Direct mass-balance diagnostics [mg/s] ---
    realArr StorageChange_spec_glob;                ///< dM/dt between output samples [mg/s]
    realArr BalanceRHS_spec_glob;                   ///< In - out + exchange + reaction + rain [mg/s]
    realArr BalanceResidual_spec_glob;              ///< StorageChange - BalanceRHS [mg/s]
    realArr BalanceResidualRel_spec_glob;           ///< Net-balance relative residual [-]
    realArr BalanceResidualRelThroughput_spec_glob; ///< Residual relative to total process throughput [-]
    realArr PrevMassTot_spec_glob;                  ///< Previous output sample mass [mg]
    real lastBalanceTime;                           ///< Previous diagnostic time [s]
    bool hasBalanceHistory;                         ///< False until first integrate() call initializes history

    // --- Temporary storage for reaction flux computation ---
    realArr tempMassBefore; ///< Mass snapshot before reaction step [mg]
                            ///< Used to compute the reaction flux as:
                            ///< QMassReaction = (m_after - m_before) / dt_react

    // Reused output storage for the fused dynamic-array reductions.
    realArr fusedAuditScratch;

    /* ------------------------------------------------------------------ */
    /*  Cell count statistics                                              */
    /* ------------------------------------------------------------------ */
    int ncellsBC;      ///< Number of boundary condition cells (local rank)
    int ncellsBC_glob; ///< Number of boundary condition cells (global, MPI-reduced)
    int ncellsSS;      ///< Number of source/sink cells (local rank)
    int ncellsSS_glob; ///< Number of source/sink cells (global, MPI-reduced)

    /* ------------------------------------------------------------------ */
    /*  Pointers to external objects                                       */
    /* ------------------------------------------------------------------ */
    RTStateSW *rtsw;          ///< Pointer to the surface water RT state
    State *state;             ///< Pointer to the surface water hydraulic state
    Domain *dom;              ///< Pointer to the surface water domain geometry
    SourceSinkData const *ss; ///< Pointer to the source/sink data (const)
    ExternalBoundaries *ebc;  ///< Pointer to the external boundary conditions

    /* ------------------------------------------------------------------ */
    /*  Pre-allocated MPI communication buffers                            */
    /*  Promoted to class members to avoid dynamic allocation every        */
    /*  time step.  Each species contributes 8 double values; 2 extra      */
    /*  integers are appended for boundary and source/sink cell counts.    */
    /* ------------------------------------------------------------------ */
    std::vector<real> send_buf; ///< Send buffer for MPI_Allreduce [total_mpi_vars]
    std::vector<real> recv_buf; ///< Receive buffer for MPI_Allreduce [total_mpi_vars]
    int total_mpi_vars;         ///< Total number of variables in the MPI buffer

    bool hasReaction() const
    {
        if (rtsw == nullptr)
            return false;
        return rtsw->Nitrogen_Cycle_Simulation == 1.0 && rtsw->ReactionModule != 0.0;
    }

    bool hasRainInput() const
    {
        if (rtsw == nullptr || dom == nullptr || !dom->isRain)
            return false;
        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            if (rtsw->RainCon_arr(iSpec) > 0.0)
                return true;
        }
        return false;
    }

    bool hasSurfaceSubsurfaceExchange() const
    {
#if SERGHEI_RE_MODEL
        return true;
#else
        return false;
#endif
    }

    /* ================================================================== */
    /*  Memory allocation                                                  */
    /* ================================================================== */

    /**
     * @brief Allocate all Kokkos views and MPI buffers for the given number
     *        of solute species.
     *
     * All views are initialised to zero.  The MPI send/receive buffers are
     * sized to hold 9 doubles per species plus 2 extra for cell counts.
     *
     * @param n_mass_  Number of solute species to track.
     */
    void allocate(int n_mass_)
    {
        n_mass = n_mass_;

        // Total mass inventory [mg]
        MassTot_spec = realArr("MassTot_spec", n_mass);
        MassTot_spec_glob = realArr("MassTot_spec_glob", n_mass);
        MassExch_spec = realArr("MassExch_spec", n_mass);
        MassExch_spec_glob = realArr("MassExch_spec_glob", n_mass);
        ExchangeIntoSW_spec_glob = realArr("ExchangeIntoSW_spec_glob", n_mass);

        // Boundary fluxes [mg/s]
        QMassInBC_spec = realArr("QMassInBC_spec", n_mass);
        QMassOutBC_spec = realArr("QMassOutBC_spec", n_mass);
        QMassInBC_spec_glob = realArr("QMassInBC_spec_glob", n_mass);
        QMassOutBC_spec_glob = realArr("QMassOutBC_spec_glob", n_mass);

        // Source/sink fluxes [mg/s]
        QMassInSS_spec = realArr("QMassInSS_spec", n_mass);
        QMassOutSS_spec = realArr("QMassOutSS_spec", n_mass);
        QMassInSS_spec_glob = realArr("QMassInSS_spec_glob", n_mass);
        QMassOutSS_spec_glob = realArr("QMassOutSS_spec_glob", n_mass);

        // Reaction flux [mg/s]
        QMassReaction_spec = realArr("QMassReaction_spec", n_mass);
        QMassReaction_spec_glob = realArr("QMassReaction_spec_glob", n_mass);

        // Rainfall solute input [mg/s]
        QMassRain_spec = realArr("QMassRain_spec", n_mass);
        QMassRain_spec_glob = realArr("QMassRain_spec_glob", n_mass);
        AccumMassExch_spec = realArr("AccumMassExch_spec", n_mass);
        AccumMassReaction_spec = realArr("AccumMassReaction_spec", n_mass);
        AccumMassRain_spec = realArr("AccumMassRain_spec", n_mass);
        AccumMassInBC_spec = realArr("AccumMassInBC_spec", n_mass);
        AccumMassOutBC_spec = realArr("AccumMassOutBC_spec", n_mass);
        AccumBalanceMassInBC_spec = realArr("AccumBalanceMassInBC_spec", n_mass);
        AccumBalanceMassOutBC_spec = realArr("AccumBalanceMassOutBC_spec", n_mass);
        AccumMassInSS_spec = realArr("AccumMassInSS_spec", n_mass);
        AccumMassOutSS_spec = realArr("AccumMassOutSS_spec", n_mass);
        AccumMassDryLoss_spec = realArr("AccumMassDryLoss_spec", n_mass);
        QMassDryLoss_spec_glob = realArr("QMassDryLoss_spec_glob", n_mass);
        StorageChange_spec_glob = realArr("StorageChange_spec_glob", n_mass);
        BalanceRHS_spec_glob = realArr("BalanceRHS_spec_glob", n_mass);
        BalanceResidual_spec_glob = realArr("BalanceResidual_spec_glob", n_mass);
        BalanceResidualRel_spec_glob = realArr("BalanceResidualRel_spec_glob", n_mass);
        BalanceResidualRelThroughput_spec_glob = realArr("BalanceResidualRelThroughput_spec_glob", n_mass);
        PrevMassTot_spec_glob = realArr("PrevMassTot_spec_glob", n_mass);

        // Temporary storage for reaction before/after comparison [mg]
        tempMassBefore = realArr("tempMassBefore", n_mass_);
        fusedAuditScratch = realArr("RTSW_FusedAuditScratch", 3 * n_mass_);

        // Zero-initialise all Kokkos views
        Kokkos::deep_copy(MassTot_spec, 0.0);
        Kokkos::deep_copy(MassTot_spec_glob, 0.0);
        Kokkos::deep_copy(MassExch_spec, 0.0);
        Kokkos::deep_copy(MassExch_spec_glob, 0.0);
        Kokkos::deep_copy(ExchangeIntoSW_spec_glob, 0.0);
        Kokkos::deep_copy(QMassInBC_spec, 0.0);
        Kokkos::deep_copy(QMassOutBC_spec, 0.0);
        Kokkos::deep_copy(QMassInBC_spec_glob, 0.0);
        Kokkos::deep_copy(QMassOutBC_spec_glob, 0.0);
        Kokkos::deep_copy(QMassInSS_spec, 0.0);
        Kokkos::deep_copy(QMassOutSS_spec, 0.0);
        Kokkos::deep_copy(QMassInSS_spec_glob, 0.0);
        Kokkos::deep_copy(QMassOutSS_spec_glob, 0.0);
        Kokkos::deep_copy(QMassReaction_spec, 0.0);
        Kokkos::deep_copy(QMassReaction_spec_glob, 0.0);
        Kokkos::deep_copy(QMassRain_spec, 0.0);
        Kokkos::deep_copy(QMassRain_spec_glob, 0.0);
        Kokkos::deep_copy(AccumMassExch_spec, 0.0);
        Kokkos::deep_copy(AccumMassReaction_spec, 0.0);
        Kokkos::deep_copy(AccumMassRain_spec, 0.0);
        Kokkos::deep_copy(AccumMassInBC_spec, 0.0);
        Kokkos::deep_copy(AccumMassOutBC_spec, 0.0);
        Kokkos::deep_copy(AccumBalanceMassInBC_spec, 0.0);
        Kokkos::deep_copy(AccumBalanceMassOutBC_spec, 0.0);
        Kokkos::deep_copy(AccumMassInSS_spec, 0.0);
        Kokkos::deep_copy(AccumMassOutSS_spec, 0.0);
        // Dry-cell loss is accumulated over the current output interval,
        // just like the other transport-process mass integrals.
        Kokkos::deep_copy(AccumMassDryLoss_spec, 0.0);
        Kokkos::deep_copy(QMassDryLoss_spec_glob, 0.0);
        Kokkos::deep_copy(StorageChange_spec_glob, 0.0);
        Kokkos::deep_copy(BalanceRHS_spec_glob, 0.0);
        Kokkos::deep_copy(BalanceResidual_spec_glob, 0.0);
        Kokkos::deep_copy(BalanceResidualRel_spec_glob, 0.0);
        Kokkos::deep_copy(BalanceResidualRelThroughput_spec_glob, 0.0);
        Kokkos::deep_copy(PrevMassTot_spec_glob, 0.0);
        Kokkos::deep_copy(tempMassBefore, 0.0);
        Kokkos::deep_copy(fusedAuditScratch, 0.0);
        lastBalanceTime = 0.0;
        hasBalanceHistory = false;

        // Pre-allocate MPI buffers: 8 double values per species + 2 cell counts
        total_mpi_vars = n_mass * 9 + 2;
        send_buf.resize(total_mpi_vars, 0.0);
        recv_buf.resize(total_mpi_vars, 0.0);
    }

    /* ================================================================== */
    /*  Initialization                                                     */
    /* ================================================================== */

    /**
     * @brief Initialize the integrator by binding external objects and
     *        allocating storage.
     *
     * @param rtsw_  Surface water reactive-transport state.
     * @param state_ Surface water hydraulic state.
     * @param dom_   Surface water domain geometry.
     * @param ss_    Source/sink data.
     * @param ebc_   External boundary conditions.
     */
    void initialize(RTStateSW &rtsw_, State &state_, Domain &dom_,
                    SourceSinkData const &ss_, ExternalBoundaries &ebc_)
    {
        rtsw = &rtsw_;
        state = &state_;
        dom = &dom_;
        ss = &ss_;
        ebc = &ebc_;
        allocate(rtsw_.n_mass);

        const size_t nSourceSlots = ss_.rtswss.size() * static_cast<size_t>(n_mass);
        if (nSourceSlots > 0)
        {
            dissolvedMassScratch = realArr("RTSW_DissolvedMassScratch", nSourceSlots);
            Kokkos::deep_copy(dissolvedMassScratch, 0.0);
            dissolvedMassHost = Kokkos::View<real *, Kokkos::HostSpace>(
                "RTSW_DissolvedMassHost", nSourceSlots);
            Kokkos::deep_copy(dissolvedMassHost, 0.0);
            dissolvedMassScratchReady = true;
        }
        else
        {
            dissolvedMassScratch = realArr();
            dissolvedMassHost = Kokkos::View<real *, Kokkos::HostSpace>();
            dissolvedMassScratchReady = false;
        }
    }

    void clearAccumulatedMasses()
    {
        Kokkos::deep_copy(AccumMassExch_spec, 0.0);
        Kokkos::deep_copy(AccumMassReaction_spec, 0.0);
        Kokkos::deep_copy(AccumMassRain_spec, 0.0);
        Kokkos::deep_copy(AccumMassInBC_spec, 0.0);
        Kokkos::deep_copy(AccumMassOutBC_spec, 0.0);
        Kokkos::deep_copy(AccumBalanceMassInBC_spec, 0.0);
        Kokkos::deep_copy(AccumBalanceMassOutBC_spec, 0.0);
        Kokkos::deep_copy(AccumMassInSS_spec, 0.0);
        Kokkos::deep_copy(AccumMassOutSS_spec, 0.0);
        Kokkos::deep_copy(AccumMassDryLoss_spec, 0.0);
    }

    void closeExchangeWithStorageResidual()
    {
        /*
         * Strict diagnostic mode: do not force the time-series mass budget to
         * close by assigning the residual to the SW-GW exchange term.  The
         * exchange reported in ExchangeIntoSW_spec_glob is the raw ConQss
         * integral accumulated in integrate()/accumulateTransportFluxes().
         */
    }

    /**
     * Accumulate boundary mass for the conservation audit only.
     *
     * applyrtbc() is executed once per surface transport half-step and its
     * spec_mass_*_rate values correspond to the mass update performed in that
     * half-step.  Keep this integral separate from AccumMass*BC, whose output
     * is the step-end discharge diagnostic used by outlet hydrographs.
     */
    void accumulateBalanceBoundaryMass(ExternalBoundaries const &ebc, real dt_part)
    {
        if (dt_part <= 0.0)
            return;

        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            real qin_bc = 0.0;
            real qout_bc = 0.0;
            for (int k = 0; k < (int)ebc.rtbc.size(); ++k)
            {
                if (iSpec < (int)ebc.rtbc[k].spec_mass_in_rate.size())
                {
                    qin_bc += ebc.rtbc[k].spec_mass_in_rate[iSpec];
                    qout_bc += ebc.rtbc[k].spec_mass_out_rate[iSpec];
                }
            }
            AccumBalanceMassInBC_spec(iSpec) += qin_bc * dt_part;
            AccumBalanceMassOutBC_spec(iSpec) += qout_bc * dt_part;
        }
    }

    /**
     * Accumulate the exact SW-side exchange, rainfall and dry-cell masses
     * applied by one transport half-step. The values are produced in the
     * transport kernel, so the audit uses the same wet/dry decisions and mass
     * increments as the solver.
     */
    void accumulateAppliedTransportMass(RTStateSW const &rtsw, Domain const &dom)
    {
        FusedAppliedTransportMassFunctor functor(rtsw);
        auto result = Kokkos::subview(fusedAuditScratch,
                                      std::make_pair(0, 3 * n_mass));
        Kokkos::parallel_reduce("RTSW_AccumulateAppliedMass_AllSpecies",
                                dom.nCell, functor, result);
        Kokkos::fence("RTSW_AccumulateAppliedMass_AllSpecies_done");

        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            AccumMassExch_spec(iSpec) += fusedAuditScratch(iSpec);
            AccumMassDryLoss_spec(iSpec) += fusedAuditScratch(n_mass + iSpec);
            AccumMassRain_spec(iSpec) += fusedAuditScratch(2 * n_mass + iSpec);
        }
    }

    /* ================================================================== */
    /*  Reaction flux computation (two-phase snapshot approach)            */
    /* ================================================================== */

    /**
     * @brief Record the total mass of each species BEFORE the chemical
     *        reaction step.
     *
     * This method computes the domain-wide dissolved mass using the current
     * water depth and concentration fields.  The result is stored in
     * tempMassBefore for later comparison with the post-reaction mass
     * in finalizeReactionFlux().
     *
     * Mass computation per cell:
     *   m = water_depth [m] * concentration [mg/L] * cell_area [m^2] * 1000
     *
     * @param rtsw   Surface water reactive-transport state (concentrations).
     * @param state  Surface water hydraulic state (water depth).
     * @param dom    Surface water domain (cell geometry).
     */
    void recordMassBefore(RTStateSW const &rtsw, State const &state, Domain const &dom)
    {
        FusedReactionMassFunctor functor(state, rtsw, dom);
        auto reactionResult = Kokkos::subview(fusedAuditScratch,
                                              std::make_pair(0, n_mass));
        Kokkos::parallel_reduce("RT_MassBeforeReact_AllSpecies", dom.nCell,
                                functor, reactionResult);
        Kokkos::fence("RT_MassBeforeReact_AllSpecies_done");
        // Keep the complete before snapshot, but transfer it as one view
        // operation instead of issuing one host-side element access per
        // species.  The explicit fence above is the CPU-consumption barrier.
        Kokkos::deep_copy(tempMassBefore, reactionResult);
    }

    /**
     * @brief Compute the reaction-induced mass change rate for each species
     *        AFTER the chemical reaction step.
     *
     * The reaction flux is computed as:
     *   QMassReaction = (m_after - m_before) / dt_react
     *
     * where m_before was recorded by recordMassBefore() and m_after is the
     * current post-reaction mass.  A positive value indicates net production;
     * a negative value indicates net consumption.
     *
     * @param rtsw      Surface water reactive-transport state (post-reaction concentrations).
     * @param state     Surface water hydraulic state (water depth).
     * @param dom       Surface water domain (cell geometry).
     * @param dt_react  Duration of the reaction time step [s].
     */
    void finalizeReactionFlux(RTStateSW const &rtsw, State const &state, Domain const &dom, real dt_react)
    {
        FusedReactionMassFunctor functor(state, rtsw, dom);
        auto reactionResult = Kokkos::subview(fusedAuditScratch,
                                              std::make_pair(0, n_mass));
        Kokkos::parallel_reduce("RT_MassAfterReact_AllSpecies", dom.nCell,
                                functor, reactionResult);
        Kokkos::fence("RT_MassAfterReact_AllSpecies_done");

        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            const real m_after = fusedAuditScratch(iSpec);

            // Compute the mass change rate for this rank [mg/s]
            QMassReaction_spec(iSpec) = (m_after - tempMassBefore(iSpec)) / dt_react;
            AccumMassReaction_spec(iSpec) += (m_after - tempMassBefore(iSpec));
        }
    }

    void accumulateTransportFluxes(RTStateSW const &rtsw, State const &state, Domain const &dom,
                                   SourceSinkData const &ss, ExternalBoundaries &ebc, real dt_part)
    {
        if (dt_part <= 0.0)
            return;

        // Consume the actual mass dissolved during both Strang transport
        // half-steps.  The surface inventory contains liquid mass only, so the
        // matching source term must be dissolution into the liquid phase rather
        // than external fertilizer added to the undissolved solid pool.
        // Consume and reset every source accumulator on the device first.
        // All source/species results are then copied to Host in one transfer.
        if (dissolvedMassScratchReady)
        {
            for (size_t k = 0; k < ss.rtswss.size(); ++k)
            {
                RTSwSS const &sw_ss_obj = ss.rtswss[k];
                auto result = Kokkos::subview(
                    dissolvedMassScratch,
                    std::make_pair(k * static_cast<size_t>(n_mass),
                                   (k + 1) * static_cast<size_t>(n_mass)));
                if (sw_ss_obj.dissolved_mass_accum.extent(0) == 0)
                {
                    Kokkos::deep_copy(result, 0.0);
                    continue;
                }
                FusedDissolvedMassFunctor functor(
                    sw_ss_obj.dissolved_mass_accum, n_mass,
                    std::min(n_mass, static_cast<int>(sw_ss_obj.dissolved_mass_accum.extent(0))));
                Kokkos::parallel_reduce("RTSW_DissolvedMassReduce", n_mass,
                                        functor, result);
                Kokkos::fence("RTSW_DissolvedMassReduce_done");
                Kokkos::deep_copy(sw_ss_obj.dissolved_mass_accum, 0.0);
            }
            Kokkos::deep_copy(dissolvedMassHost, dissolvedMassScratch);
        }

        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            real qin_bc = 0.0;
            real qout_bc = 0.0;
            for (int k = 0; k < ebc.rtbc.size(); k++)
            {
                if (iSpec < ebc.rtbc[k].spec_mass_in_rate.size())
                {
                    qin_bc += ebc.rtbc[k].spec_mass_in_rate[iSpec];
                    qout_bc += ebc.rtbc[k].spec_mass_out_rate[iSpec];
                }
            }
            AccumMassInBC_spec(iSpec) += qin_bc * dt_part;
            AccumMassOutBC_spec(iSpec) += qout_bc * dt_part;

            real dissolved_mass = 0.0;
            real qout_tmp = 0.0;
            if (dissolvedMassScratchReady)
            {
                for (size_t k = 0; k < ss.rtswss.size(); ++k)
                    dissolved_mass += dissolvedMassHost(k * static_cast<size_t>(n_mass) + iSpec);
            }

            QMassInSS_spec(iSpec) = dissolved_mass / dt_part;
            QMassOutSS_spec(iSpec) = qout_tmp;
            AccumMassInSS_spec(iSpec) += dissolved_mass;
            AccumMassOutSS_spec(iSpec) += qout_tmp * dt_part;
        }
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
     *   1. Total mass:    sum over all cells of h * c * Area * 1000 [mg]
     *   2. Exchange flux: sum of ConQss * Area * 1000 [mg/s]
     *      (only if SERGHEI_RE_MODEL is enabled)
     *   3. Reaction flux: pre-computed by finalizeReactionFlux() [mg/s]
     *   4. Rainfall flux: sum of rainRate * rainConcentration * Area * 1000 [mg/s]
     *   5. Boundary fluxes: accumulated from boundary condition objects [mg/s]
     *   6. Source/sink fluxes: accumulated from source/sink objects [mg/s]
     *
     * All local values are packed into a contiguous send buffer (9 doubles
     * per species + 2 cell counts) and reduced globally via a single
     * MPI_Allreduce call.
     *
     * @param rtsw   Surface water reactive-transport state (concentrations).
     * @param state  Surface water hydraulic state (water depth, rain status).
     * @param dom    Surface water domain (cell geometry, elapsed time).
     * @param ss     Source/sink data.
     * @param ebc    External boundary conditions.
     */
    void integrate(RTStateSW const &rtsw, State const &state, Domain const &dom,
                   SourceSinkData const &ss, ExternalBoundaries &ebc)
    {
        const bool rain_input_active = hasRainInput();

        // Zero out the send buffer for this time step
        std::fill(send_buf.begin(), send_buf.end(), 0.0);

        // ==================================================================
        // Phase 1: one fused Kokkos reduction for all species and spatial
        // integrals (inventory, exchange and rainfall).
        // ==================================================================
        FusedTimeseriesMassFunctor functor(rtsw, state, dom, ss, rain_input_active);
        Kokkos::parallel_reduce("RT_MassAudit_AllSpecies", dom.nCell,
                                functor, fusedAuditScratch);
        Kokkos::fence("RT_MassAudit_AllSpecies_done");

        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            // ==================================================================
            // 1. Total dissolved mass [mg]
            //    m_tot = sum( h * c * cellArea * 1000 )
            //    h = water depth [m], c = concentration [mg/L]
            // ==================================================================
            const real m_tot = fusedAuditScratch(iSpec);
            MassTot_spec(iSpec) = m_tot;
            send_buf[iSpec * 9 + 0] = m_tot;

            // ==================================================================
            // 2. Surface-subsurface exchange flux [mg/s]
            //    ConQss is the concentration exchange rate [mg/L/s]
            //    multiplied by cell area and unit conversion.
            //    Only computed when coupled to the subsurface module.
            // ==================================================================
#if SERGHEI_RE_MODEL // Renamed from SERGHEI_SUBSURFACE_MODEL to SERGHEI_RE_MODEL
            const real m_exch = fusedAuditScratch(n_mass + iSpec);
            MassExch_spec(iSpec) = m_exch;
            send_buf[iSpec * 9 + 1] = AccumMassExch_spec(iSpec);
#endif

            // ==================================================================
            // 3. Reaction-induced mass change rate [mg/s]
            //    Pre-computed by finalizeReactionFlux().
            // ==================================================================
            send_buf[iSpec * 9 + 2] = hasReaction() ? AccumMassReaction_spec(iSpec) : 0.0;

            // ==================================================================
            // 4. Rainfall solute input [mg/s]
            //    m_rain = sum( rainRate * rainConcentration * cellArea * 1000 )
            //    Only computed when rainfall is active and the species has
            //    a non-zero rain concentration.
            // ==================================================================
            const real m_rain = fusedAuditScratch(2 * n_mass + iSpec);
            QMassRain_spec(iSpec) = m_rain;
            send_buf[iSpec * 9 + 3] = rain_input_active ? AccumMassRain_spec(iSpec) : 0.0;
        }

        // ==================================================================
        // Phase 2: Boundary and source/sink fluxes (CPU-side accumulation)
        // ==================================================================

        // Count boundary condition cells on this rank
        ncellsBC = 0;
        for (int k = 0; k < ebc.rtbc.size(); k++)
            ncellsBC += ebc.rtbc[k].ncellsBC;

        // Accumulate boundary inflows and outflows for each species
        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            real qin_bc = 0.0, qout_bc = 0.0;

            // Sum contributions from all boundary segments
            for (int k = 0; k < ebc.rtbc.size(); k++)
            {
                if (iSpec < ebc.rtbc[k].spec_mass_in_rate.size())
                {
                    qin_bc += ebc.rtbc[k].spec_mass_in_rate[iSpec];
                    qout_bc += ebc.rtbc[k].spec_mass_out_rate[iSpec];
                }
            }

            QMassInBC_spec(iSpec) = qin_bc;
            QMassOutBC_spec(iSpec) = qout_bc;
            send_buf[iSpec * 9 + 4] = AccumMassInBC_spec(iSpec);
            send_buf[iSpec * 9 + 5] = AccumMassOutBC_spec(iSpec);
        }

        // ==================================================================
        // Phase 3: Source/sink fluxes
        // ==================================================================

        // Count source/sink cells on this rank
        ncellsSS = 0;
        for (int k = 0; k < ss.rtswss.size(); k++)
            ncellsSS += ss.rtswss[k].ncellsIT;

        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            send_buf[iSpec * 9 + 6] = AccumMassInSS_spec(iSpec);
            send_buf[iSpec * 9 + 7] = AccumMassOutSS_spec(iSpec);
            send_buf[iSpec * 9 + 8] = AccumMassDryLoss_spec(iSpec);
        }

        // Pack cell counts at the end of the buffer
        send_buf[n_mass * 9 + 0] = (real)ncellsBC;
        send_buf[n_mass * 9 + 1] = (real)ncellsSS;

        // ==================================================================
        // Single MPI_Allreduce for all species and quantities
        // ==================================================================
        SERGHEI_MPI_CHECK(MPI_Allreduce(send_buf.data(), recv_buf.data(), total_mpi_vars, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD));

        // Boundary rates written to the time series remain the existing
        // step-end diagnostics.  Reduce a separate half-step mass integral
        // for the conservation equation so changing the audit cannot alter
        // outlet discharge curves.
        std::vector<real> balance_bc_send(2 * n_mass, 0.0);
        std::vector<real> balance_bc_recv(2 * n_mass, 0.0);
        for (int i = 0; i < n_mass; ++i)
        {
            balance_bc_send[i] = AccumBalanceMassInBC_spec(i);
            balance_bc_send[n_mass + i] = AccumBalanceMassOutBC_spec(i);
        }
        SERGHEI_MPI_CHECK(MPI_Allreduce(balance_bc_send.data(), balance_bc_recv.data(), 2 * n_mass,
                                        SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD));

        // Unpack the global results from the receive buffer
        for (int i = 0; i < n_mass; ++i)
        {
            MassTot_spec_glob(i) = recv_buf[i * 9 + 0];       // Total mass [mg]
            MassExch_spec_glob(i) = recv_buf[i * 9 + 1];      // Accumulated exchange mass [mg]
            QMassReaction_spec_glob(i) = recv_buf[i * 9 + 2]; // Accumulated reaction change [mg]
            QMassRain_spec_glob(i) = recv_buf[i * 9 + 3];     // Accumulated rainfall mass [mg]
            QMassInBC_spec_glob(i) = recv_buf[i * 9 + 4];     // Accumulated boundary inflow [mg]
            QMassOutBC_spec_glob(i) = recv_buf[i * 9 + 5];    // Accumulated boundary outflow [mg]
            QMassInSS_spec_glob(i) = recv_buf[i * 9 + 6];     // Accumulated source inflow [mg]
            QMassOutSS_spec_glob(i) = recv_buf[i * 9 + 7];    // Accumulated sink outflow [mg]
            QMassDryLoss_spec_glob(i) = recv_buf[i * 9 + 8];  // Accumulated dry-cell loss [mg]
        }
        ncellsBC_glob = (int)recv_buf[n_mass * 9 + 0];
        ncellsSS_glob = (int)recv_buf[n_mass * 9 + 1];

        real current_time = dom.etime;
        real dt_balance = current_time - lastBalanceTime;
        for (int i = 0; i < n_mass; ++i)
        {
            real storage_change = 0.0;
            bool has_interval = hasBalanceHistory && dt_balance > 0.0;
            real inv_dt = has_interval ? (1.0 / dt_balance) : 0.0;
            real exchange_into_sw = MassExch_spec_glob(i) * inv_dt;
            real reaction_rate = QMassReaction_spec_glob(i) * inv_dt;
            real rain_rate = QMassRain_spec_glob(i) * inv_dt;
            real qin_bc_rate = QMassInBC_spec_glob(i) * inv_dt;
            real qout_bc_rate = QMassOutBC_spec_glob(i) * inv_dt;
            real balance_qin_bc_rate = balance_bc_recv[i] * inv_dt;
            real balance_qout_bc_rate = balance_bc_recv[n_mass + i] * inv_dt;
            real qin_ss_rate = QMassInSS_spec_glob(i) * inv_dt;
            real qout_ss_rate = QMassOutSS_spec_glob(i) * inv_dt;
            // recv_buf stores dry-cell loss as an accumulated mass [mg] for
            // the current output interval. Convert it to the same rate unit
            // used by the other budget terms before using or writing it.
            real dry_loss_rate = has_interval
                               ? QMassDryLoss_spec_glob(i) * inv_dt
                               : 0.0;
            if (has_interval)
                storage_change = (MassTot_spec_glob(i) - PrevMassTot_spec_glob(i)) / dt_balance;

            real rhs = balance_qin_bc_rate - balance_qout_bc_rate + qin_ss_rate - qout_ss_rate + rain_rate + exchange_into_sw + reaction_rate - dry_loss_rate;
            real residual = has_interval ? storage_change - rhs : 0.0;
            real scale = std::max(std::fabs(storage_change), std::fabs(rhs));
            real throughput_scale = std::fabs(balance_qin_bc_rate) + std::fabs(balance_qout_bc_rate) + std::fabs(qin_ss_rate) + std::fabs(qout_ss_rate) + std::fabs(rain_rate) + std::fabs(exchange_into_sw) + std::fabs(reaction_rate) + std::fabs(dry_loss_rate);
            // Below this absolute rate, both numerator and denominator are
            // floating-point noise for the mg/s diagnostics. Reporting a
            // relative value would amplify roundoff after complete drying.
            constexpr real relative_rate_floor = 1.0e-12;

            MassExch_spec_glob(i) = exchange_into_sw;
            QMassReaction_spec_glob(i) = reaction_rate;
            QMassRain_spec_glob(i) = rain_rate;
            QMassInBC_spec_glob(i) = qin_bc_rate;
            QMassOutBC_spec_glob(i) = qout_bc_rate;
            QMassInSS_spec_glob(i) = qin_ss_rate;
            QMassOutSS_spec_glob(i) = qout_ss_rate;
            QMassDryLoss_spec_glob(i) = dry_loss_rate;
            ExchangeIntoSW_spec_glob(i) = exchange_into_sw;
            StorageChange_spec_glob(i) = storage_change;
            BalanceRHS_spec_glob(i) =
#if SERGHEI_OUTPUT_MASS_CONSERVATION_METRICS
                rhs;
#else
                0.0;
#endif
            BalanceResidual_spec_glob(i) = residual;
            BalanceResidualRel_spec_glob(i) =
#if SERGHEI_OUTPUT_MASS_CONSERVATION_METRICS
                (has_interval && scale > relative_rate_floor) ? residual / scale : 0.0;
#else
                0.0;
#endif
            BalanceResidualRelThroughput_spec_glob(i) =
                (has_interval && throughput_scale > relative_rate_floor) ? residual / throughput_scale : 0.0;
            PrevMassTot_spec_glob(i) = MassTot_spec_glob(i);
        }
        lastBalanceTime = current_time;
        hasBalanceHistory = true;
        clearAccumulatedMasses();
    }
};

#endif // SERGHEI_SURFACE_TRANSPORT
#endif // _RT_INTEGRATOR_SW_H_
