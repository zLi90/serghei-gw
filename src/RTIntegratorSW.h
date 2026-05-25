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

#if SERGHEI_RE_MODEL // [MERGED] Renamed from SERGHEI_SUBSURFACE_MODEL to SERGHEI_RE_MODEL
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
    int n_mass; ///< Number of solute species (constituents) being transported

    /* ------------------------------------------------------------------ */
    /*  Per-species integral quantities (local to this MPI rank)           */
    /* ------------------------------------------------------------------ */

    // --- Total mass inventory [mg] ---
    realArr MassTot_spec;           ///< Total dissolved mass per species (local rank) [mg]
    realArr MassTot_spec_glob;      ///< Total dissolved mass per species (global, MPI-reduced) [mg]

    // --- Surface-subsurface exchange flux [mg/s] ---
    realArr MassExch_spec;          ///< Net exchange flux per species (local rank) [mg/s]
    realArr MassExch_spec_glob;     ///< Net exchange flux per species (global, MPI-reduced) [mg/s]

    // --- Boundary fluxes [mg/s] ---
    realArr QMassInBC_spec;         ///< Boundary mass inflow per species (local rank) [mg/s]
    realArr QMassInBC_spec_glob;    ///< Boundary mass inflow per species (global) [mg/s]
    realArr QMassOutBC_spec;        ///< Boundary mass outflow per species (local rank) [mg/s]
    realArr QMassOutBC_spec_glob;   ///< Boundary mass outflow per species (global) [mg/s]

    // --- Source/sink fluxes [mg/s] ---
    realArr QMassInSS_spec;         ///< Source mass inflow per species (local rank) [mg/s]
    realArr QMassInSS_spec_glob;    ///< Source mass inflow per species (global) [mg/s]
    realArr QMassOutSS_spec;        ///< Sink mass outflow per species (local rank) [mg/s]
    realArr QMassOutSS_spec_glob;   ///< Sink mass outflow per species (global) [mg/s]

    // --- Reaction flux [mg/s] ---
    realArr QMassReaction_spec;       ///< Reaction-induced mass change rate per species (local) [mg/s]
    realArr QMassReaction_spec_glob;  ///< Reaction-induced mass change rate per species (global) [mg/s]

    // --- Rainfall solute input [mg/s] ---
    realArr QMassRain_spec;         ///< Rainfall solute mass input per species (local rank) [mg/s]
    realArr QMassRain_spec_glob;    ///< Rainfall solute mass input per species (global) [mg/s]

    // --- Temporary storage for reaction flux computation ---
    realArr tempMassBefore;         ///< Mass snapshot before reaction step [mg]
                                    ///< Used to compute the reaction flux as:
                                    ///< QMassReaction = (m_after - m_before) / dt_react

    /* ------------------------------------------------------------------ */
    /*  Cell count statistics                                              */
    /* ------------------------------------------------------------------ */
    int ncellsBC;       ///< Number of boundary condition cells (local rank)
    int ncellsBC_glob;  ///< Number of boundary condition cells (global, MPI-reduced)
    int ncellsSS;       ///< Number of source/sink cells (local rank)
    int ncellsSS_glob;  ///< Number of source/sink cells (global, MPI-reduced)

    /* ------------------------------------------------------------------ */
    /*  Pointers to external objects                                       */
    /* ------------------------------------------------------------------ */
    RTStateSW *rtsw;                   ///< Pointer to the surface water RT state
    State *state;                      ///< Pointer to the surface water hydraulic state
    Domain *dom;                       ///< Pointer to the surface water domain geometry
    SourceSinkData const *ss;          ///< Pointer to the source/sink data (const)
    ExternalBoundaries *ebc;           ///< Pointer to the external boundary conditions

    /* ------------------------------------------------------------------ */
    /*  Pre-allocated MPI communication buffers                            */
    /*  Promoted to class members to avoid dynamic allocation every        */
    /*  time step.  Each species contributes 8 double values; 2 extra      */
    /*  integers are appended for boundary and source/sink cell counts.    */
    /* ------------------------------------------------------------------ */
    std::vector<real> send_buf;   ///< Send buffer for MPI_Allreduce [total_mpi_vars]
    std::vector<real> recv_buf;   ///< Receive buffer for MPI_Allreduce [total_mpi_vars]
    int total_mpi_vars;           ///< Total number of variables in the MPI buffer

    /* ================================================================== */
    /*  Memory allocation                                                  */
    /* ================================================================== */

    /**
     * @brief Allocate all Kokkos views and MPI buffers for the given number
     *        of solute species.
     *
     * All views are initialised to zero.  The MPI send/receive buffers are
     * sized to hold 8 doubles per species plus 2 extra for cell counts.
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

        // Temporary storage for reaction before/after comparison [mg]
        tempMassBefore = realArr("tempMassBefore", n_mass_);

        // Zero-initialise all Kokkos views
        Kokkos::deep_copy(MassTot_spec, 0.0);
        Kokkos::deep_copy(MassTot_spec_glob, 0.0);
        Kokkos::deep_copy(MassExch_spec, 0.0);
        Kokkos::deep_copy(MassExch_spec_glob, 0.0);
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
        Kokkos::deep_copy(tempMassBefore, 0.0);

        // Pre-allocate MPI buffers: 8 double values per species + 2 cell counts
        total_mpi_vars = n_mass * 8 + 2;
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
        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            real m_local = 0.0;
            Kokkos::parallel_reduce("RT_MassBeforeReact", dom.nCell, KOKKOS_LAMBDA(int iGlob, real &tmp) {
                int ii = dom.getIndex(iGlob);
                if(!state.isnodata(ii)) {
                    // Mass = water_depth * concentration * cell_area * unit_conversion
                    // h4rtsw is the water depth used for transport [m]
                    tmp += state.h4rtsw(ii, 0) * rtsw.c(iSpec, ii, 1) * dom.cellArea() * 1000.0;
                } }, Kokkos::Sum<real>(m_local));
            tempMassBefore(iSpec) = m_local;
        }
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
        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            real m_after = 0.0;
            Kokkos::parallel_reduce("RT_MassAfterReact", dom.nCell, KOKKOS_LAMBDA(int iGlob, real &tmp) {
                int ii = dom.getIndex(iGlob);
                if(!state.isnodata(ii)) {
                    tmp += state.h4rtsw(ii, 0) * rtsw.c(iSpec, ii, 1) * dom.cellArea() * 1000.0;
                } }, Kokkos::Sum<real>(m_after));

            // Compute the mass change rate for this rank [mg/s]
            QMassReaction_spec(iSpec) = (m_after - tempMassBefore(iSpec)) / dt_react;
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
     * All local values are packed into a contiguous send buffer (8 doubles
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
        bool local_isRain = dom.isRain;

        // Zero out the send buffer for this time step
        std::fill(send_buf.begin(), send_buf.end(), 0.0);

        // ==================================================================
        // Phase 1: Kokkos-parallel reductions for spatial integrals
        // Uses standard parallel_reduce instead of atomic operations to
        // avoid CPU false sharing and lock contention.
        // ==================================================================
        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            // ==================================================================
            // 1. Total dissolved mass [mg]
            //    m_tot = sum( h * c * cellArea * 1000 )
            //    h = water depth [m], c = concentration [mg/L]
            // ==================================================================
            real m_tot = 0.0;
            Kokkos::parallel_reduce("RT_MassTot", dom.nCell, KOKKOS_LAMBDA(int iGlob, real &tmp) {
                int ii = dom.getIndex(iGlob);
                if(!state.isnodata(ii)) {
                    tmp += state.h(ii) * rtsw.c(iSpec, ii, 1) * dom.cellArea() * 1000.0;
                } }, Kokkos::Sum<real>(m_tot));
            MassTot_spec(iSpec) = m_tot;
            send_buf[iSpec * 8 + 0] = m_tot;

            // ==================================================================
            // 2. Surface-subsurface exchange flux [mg/s]
            //    ConQss is the concentration exchange rate [mg/L/s]
            //    multiplied by cell area and unit conversion.
            //    Only computed when coupled to the subsurface module.
            // ==================================================================
#if SERGHEI_RE_MODEL // [MERGED] Renamed from SERGHEI_SUBSURFACE_MODEL to SERGHEI_RE_MODEL
            real m_exch = 0.0;
            Kokkos::parallel_reduce("RT_MassExch", dom.nCell, KOKKOS_LAMBDA(int iGlob, real &tmp) {
                int ii = dom.getIndex(iGlob);
                if(!state.isnodata(ii)) {
                    tmp += rtsw.ConQss(iSpec, iGlob) * dom.cellArea() * 1000.0;
                } }, Kokkos::Sum<real>(m_exch));
            MassExch_spec(iSpec) = m_exch;
            send_buf[iSpec * 8 + 1] = m_exch;
#endif

            // ==================================================================
            // 3. Reaction-induced mass change rate [mg/s]
            //    Pre-computed by finalizeReactionFlux().
            // ==================================================================
            send_buf[iSpec * 8 + 2] = QMassReaction_spec(iSpec);

            // ==================================================================
            // 4. Rainfall solute input [mg/s]
            //    m_rain = sum( rainRate * rainConcentration * cellArea * 1000 )
            //    Only computed when rainfall is active and the species has
            //    a non-zero rain concentration.
            // ==================================================================
            real m_rain = 0.0;
            if (local_isRain)
            {
                real rain_con = rtsw.RainCon_arr(iSpec);  // Rain concentration for this species [mg/L]
                if (rain_con > 0.0)
                {
                    Kokkos::parallel_reduce("RT_MassRain", dom.nCell, KOKKOS_LAMBDA(int iGlob, real &tmp) {
                        int ii = dom.getIndex(iGlob);
                        if(!state.isnodata(ii)) {
                            // rainRate [m/s], rain_con [mg/L], area [m^2], 1000 [unit conversion]
                            tmp += ss.rainRate(iGlob) * rain_con * dom.cellArea() * 1000.0;
                        } }, Kokkos::Sum<real>(m_rain));
                }
            }
            QMassRain_spec(iSpec) = m_rain;
            send_buf[iSpec * 8 + 3] = m_rain;
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
            send_buf[iSpec * 8 + 4] = qin_bc;
            send_buf[iSpec * 8 + 5] = qout_bc;
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
            real qin_tmp = 0.0, qout_tmp = 0.0;

            for (int k = 0; k < ss.rtswss.size(); k++)
            {
                const RTSwSS &sw_ss_obj = ss.rtswss[k];
                if (sw_ss_obj.ncellsIT <= 0)
                    continue;

                // Determine the injection concentration for this species
                real C_inject = 0.0;
                if (sw_ss_obj.spec_sstype[iSpec] == 0)
                {
                    // Type 0: constant concentration [mg/L]
                    C_inject = sw_ss_obj.spec_ssval_const[iSpec];
                }
                else if (sw_ss_obj.spec_sstype[iSpec] == 1 || sw_ss_obj.spec_sstype[iSpec] == 2)
                {
                    // Type 1: linearly interpolated time series
                    // Type 2: step-wise constant time series
                    // Use const reference to avoid deep-copying the TimeSeries object
                    const auto &ts = sw_ss_obj.spec_ts[iSpec];
                    real t = dom.etime;
                    int ii_ts = 0;
                    while (ii_ts < ts.np - 1 && t >= ts.time(ii_ts + 1))
                        ii_ts++;

                    if (sw_ss_obj.spec_sstype[iSpec] == 2)
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

                // For type 3 sources (area-based), compute the injection flux
                // as concentration * number_of_cells * cell_area
                if (C_inject > 0.0 && sw_ss_obj.rtsstype == 3)
                {
                    qin_tmp += C_inject * sw_ss_obj.ncellsIT * dom.cellArea();
                }
            }
            QMassInSS_spec(iSpec) = qin_tmp;
            QMassOutSS_spec(iSpec) = qout_tmp;
            send_buf[iSpec * 8 + 6] = qin_tmp;
            send_buf[iSpec * 8 + 7] = qout_tmp;
        }

        // Pack cell counts at the end of the buffer
        send_buf[n_mass * 8 + 0] = (real)ncellsBC;
        send_buf[n_mass * 8 + 1] = (real)ncellsSS;

        // ==================================================================
        // Single MPI_Allreduce for all species and quantities
        // ==================================================================
        MPI_Allreduce(send_buf.data(), recv_buf.data(), total_mpi_vars, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);

        // Unpack the global results from the receive buffer
        for (int i = 0; i < n_mass; ++i)
        {
            MassTot_spec_glob(i) = recv_buf[i * 8 + 0];           // Total mass [mg]
            MassExch_spec_glob(i) = recv_buf[i * 8 + 1];          // Exchange flux [mg/s]
            QMassReaction_spec_glob(i) = recv_buf[i * 8 + 2];     // Reaction flux [mg/s]
            QMassRain_spec_glob(i) = recv_buf[i * 8 + 3];         // Rainfall flux [mg/s]
            QMassInBC_spec_glob(i) = recv_buf[i * 8 + 4];         // Boundary inflow [mg/s]
            QMassOutBC_spec_glob(i) = recv_buf[i * 8 + 5];        // Boundary outflow [mg/s]
            QMassInSS_spec_glob(i) = recv_buf[i * 8 + 6];         // Source inflow [mg/s]
            QMassOutSS_spec_glob(i) = recv_buf[i * 8 + 7];        // Sink outflow [mg/s]
        }
        ncellsBC_glob = (int)recv_buf[n_mass * 8 + 0];
        ncellsSS_glob = (int)recv_buf[n_mass * 8 + 1];
    }
};

#endif // SERGHEI_SURFACE_TRANSPORT
#endif // _RT_INTEGRATOR_SW_H_
