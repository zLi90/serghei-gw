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
    realArr LiquidMassTot_spec;       ///< Total dissolved mass per species (local rank) [mg]
    realArr LiquidMassTot_spec_glob;  ///< Total dissolved mass per species (global, MPI-reduced) [mg]
    realArr SolidMassTot_spec;        ///< Total adsorbed mass per species (local rank) [mg]
    realArr SolidMassTot_spec_glob;   ///< Total adsorbed mass per species (global, MPI-reduced) [mg]

    // --- Surface-subsurface exchange flux [mg/s] ---
    realArr MassExch_spec;           ///< Net exchange flux per species (local rank) [mg/s]
    realArr MassExch_spec_glob;      ///< Net exchange flux per species (global, MPI-reduced) [mg/s]

    // --- Boundary fluxes [mg/s] ---
    realArr QMassInBC_spec;          ///< Boundary mass inflow per species (local rank) [mg/s]
    realArr QMassOutBC_spec;         ///< Boundary mass outflow per species (local rank) [mg/s]
    realArr QMassInBC_spec_glob;     ///< Boundary mass inflow per species (global) [mg/s]
    realArr QMassOutBC_spec_glob;    ///< Boundary mass outflow per species (global) [mg/s]

    // --- Source/sink fluxes [mg/s] ---
    realArr QMassInSS_spec;          ///< Source mass inflow per species (local rank) [mg/s]
    realArr QMassOutSS_spec;         ///< Sink mass outflow per species (local rank) [mg/s]
    realArr QMassInSS_spec_glob;     ///< Source mass inflow per species (global) [mg/s]
    realArr QMassOutSS_spec_glob;    ///< Sink mass outflow per species (global) [mg/s]

    // --- Reaction flux [mg/s] ---
    realArr QMassReaction_spec;       ///< Reaction-induced mass change rate per species (local) [mg/s]
    realArr QMassReaction_spec_glob;  ///< Reaction-induced mass change rate per species (global) [mg/s]

    // --- Temporary storage for reaction flux computation ---
    realArr tempLiquidMassBefore;     ///< Liquid mass snapshot before reaction step [mg]
                                     ///< Used to compute the reaction flux as:
                                     ///< QMassReaction = (m_after - m_before) / dt_react

    /* ------------------------------------------------------------------ */
    /*  Pointers to external objects                                       */
    /* ------------------------------------------------------------------ */
    SourceSink *ss;              ///< Pointer to the source/sink manager
    RTStateGW *rt;              ///< Pointer to the groundwater reactive-transport state
    GwDomain *gdom;             ///< Pointer to the groundwater domain geometry
    std::vector<RTBCGW> *rtbc;  ///< Pointer to the vector of boundary condition objects

    /* ------------------------------------------------------------------ */
    /*  Cell count statistics                                              */
    /* ------------------------------------------------------------------ */
    int ncellsBC;       ///< Number of boundary condition cells (local rank)
    int ncellsBC_glob;  ///< Number of boundary condition cells (global, MPI-reduced)
    int ncellsIT;       ///< Number of interior (integration) cells (local rank)
    int ncellsIT_glob;  ///< Number of interior cells (global, MPI-reduced)

    /* ------------------------------------------------------------------ */
    /*  Pre-allocated MPI communication buffers                            */
    /*  Promoted to class members to avoid dynamic allocation every        */
    /*  time step.  Each species contributes 8 double values; one extra    */
    /*  integer is appended for the boundary cell count.                   */
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
        MassExch_spec = realArr("MassExch_spec", n_mass);
        MassExch_spec_glob = realArr("MassExch_spec_glob", n_mass);

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

        // Temporary storage for reaction before/after comparison [mg]
        tempLiquidMassBefore = realArr("tempLiquidMassBefore", n_mass);

        // Zero-initialise all Kokkos views
        Kokkos::deep_copy(LiquidMassTot_spec, 0.0);
        Kokkos::deep_copy(LiquidMassTot_spec_glob, 0.0);
        Kokkos::deep_copy(SolidMassTot_spec, 0.0);
        Kokkos::deep_copy(SolidMassTot_spec_glob, 0.0);
        Kokkos::deep_copy(MassExch_spec, 0.0);
        Kokkos::deep_copy(MassExch_spec_glob, 0.0);
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
        Kokkos::deep_copy(tempLiquidMassBefore, 0.0);

        // Pre-allocate MPI buffers: 8 double values per species + 1 cell count
        total_mpi_vars = n_mass * 8 + 1;
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
    }

    /* ================================================================== */
    /*  Reaction flux computation (two-phase snapshot approach)            */
    /* ================================================================== */

    /**
     * @brief Record the total liquid-phase mass of each species BEFORE the
     *        chemical reaction step.
     *
     * This method computes the domain-wide dissolved mass using the current
     * water content and concentration fields.  The result is stored in
     * tempLiquidMassBefore for later comparison with the post-reaction mass
     * in finalizeReactionFlux().
     *
     * Mass computation per cell:
     *   m = water_content [-] * concentration [mg/L] * cell_volume [m^3] * 1000 [mg/mg]
     *
     * @param rt    Groundwater reactive-transport state (concentrations).
     * @param gw    Groundwater state (water content).
     * @param gdom  Groundwater domain (cell geometry).
     */
    void recordLiquidMassBefore(RTStateGW const &rt, GwState const &gw, GwDomain const &gdom)
    {
        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            real m_local = 0.0;
            Kokkos::parallel_reduce("RTGW_MassBeforeReact", gdom.nCell, KOKKOS_LAMBDA(int idx, real &tmp) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idx, kk, jj, ii);
                iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
                real Vol = gdom.dx * gdom.dy * gdom.dz(iGlob);
                // Liquid mass = water_content * concentration * volume * unit conversion
                tmp += gw.wc(iGlob, 1) * rt.c(iSpec, iGlob, 1) * Vol * 1000.0; }, Kokkos::Sum<real>(m_local));
            tempLiquidMassBefore(iSpec) = m_local;
        }
    }

    /**
     * @brief Compute the reaction-induced mass change rate for each species
     *        AFTER the chemical reaction step.
     *
     * The reaction flux is computed as:
     *   QMassReaction = (m_after - m_before) / dt_react
     *
     * where m_before was recorded by recordLiquidMassBefore() and m_after is
     * the current post-reaction liquid mass.  A positive value indicates net
     * production; a negative value indicates net consumption.
     *
     * @param rt        Groundwater reactive-transport state (post-reaction concentrations).
     * @param gw        Groundwater state (water content).
     * @param gdom      Groundwater domain (cell geometry).
     * @param dt_react  Duration of the reaction time step [s].
     */
    void finalizeReactionFlux(RTStateGW const &rt, GwState const &gw, GwDomain const &gdom, real dt_react)
    {
        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            real m_after = 0.0;
            Kokkos::parallel_reduce("RTGW_MassAfterReact", gdom.nCell, KOKKOS_LAMBDA(int idx, real &tmp) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idx, kk, jj, ii);
                iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
                real Vol = gdom.dx * gdom.dy * gdom.dz(iGlob);
                tmp += gw.wc(iGlob, 1) * rt.c(iSpec, iGlob, 1) * Vol * 1000.0; }, Kokkos::Sum<real>(m_after));

            // Compute the mass change rate [mg/s]
            QMassReaction_spec(iSpec) = (m_after - tempLiquidMassBefore(iSpec)) / dt_react;
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
    void integrate(GwState const &gw, RTStateGW const &rt, GwDomain const &gdom, std::vector<RTBCGW> &rtbc)
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
            //    m_liquid = sum( wc * c * dx * dy * dz * 1000 )
            //    where wc = volumetric water content [-], c = concentration [mg/L]
            // ==================================================================
            real liquid_tmp = 0.0;
            Kokkos::parallel_reduce("RTGW_Liquid", gdom.nCell, KOKKOS_LAMBDA(int idx, real &tmp) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idx, kk, jj, ii);
                iGlob = (gdom.hc+kk)*gdom.nxhc*gdom.nyhc + (gdom.hc+jj)*gdom.nxhc + ii + gdom.hc;
                real Vol = gdom.dx * gdom.dy * gdom.dz(iGlob);
                tmp += gw.wc(iGlob, 1) * rt.c(iSpec, iGlob, 1) * Vol * 1000.0; }, Kokkos::Sum<real>(liquid_tmp));
            LiquidMassTot_spec(iSpec) = liquid_tmp;
            send_buf[iSpec * 8 + 0] = liquid_tmp;

            // ==================================================================
            // 2. Total solid-phase (adsorbed) mass [mg]
            //    m_solid = sum( c_solid * dx * dy * dz * 1000 )
            // ==================================================================
            real solid_tmp = 0.0;
            Kokkos::parallel_reduce("RTGW_Solid", gdom.nCell, KOKKOS_LAMBDA(int idx, real &tmp) {
                int ii, jj, kk, iGlob;
                gdom.unpackIndices(idx, kk, jj, ii);
                iGlob = (gdom.hc + kk) * gdom.nxhc * gdom.nyhc + (gdom.hc + jj) * gdom.nxhc + ii + gdom.hc;
                real Vol = gdom.dx * gdom.dy * gdom.dz(iGlob);
                tmp += rt.c_solid(iSpec, iGlob, 1) * Vol * 1000.0; }, Kokkos::Sum<real>(solid_tmp));
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
            send_buf[iSpec * 8 + 2] = exch_tmp;

            // ==================================================================
            // 4. Boundary mass inflow and outflow [mg/s]
            //    Accumulated from all boundary condition objects for this
            //    species.  Only boundaries that have active cells on this
            //    rank contribute.
            // ==================================================================
            real qin_bc = 0.0;
            real qout_bc = 0.0;
            for (int k = 0; k < rtbc.size(); k++)
            {
                // Only include contributions from boundaries that have cells
                // on this MPI rank
                if (rtbc[k].ncellsBC > 0)
                {
                    qin_bc += rtbc[k].QMassInflow;
                    qout_bc += rtbc[k].QMassOutflow;
                }
            }
            QMassInBC_spec(iSpec) = qin_bc;
            QMassOutBC_spec(iSpec) = qout_bc;
            send_buf[iSpec * 8 + 3] = qin_bc;
            send_buf[iSpec * 8 + 4] = qout_bc;

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
            real tscf = DEFAULT_TSCF;

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
                if (transp_rate < 0.0) {
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
            send_buf[iSpec * 8 + 5] = ss_in_tmp;
            send_buf[iSpec * 8 + 6] = ss_out_tmp;

            // ==================================================================
            // 6. Reaction-induced mass change rate [mg/s]
            //    Pre-computed by finalizeReactionFlux().
            // ==================================================================
            send_buf[iSpec * 8 + 7] = QMassReaction_spec(iSpec);
        }

        // ==================================================================
        // MPI global reduction: single MPI_Allreduce replaces the previous
        // per-variable aggregation, dramatically reducing communication
        // overhead for large numbers of species.
        // ==================================================================
        MPI_Allreduce(send_buf.data(), recv_buf.data(), total_mpi_vars, SERGHEI_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);

        // Unpack the global results from the receive buffer
        for (int iSpec = 0; iSpec < n_mass; ++iSpec)
        {
            LiquidMassTot_spec_glob(iSpec) = recv_buf[iSpec * 8 + 0];   // Total liquid mass [mg]
            SolidMassTot_spec_glob(iSpec) = recv_buf[iSpec * 8 + 1];    // Total solid mass [mg]
            MassExch_spec_glob(iSpec) = recv_buf[iSpec * 8 + 2];        // Exchange flux [mg/s]
            QMassInBC_spec_glob(iSpec) = recv_buf[iSpec * 8 + 3];       // Boundary inflow [mg/s]
            QMassOutBC_spec_glob(iSpec) = recv_buf[iSpec * 8 + 4];      // Boundary outflow [mg/s]
            QMassInSS_spec_glob(iSpec) = recv_buf[iSpec * 8 + 5];       // Source inflow [mg/s]
            QMassOutSS_spec_glob(iSpec) = recv_buf[iSpec * 8 + 6];      // Sink outflow [mg/s]
            QMassReaction_spec_glob(iSpec) = recv_buf[iSpec * 8 + 7];   // Reaction flux [mg/s]
        }
        ncellsBC_glob = (int)recv_buf[n_mass * 8];

        // Note: MPI_Barrier removed to avoid unnecessary synchronisation stalls
    }
};
#endif
#endif
