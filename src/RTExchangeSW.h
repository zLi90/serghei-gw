/* -*- mode: c++; c-default-style: "linux" -*- */

/**
 * @file RTExchangeSW.h
 * @brief MPI halo exchange class for surface water solute transport.
 *
 * Handles packing, sending, receiving, and unpacking of multi-species
 * concentration data at MPI domain boundaries for the surface water
 * reactive transport module. Supports transmissive and reflective
 * boundary conditions at domain edges.
 */

#ifndef _RT_EXCHANGE_SW_H_
#define _RT_EXCHANGE_SW_H_

#if SERGHEI_SURFACE_TRANSPORT

#include "define.h"
#include "mpi.h"
#include "Indexing.h"
#include "Domain.h"
#include "Exchange.h"
#include "RTStateSW.h"
#include "Parallel.h"

class RTExchangeSW
{

protected:
    MPI_Request sReq[2];  ///< MPI send request handles for non-blocking communication
    MPI_Request rReq[2];  ///< MPI receive request handles for non-blocking communication
    MPI_Status sStat[2];  ///< MPI send status objects
    MPI_Status rStat[2];  ///< MPI receive status objects

    int n_spec;            ///< Number of solute species to exchange [-]
    int nBufX, nBufY;     ///< Buffer sizes for X-direction and Y-direction halo exchange

    realArr haloSendBufS, haloSendBufN, haloSendBufW, haloSendBufE;  ///< Send buffers for South, North, West, East halos
    realArr haloRecvBufS, haloRecvBufN, haloRecvBufW, haloRecvBufE;  ///< Receive buffers for South, North, West, East halos

    Kokkos::Timer timerExchange;  ///< Timer for measuring exchange duration

public:
    real exchangeTime = 0;  ///< Cumulative MPI exchange time [s]

    /// Allocate halo exchange buffers. Must be called after n_mass is known.
    /// @param dom  Domain with grid dimensions
    /// @param n_spec_  Number of solute species
    inline void allocate(Domain &dom, int n_spec_)
    {
        n_spec = n_spec_;
        nBufX = n_spec * dom.ny * dom.hc;  ///< X-direction buffer: species × rows × halo cells
        nBufY = n_spec * dom.hc * dom.nx;  ///< Y-direction buffer: species × halo cells × columns

        haloSendBufS = realArr("RTSendBufS", nBufY);
        haloSendBufN = realArr("RTSendBufN", nBufY);
        haloSendBufW = realArr("RTSendBufW", nBufX);
        haloSendBufE = realArr("RTSendBufE", nBufX);
        haloRecvBufS = realArr("RTRecvBufS", nBufY);
        haloRecvBufN = realArr("RTRecvBufN", nBufY);
        haloRecvBufW = realArr("RTRecvBufW", nBufX);
        haloRecvBufE = realArr("RTRecvBufE", nBufX);
    }

    // =====================================================================
    // X-direction pack: extract west/east boundary concentration data
    // =====================================================================
    /// Pack concentration data from the west and east halo regions for X-direction MPI exchange.
    /// Reads from rt.c(iSpec, iGlob, iTime) for the first/last hc columns.
    /// @param dom  Domain with grid layout
    /// @param rt   Surface water reactive transport state
    /// @param iTime  Time index: 0=previous, 1=current
    inline void haloPack_x_rt(Domain &dom, RTStateSW &rt, int iTime)
    {
        const int nGlobX = dom.ny * dom.hc;
        auto local_sendBufW = haloSendBufW;
        auto local_sendBufE = haloSendBufE;
        auto local_c = rt.c;
        int local_n_spec = n_spec;
        int local_ny = dom.ny;
        int local_nx = dom.nx;
        int local_hc = dom.hc;

        Kokkos::parallel_for("RT_haloPack_x", local_n_spec * nGlobX, KOKKOS_LAMBDA(const int idx) {
            int iSpec = idx / nGlobX;
            int iGlob = idx % nGlobX;
            int rx, ry;
            unpackIndicesUniformGrid(iGlob, local_ny, local_hc, ry, rx);

            // West boundary: first hc interior columns
            int ii_west = (local_hc + ry) * (local_nx + 2 * local_hc) + local_hc + rx;
            // East boundary: last hc interior columns
            int ii_east = (local_hc + ry) * (local_nx + 2 * local_hc) + local_nx + rx;

            local_sendBufW(iSpec * nGlobX + iGlob) = local_c(iSpec, ii_west, iTime);
            local_sendBufE(iSpec * nGlobX + iGlob) = local_c(iSpec, ii_east, iTime); });
    }

    // =====================================================================
    // Y-direction pack: extract north/south boundary concentration data
    // =====================================================================
    /// Pack concentration data from the north and south halo regions for Y-direction MPI exchange.
    /// @param dom  Domain with grid layout
    /// @param rt   Surface water reactive transport state
    /// @param iTime  Time index: 0=previous, 1=current
    inline void haloPack_y_rt(Domain &dom, RTStateSW &rt, int iTime)
    {
        const int nGlobY = dom.hc * dom.nx;
        auto local_sendBufN = haloSendBufN;
        auto local_sendBufS = haloSendBufS;
        auto local_c = rt.c;
        int local_n_spec = n_spec;
        int local_ny = dom.ny;
        int local_nx = dom.nx;
        int local_hc = dom.hc;

        Kokkos::parallel_for("RT_haloPack_y", local_n_spec * nGlobY, KOKKOS_LAMBDA(const int idx) {
            int iSpec = idx / nGlobY;
            int iGlob = idx % nGlobY;
            int rx, ry;
            unpackIndicesUniformGrid(iGlob, local_hc, local_nx, ry, rx);

            // North boundary: first hc interior rows
            int ii_north = (local_hc + ry) * (local_nx + 2 * local_hc) + local_hc + rx;
            // South boundary: last hc interior rows
            int ii_south = (local_ny + ry) * (local_nx + 2 * local_hc) + local_hc + rx;

            local_sendBufN(iSpec * nGlobY + iGlob) = local_c(iSpec, ii_north, iTime);
            local_sendBufS(iSpec * nGlobY + iGlob) = local_c(iSpec, ii_south, iTime); });
    }

    // =====================================================================
    // X-direction unpack: write received data into halo cells
    // =====================================================================
    /// Unpack received concentration data into the west and east halo cells for X-direction.
    /// @param dom  Domain with grid layout
    /// @param rt   Surface water reactive transport state
    /// @param iTime  Time index: 0=previous, 1=current
    inline void haloUnpack_x_rt(Domain &dom, RTStateSW &rt, int iTime)
    {
        const int nGlobX = dom.ny * dom.hc;

        auto local_recvBufW = haloRecvBufW;
        auto local_recvBufE = haloRecvBufE;
        auto local_c = rt.c;
        int local_n_spec = n_spec;
        int local_ny = dom.ny;
        int local_nx = dom.nx;
        int local_hc = dom.hc;

        Kokkos::parallel_for("RT_haloUnpack_x", local_n_spec * nGlobX, KOKKOS_LAMBDA(const int idx) {
            int iSpec = idx / nGlobX;
            int iGlob = idx % nGlobX;

            int rx, ry;
            unpackIndicesUniformGrid(iGlob, local_ny, local_hc, ry, rx);

            // West halo cells: columns before interior
            int ii_west = (local_hc + ry) * (local_nx + 2 * local_hc) + local_hc - rx - 1;
            // East halo cells: columns after interior
            int ii_east = (local_hc + ry) * (local_nx + 2 * local_hc) + local_hc + local_nx + rx;

            local_c(iSpec, ii_west, iTime) = local_recvBufW(iSpec * nGlobX + iGlob);
            local_c(iSpec, ii_east, iTime) = local_recvBufE(iSpec * nGlobX + iGlob); });
    }

    // =====================================================================
    // Y-direction unpack: write received data into halo cells
    // =====================================================================
    /// Unpack received concentration data into the north and south halo cells for Y-direction.
    /// @param dom  Domain with grid layout
    /// @param rt   Surface water reactive transport state
    /// @param iTime  Time index: 0=previous, 1=current
    inline void haloUnpack_y_rt(Domain &dom, RTStateSW &rt, int iTime)
    {
        const int nGlobY = dom.hc * dom.nx;

        auto local_recvBufN = haloRecvBufN;
        auto local_recvBufS = haloRecvBufS;
        auto local_c = rt.c;
        int local_n_spec = n_spec;
        int local_ny = dom.ny;
        int local_nx = dom.nx;
        int local_hc = dom.hc;

        Kokkos::parallel_for("RT_haloUnpack_y", local_n_spec * nGlobY, KOKKOS_LAMBDA(const int idx) {
            int iSpec = idx / nGlobY;
            int iGlob = idx % nGlobY;

            int rx, ry;
            unpackIndicesUniformGrid(iGlob, local_hc, local_nx, ry, rx);

            // North halo cells: rows before interior
            int ii_north = (local_hc - ry - 1) * (local_nx + 2 * local_hc) + local_hc + rx;
            // South halo cells: rows after interior
            int ii_south = (local_hc + local_ny + ry) * (local_nx + 2 * local_hc) + local_hc + rx;

            local_c(iSpec, ii_north, iTime) = local_recvBufN(iSpec * nGlobY + iGlob);
            local_c(iSpec, ii_south, iTime) = local_recvBufS(iSpec * nGlobY + iGlob); });
    }

    // =====================================================================
    // X-direction MPI send/receive
    // =====================================================================
    /// Perform non-blocking MPI send/receive for X-direction halo exchange.
    /// @return MPI error code (0 on success)
    inline int haloExchange_x_rt(Domain &dom, Parallel &par)
    {
        int ierr = 1;

        if (par.nproc_x > 1)
        {
            Kokkos::fence();
            ierr = MPI_Irecv(haloRecvBufW.data(), nBufX, SERGHEI_MPI_REAL, par.neigh(1, 0), 0, MPI_COMM_WORLD, &rReq[0]);
            ierr = MPI_Irecv(haloRecvBufE.data(), nBufX, SERGHEI_MPI_REAL, par.neigh(1, 2), 1, MPI_COMM_WORLD, &rReq[1]);
            ierr = MPI_Isend(haloSendBufW.data(), nBufX, SERGHEI_MPI_REAL, par.neigh(1, 0), 1, MPI_COMM_WORLD, &sReq[0]);
            ierr = MPI_Isend(haloSendBufE.data(), nBufX, SERGHEI_MPI_REAL, par.neigh(1, 2), 0, MPI_COMM_WORLD, &sReq[1]);

            ierr = MPI_Waitall(2, sReq, sStat);
            ierr = MPI_Waitall(2, rReq, rStat);
        }

        if (dom.BCtype == SWE_BC_TRANSMISSIVE)
        {
            if (dom.iW)
                haloTransmissive(nBufX, haloSendBufW, haloRecvBufW);
            if (dom.iE)
                haloTransmissive(nBufX, haloSendBufE, haloRecvBufE);
        }
        if (dom.BCtype == SWE_BC_REFLECTIVE)
        {
            if (dom.iW)
                haloReflective(nBufX, haloRecvBufW);
            if (dom.iE)
                haloReflective(nBufX, haloRecvBufE);
        }
        return ierr;
    }

    // =====================================================================
    // Y-direction MPI send/receive
    // =====================================================================
    /// Perform non-blocking MPI send/receive for Y-direction halo exchange.
    /// @return MPI error code (0 on success)
    inline int haloExchange_y_rt(Domain &dom, Parallel &par)
    {
        int ierr = 1;

        if (par.nproc_y > 1)
        {
            Kokkos::fence();
            ierr = MPI_Irecv(haloRecvBufN.data(), nBufY, SERGHEI_MPI_REAL, par.neigh(0, 1), 0, MPI_COMM_WORLD, &rReq[0]);
            ierr = MPI_Irecv(haloRecvBufS.data(), nBufY, SERGHEI_MPI_REAL, par.neigh(2, 1), 1, MPI_COMM_WORLD, &rReq[1]);
            ierr = MPI_Isend(haloSendBufN.data(), nBufY, SERGHEI_MPI_REAL, par.neigh(0, 1), 1, MPI_COMM_WORLD, &sReq[0]);
            ierr = MPI_Isend(haloSendBufS.data(), nBufY, SERGHEI_MPI_REAL, par.neigh(2, 1), 0, MPI_COMM_WORLD, &sReq[1]);

            ierr = MPI_Waitall(2, sReq, sStat);
            ierr = MPI_Waitall(2, rReq, rStat);
        }

        if (dom.BCtype == SWE_BC_TRANSMISSIVE)
        {
            if (dom.iN)
                haloTransmissive(nBufY, haloSendBufN, haloRecvBufN);
            if (dom.iS)
                haloTransmissive(nBufY, haloSendBufS, haloRecvBufS);
        }
        if (dom.BCtype == SWE_BC_REFLECTIVE)
        {
            if (dom.iN)
                haloReflective(nBufY, haloRecvBufN);
            if (dom.iS)
                haloReflective(nBufY, haloRecvBufS);
        }
        return ierr;
    }

    // =====================================================================
    // Boundary condition helpers
    // =====================================================================
    /// Apply reflective boundary: set halo concentration to zero (no-flux condition for solute)
    inline void haloReflective(const int ncells, realArr &haloRecv)
    {
        auto local_recv = haloRecv;
        Kokkos::parallel_for("RT_haloReflective", ncells, KOKKOS_LAMBDA(const int i) { local_recv(i) = 0.0; });
    }

    /// Apply transmissive boundary: copy interior concentration to halo (zero-gradient outflow)
    inline void haloTransmissive(const int ncells, realArr const &haloSend, realArr &haloRecv)
    {
        auto local_send = haloSend;
        auto local_recv = haloRecv;
        Kokkos::parallel_for("RT_haloTransmissive", ncells, KOKKOS_LAMBDA(const int i) { local_recv(i) = local_send(i); });
    }

    // =====================================================================
    // High-level exchange interface
    // =====================================================================
    /// Perform full MPI halo exchange for surface water concentrations in both X and Y directions.
    /// @param rt    Surface water reactive transport state
    /// @param dom   Domain with grid layout
    /// @param par   Parallel decomposition information
    /// @param iTime Time index: 0=previous, 1=current
    inline void exchangeMPI(RTStateSW &rt, Domain &dom, Parallel &par, int iTime)
    {
        // X-direction exchange
        haloPack_x_rt(dom, rt, iTime);
        haloExchange_x_rt(dom, par);
        haloUnpack_x_rt(dom, rt, iTime);

        // Y-direction exchange
        haloPack_y_rt(dom, rt, iTime);
        haloExchange_y_rt(dom, par);
        haloUnpack_y_rt(dom, rt, iTime);
    }
};

#endif // SERGHEI_SURFACE_TRANSPORT
#endif // _RT_EXCHANGE_SW_H_
