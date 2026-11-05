/* -*- mode: c++; c-default-style: "linux" -*- */
#ifndef _RT_EXCHANGE_SW_H_
#define _RT_EXCHANGE_SW_H_

#if SERGHEI_SURFACE_TRANSPORT

#include "define.h"
#include "mpi.h"
#include "Domain.h"
#include "RTStateSW.h"
#include "Parallel.h"
#include "MPIUtils.h"
#include <algorithm>
#include <vector>

/** Independent two-cell halo exchange for surface-water solute concentration.
 *
 * The hydrodynamic Domain keeps hc=1.  X is exchanged first over physical
 * rows; Y is then exchanged over the complete RT width, including the newly
 * received X halos.  This also propagates all four corner blocks.
 */
class RTExchangeSW
{
protected:
    MPI_Request req[4] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL,
                          MPI_REQUEST_NULL, MPI_REQUEST_NULL};
    MPI_Status stat[4];

    int n_spec = 0;
    int nBufX = 0, nBufY = 0;
    realArr haloSendBufS, haloSendBufN, haloSendBufW, haloSendBufE;
    realArr haloRecvBufS, haloRecvBufN, haloRecvBufW, haloRecvBufE;
    // Persistent host mirrors avoid allocations and redundant host copies in
    // every halo exchange. They remain valid for the lifetime of the solver.
    Kokkos::View<real *, RTMPIHostSpace> h_sendS, h_sendN, h_sendW, h_sendE;
    Kokkos::View<real *, RTMPIHostSpace> h_recvS, h_recvN, h_recvW, h_recvE;
    Kokkos::Timer timerExchange;
    Kokkos::Timer timerStage;

public:
    real exchangeTime = 0;

    inline void allocate(Domain &dom, int n_spec_)
    {
        const int requiredBufX = n_spec_ * dom.ny * RTStateSW::rt_hc;
        const int requiredBufY = n_spec_ * RTStateSW::rt_hc *
                                 (dom.nx + 2 * RTStateSW::rt_hc);
        if (n_spec == n_spec_ && nBufX == requiredBufX && nBufY == requiredBufY &&
            haloSendBufS.extent(0) == static_cast<std::size_t>(requiredBufY) &&
            haloSendBufN.extent(0) == static_cast<std::size_t>(requiredBufY) &&
            haloRecvBufS.extent(0) == static_cast<std::size_t>(requiredBufY) &&
            haloRecvBufN.extent(0) == static_cast<std::size_t>(requiredBufY) &&
            haloSendBufW.extent(0) == static_cast<std::size_t>(requiredBufX) &&
            haloSendBufE.extent(0) == static_cast<std::size_t>(requiredBufX) &&
            haloRecvBufW.extent(0) == static_cast<std::size_t>(requiredBufX) &&
            haloRecvBufE.extent(0) == static_cast<std::size_t>(requiredBufX) &&
            h_sendS.extent(0) == static_cast<std::size_t>(requiredBufY) &&
            h_sendN.extent(0) == static_cast<std::size_t>(requiredBufY) &&
            h_recvS.extent(0) == static_cast<std::size_t>(requiredBufY) &&
            h_recvN.extent(0) == static_cast<std::size_t>(requiredBufY) &&
            h_sendW.extent(0) == static_cast<std::size_t>(requiredBufX) &&
            h_sendE.extent(0) == static_cast<std::size_t>(requiredBufX) &&
            h_recvW.extent(0) == static_cast<std::size_t>(requiredBufX) &&
            h_recvE.extent(0) == static_cast<std::size_t>(requiredBufX))
            return;
        n_spec = n_spec_;
        nBufX = requiredBufX;
        nBufY = requiredBufY;
        haloSendBufS = realArr("RTSendBufS", nBufY);
        haloSendBufN = realArr("RTSendBufN", nBufY);
        haloSendBufW = realArr("RTSendBufW", nBufX);
        haloSendBufE = realArr("RTSendBufE", nBufX);
        haloRecvBufS = realArr("RTRecvBufS", nBufY);
        haloRecvBufN = realArr("RTRecvBufN", nBufY);
        haloRecvBufW = realArr("RTRecvBufW", nBufX);
        haloRecvBufE = realArr("RTRecvBufE", nBufX);
        h_sendS = Kokkos::View<real *, RTMPIHostSpace>("RT_h_sendS", nBufY);
        h_sendN = Kokkos::View<real *, RTMPIHostSpace>("RT_h_sendN", nBufY);
        h_recvS = Kokkos::View<real *, RTMPIHostSpace>("RT_h_recvS", nBufY);
        h_recvN = Kokkos::View<real *, RTMPIHostSpace>("RT_h_recvN", nBufY);
        h_sendW = Kokkos::View<real *, RTMPIHostSpace>("RT_h_sendW", nBufX);
        h_sendE = Kokkos::View<real *, RTMPIHostSpace>("RT_h_sendE", nBufX);
        h_recvW = Kokkos::View<real *, RTMPIHostSpace>("RT_h_recvW", nBufX);
        h_recvE = Kokkos::View<real *, RTMPIHostSpace>("RT_h_recvE", nBufX);
    }

    inline void haloPack_x_rt(Domain &dom, RTStateSW &rt, int iTime)
    {
        const int hc = RTStateSW::rt_hc;
        const int perSpec = dom.ny * hc;
        auto sendW = haloSendBufW; auto sendE = haloSendBufE; auto c = rt.c;
        const int nx = dom.nx, nxhc = rt.rt_nxhc, ns = n_spec;
        Kokkos::parallel_for("RT_haloPack_x2", ns * perSpec, KOKKOS_LAMBDA(int idx) {
            const int s = idx / perSpec;
            const int q = idx % perSpec;
            const int j = q / hc;
            const int l = q % hc;
            sendW(idx) = c(s, (j + hc) * nxhc + hc + l, iTime);
            sendE(idx) = c(s, (j + hc) * nxhc + hc + nx - hc + l, iTime);
        });
    }

    inline void haloUnpack_x_rt(Domain &dom, RTStateSW &rt, int iTime)
    {
        const int hc = RTStateSW::rt_hc;
        const int perSpec = dom.ny * hc;
        auto recvW = haloRecvBufW; auto recvE = haloRecvBufE; auto c = rt.c;
        const int nx = dom.nx, nxhc = rt.rt_nxhc, ns = n_spec;
        Kokkos::parallel_for("RT_haloUnpack_x2", ns * perSpec, KOKKOS_LAMBDA(int idx) {
            const int s = idx / perSpec;
            const int q = idx % perSpec;
            const int j = q / hc;
            const int l = q % hc;
            c(s, (j + hc) * nxhc + l, iTime) = recvW(idx);
            c(s, (j + hc) * nxhc + hc + nx + l, iTime) = recvE(idx);
        });
    }

    inline void haloPack_y_rt(Domain &dom, RTStateSW &rt, int iTime)
    {
        const int hc = RTStateSW::rt_hc;
        const int width = rt.rt_nxhc;
        const int perSpec = hc * width;
        auto sendN = haloSendBufN; auto sendS = haloSendBufS; auto c = rt.c;
        const int ny = dom.ny, ns = n_spec;
        Kokkos::parallel_for("RT_haloPack_y2", ns * perSpec, KOKKOS_LAMBDA(int idx) {
            const int s = idx / perSpec;
            const int q = idx % perSpec;
            const int l = q / width;
            const int i = q % width;
            sendN(idx) = c(s, (hc + l) * width + i, iTime);
            sendS(idx) = c(s, (hc + ny - hc + l) * width + i, iTime);
        });
    }

    inline void haloUnpack_y_rt(Domain &dom, RTStateSW &rt, int iTime)
    {
        const int hc = RTStateSW::rt_hc;
        const int width = rt.rt_nxhc;
        const int perSpec = hc * width;
        auto recvN = haloRecvBufN; auto recvS = haloRecvBufS; auto c = rt.c;
        const int ny = dom.ny, ns = n_spec;
        Kokkos::parallel_for("RT_haloUnpack_y2", ns * perSpec, KOKKOS_LAMBDA(int idx) {
            const int s = idx / perSpec;
            const int q = idx % perSpec;
            const int l = q / width;
            const int i = q % width;
            c(s, l * width + i, iTime) = recvN(idx);
            c(s, (hc + ny + l) * width + i, iTime) = recvS(idx);
        });
    }

    inline void fillPhysicalX(Domain &dom, RTStateSW &rt, int iTime)
    {
        const int hc = RTStateSW::rt_hc, perSpec = dom.ny * hc;
        auto rw = haloRecvBufW; auto re = haloRecvBufE; auto c = rt.c;
        const int nx = dom.nx, width = rt.rt_nxhc, ns = n_spec;
        const bool west = dom.iW, east = dom.iE;
        Kokkos::parallel_for("RT_physical_x2", ns * perSpec, KOKKOS_LAMBDA(int idx) {
            const int s=idx/perSpec, q=idx%perSpec, j=q/hc;
            // Solute concentration is a scalar: both transmissive and
            // reflective hydraulic walls use an even reflection (zero normal
            // gradient), never a zero-concentration ghost value.
            if (west) rw(idx)=c(s,(j+hc)*width+hc,iTime);
            if (east) re(idx)=c(s,(j+hc)*width+hc+nx-1,iTime);
        });
    }

    inline void fillPhysicalY(Domain &dom, RTStateSW &rt, int iTime)
    {
        const int hc = RTStateSW::rt_hc, width = rt.rt_nxhc, perSpec=hc*width;
        auto rn=haloRecvBufN; auto rs=haloRecvBufS; auto c=rt.c;
        const int ny=dom.ny, ns=n_spec; const bool north=dom.iN, south=dom.iS;
        Kokkos::parallel_for("RT_physical_y2", ns*perSpec, KOKKOS_LAMBDA(int idx) {
            const int s=idx/perSpec, i=idx%width;
            if (north) rn(idx)=c(s,hc*width+i,iTime);
            if (south) rs(idx)=c(s,(hc+ny-1)*width+i,iTime);
        });
    }

    /// Single-rank / per-direction replacement for fillPhysicalX + haloUnpack_x_rt.
    /// Writes the zero-normal-gradient physical-boundary ghosts directly,
    /// skipping the recv-buffer round trip.  Bit-for-bit identical to the
    /// fillPhysicalX -> haloUnpack_x_rt composition it replaces:
    /// ghost[l] = recvW(idx) = c(interior first column), l = 0..hc-1.
    inline void fillPhysicalXDirect(Domain &dom, RTStateSW &rt, int iTime)
    {
        const int hc = RTStateSW::rt_hc;
        const int perSpec = dom.ny * hc;
        auto c = rt.c;
        const int nx = dom.nx, width = rt.rt_nxhc, ns = n_spec;
        const bool west = dom.iW, east = dom.iE;
        Kokkos::parallel_for("RT_physical_x_direct", ns * perSpec, KOKKOS_LAMBDA(int idx) {
            const int s = idx / perSpec;
            const int q = idx % perSpec;
            const int j = q / hc;
            const int l = q % hc;
            const int row = (j + hc) * width;
            if (west)
                c(s, row + l, iTime) = c(s, row + hc, iTime);
            if (east)
                c(s, row + hc + nx + l, iTime) = c(s, row + hc + nx - 1, iTime);
        });
    }

    /// Single-rank / per-direction replacement for fillPhysicalY + haloUnpack_y_rt.
    /// Covers the full RT width INCLUDING the X ghost columns so the corner
    /// blocks are propagated exactly as the X-then-Y buffer path does.
    inline void fillPhysicalYDirect(Domain &dom, RTStateSW &rt, int iTime)
    {
        const int hc = RTStateSW::rt_hc;
        const int width = rt.rt_nxhc;
        const int perSpec = hc * width;
        auto c = rt.c;
        const int ny = dom.ny, ns = n_spec;
        const bool north = dom.iN, south = dom.iS;
        Kokkos::parallel_for("RT_physical_y_direct", ns * perSpec, KOKKOS_LAMBDA(int idx) {
            const int s = idx / perSpec;
            const int q = idx % perSpec;
            const int l = q / width;
            const int i = q % width; // full width: includes X ghost columns (corners)
            if (north)
                c(s, l * width + i, iTime) = c(s, hc * width + i, iTime);
            if (south)
                c(s, (hc + ny + l) * width + i, iTime) = c(s, (hc + ny - 1) * width + i, iTime);
        });
    }

    inline void haloExchange_x_rt(Domain &dom, Parallel &par)
    {
        Kokkos::fence();
        if (par.nproc_x > 1) {
#if SERGHEI_RT_CUDA_AWARE_MPI
            // Requires a CUDA-aware MPI implementation.  The device/UVM
            // views stay alive until Waitall completes.
            timerStage.reset();
            SERGHEI_MPI_CHECK(MPI_Irecv(haloRecvBufW.data(),nBufX,SERGHEI_MPI_REAL,par.neigh(1,0),0,MPI_COMM_WORLD,&req[0]));
            SERGHEI_MPI_CHECK(MPI_Irecv(haloRecvBufE.data(),nBufX,SERGHEI_MPI_REAL,par.neigh(1,2),1,MPI_COMM_WORLD,&req[1]));
            SERGHEI_MPI_CHECK(MPI_Isend(haloSendBufW.data(),nBufX,SERGHEI_MPI_REAL,par.neigh(1,0),1,MPI_COMM_WORLD,&req[2]));
            SERGHEI_MPI_CHECK(MPI_Isend(haloSendBufE.data(),nBufX,SERGHEI_MPI_REAL,par.neigh(1,2),0,MPI_COMM_WORLD,&req[3]));
            SERGHEI_MPI_CHECK(MPI_Waitall(4, req, stat));
            dom.timers.rtsw.commMPIWait += timerStage.seconds();
#else
            timerStage.reset();
            Kokkos::deep_copy(h_sendW, haloSendBufW);
            Kokkos::deep_copy(h_sendE, haloSendBufE);
            dom.timers.rtsw.commHostCopy += timerStage.seconds();
            timerStage.reset();
            SERGHEI_MPI_CHECK(MPI_Irecv(h_recvW.data(),nBufX,SERGHEI_MPI_REAL,par.neigh(1,0),0,MPI_COMM_WORLD,&req[0]));
            SERGHEI_MPI_CHECK(MPI_Irecv(h_recvE.data(),nBufX,SERGHEI_MPI_REAL,par.neigh(1,2),1,MPI_COMM_WORLD,&req[1]));
            SERGHEI_MPI_CHECK(MPI_Isend(h_sendW.data(),nBufX,SERGHEI_MPI_REAL,par.neigh(1,0),1,MPI_COMM_WORLD,&req[2]));
            SERGHEI_MPI_CHECK(MPI_Isend(h_sendE.data(),nBufX,SERGHEI_MPI_REAL,par.neigh(1,2),0,MPI_COMM_WORLD,&req[3]));
            SERGHEI_MPI_CHECK(MPI_Waitall(4,req,stat));
            dom.timers.rtsw.commMPIWait += timerStage.seconds();
            timerStage.reset();
            Kokkos::deep_copy(haloRecvBufW, h_recvW);
            Kokkos::deep_copy(haloRecvBufE, h_recvE);
            dom.timers.rtsw.commHostCopy += timerStage.seconds();
#endif
        }
    }

    inline void haloExchange_y_rt(Domain &dom, Parallel &par)
    {
        // No fence here: the only caller (exchangeMPI) fences right after the
        // X-stage unpack / direct fill, which is the actual pack_y dependency.
        if (par.nproc_y > 1) {
#if SERGHEI_RT_CUDA_AWARE_MPI
            timerStage.reset();
            SERGHEI_MPI_CHECK(MPI_Irecv(haloRecvBufN.data(),nBufY,SERGHEI_MPI_REAL,par.neigh(0,1),0,MPI_COMM_WORLD,&req[0]));
            SERGHEI_MPI_CHECK(MPI_Irecv(haloRecvBufS.data(),nBufY,SERGHEI_MPI_REAL,par.neigh(2,1),1,MPI_COMM_WORLD,&req[1]));
            SERGHEI_MPI_CHECK(MPI_Isend(haloSendBufN.data(),nBufY,SERGHEI_MPI_REAL,par.neigh(0,1),1,MPI_COMM_WORLD,&req[2]));
            SERGHEI_MPI_CHECK(MPI_Isend(haloSendBufS.data(),nBufY,SERGHEI_MPI_REAL,par.neigh(2,1),0,MPI_COMM_WORLD,&req[3]));
            SERGHEI_MPI_CHECK(MPI_Waitall(4, req, stat));
            dom.timers.rtsw.commMPIWait += timerStage.seconds();
#else
            timerStage.reset();
            Kokkos::deep_copy(h_sendN, haloSendBufN);
            Kokkos::deep_copy(h_sendS, haloSendBufS);
            dom.timers.rtsw.commHostCopy += timerStage.seconds();
            timerStage.reset();
            SERGHEI_MPI_CHECK(MPI_Irecv(h_recvN.data(),nBufY,SERGHEI_MPI_REAL,par.neigh(0,1),0,MPI_COMM_WORLD,&req[0]));
            SERGHEI_MPI_CHECK(MPI_Irecv(h_recvS.data(),nBufY,SERGHEI_MPI_REAL,par.neigh(2,1),1,MPI_COMM_WORLD,&req[1]));
            SERGHEI_MPI_CHECK(MPI_Isend(h_sendN.data(),nBufY,SERGHEI_MPI_REAL,par.neigh(0,1),1,MPI_COMM_WORLD,&req[2]));
            SERGHEI_MPI_CHECK(MPI_Isend(h_sendS.data(),nBufY,SERGHEI_MPI_REAL,par.neigh(2,1),0,MPI_COMM_WORLD,&req[3]));
            SERGHEI_MPI_CHECK(MPI_Waitall(4,req,stat));
            dom.timers.rtsw.commMPIWait += timerStage.seconds();
            timerStage.reset();
            Kokkos::deep_copy(haloRecvBufN, h_recvN);
            Kokkos::deep_copy(haloRecvBufS, h_recvS);
            dom.timers.rtsw.commHostCopy += timerStage.seconds();
#endif
        }
    }

    inline void exchangeMPI(RTStateSW &rt, Domain &dom, Parallel &par, int iTime)
    {
        timerExchange.reset();
        dom.timers.rtsw.exchangeCount += 1.0;
        dom.timers.rtsw.exchangeBytes +=
            (par.nproc_x > 1 ? 2.0 * nBufX * sizeof(real) : 0.0) +
            (par.nproc_y > 1 ? 2.0 * nBufY * sizeof(real) : 0.0);

        // ---- X direction ----
        if (par.nproc_x > 1) {
            timerStage.reset();
            haloPack_x_rt(dom, rt, iTime);
            Kokkos::fence();
            dom.timers.rtsw.commPack += timerStage.seconds();
            timerStage.reset();
            haloExchange_x_rt(dom, par);     // MPI + Waitall
            timerStage.reset();
            fillPhysicalX(dom, rt, iTime);   // physical BC overrides periodic-neighbor recv:
                                             // must stay AFTER Waitall, BEFORE unpack
            haloUnpack_x_rt(dom, rt, iTime);
            dom.timers.rtsw.commUnpack += timerStage.seconds();
        } else {
            // No MPI neighbour in x: pack -> buffer -> unpack would be a pure
            // local round trip.  Write the physical ghosts directly instead.
            fillPhysicalXDirect(dom, rt, iTime);
        }

        // pack_y / fillPhysicalYDirect read the FULL width, i.e. the X ghost
        // (corner) columns just written above -> a barrier is required between
        // the X and Y stages on asynchronous backends (cf. Exchange.h:567).
        Kokkos::fence();

        // ---- Y direction ----
        if (par.nproc_y > 1) {
            timerStage.reset();
            haloPack_y_rt(dom, rt, iTime);
            Kokkos::fence();
            dom.timers.rtsw.commPack += timerStage.seconds();
            timerStage.reset();
            haloExchange_y_rt(dom, par);
            timerStage.reset();
            fillPhysicalY(dom, rt, iTime);
            haloUnpack_y_rt(dom, rt, iTime);
            dom.timers.rtsw.commUnpack += timerStage.seconds();
        } else {
            fillPhysicalYDirect(dom, rt, iTime);
        }

        Kokkos::fence();
        exchangeTime += timerExchange.seconds();
    }
};

#endif
#endif
