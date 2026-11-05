#pragma once

#include <fstream>
#include <iostream>
#include <cmath>
#include <type_traits>
#include "Parallel.h"
#include "MPIUtils.h"
#include <mpi.h>

// ============================================
// Granular timer sub-structs for crop, RTGW, RTSW
// These integrate Code1's detailed timing with Code2's structured timer approach
// ============================================

#if CROP_GROWTH_MODEL
struct CropTimers {
    static constexpr unsigned short n = 4;
    double total = 0;
    double rzm = 0;    // root zone moisture computation
    double calc = 0;   // crop growth calculation
    double sync = 0;   // LAI/RD data synchronization
};
#endif

#if SERGHEI_SUBSURFACE_TRANSPORT
struct RTGwTimers {
    static constexpr unsigned short n = 17;
    double total = 0;
    double prepare = 0;  // water-content synchronisation and midpoint preparation
    double bc = 0;       // boundary condition application
    double boundaryAudit = 0; // boundary mass-flux diagnostics
    double comm = 0;     // MPI communication
    double disp = 0;     // dispersion tensor computation
    double linsys = 0;   // linear system construction
    double linsol = 0;   // linear solver
    double reactKernel = 0; // chemical reaction integration kernel
    double reactAudit = 0;  // pre/post-reaction mass diagnostics
    double solidUpdate = 0; // solid-phase update after reaction
    double update = 0;   // concentration & solid phase update
    double convergence = 0; // Picard convergence reductions
    double massAudit = 0;   // transport mass-balance diagnostics
    double postprocess = 0; // final state synchronisation
    double other = 0;       // unclassified time inside total
    double timeseriesIntegrate = 0; // outside solver total
};
#endif

#if SERGHEI_SURFACE_TRANSPORT
struct RTSwTimers {
    static constexpr unsigned short n = 20;
    double total = 0;
    double copyExch = 0;   // copy subsurface exchange to surface transport
    double gwExch = 0;     // surface-subsurface solute exchange time
    double transportKernel = 0; // physical transport kernel
    double transportLaunch = 0; // host launch/enqueue time before device fence
    double transportFence = 0;  // explicit device synchronization time
    double transportAudit = 0;  // applied-mass audit and its synchronization
    double boundaryAudit = 0;   // boundary mass bookkeeping
    double reactKernel = 0;     // chemical reaction kernel
    double reactAudit = 0;      // pre/post-reaction mass diagnostics
    double fluxAudit = 0;       // transport mass-balance diagnostics
    double comm = 0;       // MPI communication time
    double commPack = 0;   // halo pack kernel and its existing synchronization
    double commHostCopy = 0; // host/device mirror copies on non-CUDA-aware MPI
    double commMPIWait = 0;  // MPI initiation/wait portion
    double commUnpack = 0; // halo unpack and physical-boundary kernels
    double exchangeCount = 0; // number of RTSW exchangeMPI calls
    double exchangeBytes = 0; // payload bytes issued by RTSW exchanges
    double other = 0;      // unclassified time inside total
    double timeseriesIntegrate = 0; // outside solver total
};
#endif

#if SERGHEI_SUBSURFACE_TRANSPORT
static_assert(std::is_standard_layout<RTGwTimers>::value,
              "RTGwTimers must remain a flat MPI-gatherable structure");
static_assert(sizeof(RTGwTimers) == RTGwTimers::n * sizeof(double),
              "RTGwTimers::n must match the number of double fields");
#endif
#if SERGHEI_SURFACE_TRANSPORT
static_assert(std::is_standard_layout<RTSwTimers>::value,
              "RTSwTimers must remain a flat MPI-gatherable structure");
static_assert(sizeof(RTSwTimers) == RTSwTimers::n * sizeof(double),
              "RTSwTimers::n must match the number of double fields");
#endif

class SergheiTimers {
    struct init_t {
        static constexpr unsigned short n = 1;
        double total = 0.0;
    };

    struct io_t {
        static constexpr unsigned short n = 3;
        double in = 0.0;
        double out = 0.0;
        double mpi = 0.0;
    };

    struct ss_t {
        static constexpr unsigned short n = 1;
        double raininf = 0.0;
    };

    struct flux_t {
        static constexpr unsigned short n = 3;
        double total = 0.0;
        double dt = 0.0;
        double mpi = 0.0;
    };

    struct update_t {
        static constexpr unsigned short n = 1;
        double total = 0.0;
    };
    
    struct wdc_t {
        static constexpr unsigned short n = 1;
        double total = 0.0;
    };

    struct halo_t {
        static constexpr unsigned short n = 2;
        double total = 0.0;
        double mpi = 0.0;
    };

    struct bc_t {
        static constexpr unsigned short n = 3;
        double total = 0.0;
        double integrate = 0.0;
        double mpi = 0.0;
    };

    struct dt_t {
        static constexpr unsigned short n = 2;
        double total = 0.0;
        double mpi = 0.0;
    };

    struct integrate_t {
        static constexpr unsigned short n = 2;
        double total = 0.0;
        double mpi = 0.0;
    };

    struct coupling_t {
        static constexpr unsigned short n = 1;
        double total = 0.0;
    };

    struct io_loop_t {
        static constexpr unsigned short n = 2;
        double fields = 0.0;
        double timeseries = 0.0;
    };

    struct phase_t {
        static constexpr unsigned short n = 4;
        double init = 0.0;
        double compute = 0.0;
        double finalise = 0.0;
        double other = 0.0;
    };

    struct swe_t {
        // WARNING: make sure to account for the sizes of the nested structs
        static constexpr unsigned short n = 4 + init_t::n + io_t::n + ss_t::n + 
                                         flux_t::n + update_t::n + wdc_t::n + halo_t::n + 
                                         integrate_t::n + bc_t::n + dt_t::n;
        double total = 0.0;
        double solve = 0.0;
        double other = 0.0;
        double mpi = 0.0;
        init_t init;
        io_t io;
        flux_t flux;
        update_t update;
        wdc_t wetdrycorr;
        halo_t halo;
        dt_t dt;
        bc_t bc;
        ss_t ss;
        integrate_t integrate;
    };
    
    #if SERGHEI_LPT || SERGHEI_LPT_RK || SERGHEI_LPT_RK_OFFLINE
    struct lpt_t {
        static constexpr unsigned short n = 3;
        double update = 0;
        double init   = 0;
        double out    = 0;
    };
    #endif

    #if SERGHEI_RE_MODEL
    struct re_t {
        static constexpr unsigned short n = 10;
        double gw           = 0;
        double gwBC         = 0;
        double gwlinsys     = 0;
        double gwlinsol     = 0;
        double gwUpdateK    = 0;
        double gwUpdateQ    = 0;
        double gwUpdateWC   = 0;
        double gwMPI        = 0;
        double gwIntegrate  = 0;
        double out          = 0;     
    };
    #endif

    #if SERGHEI_SUSPENDED_SEDIMENT
    struct st_t {
        static constexpr unsigned short n = 3;
        double total = 0.0;
        double halo = 0.0;
        double mpi = 0.0;
    };
    #endif

    #if SERGHEI_SCALAR_TRANSPORT
    struct ade_t {
        static constexpr unsigned short n = 3;
        double total = 0.0;
        double halo = 0.0;
        double mpi = 0.0;
    };
    #endif

public:
    double total = 0.0;
    phase_t phase;
    coupling_t coupling;
    io_loop_t io_loop;
    swe_t swe;
    #if SERGHEI_LPT || SERGHEI_LPT_RK || SERGHEI_LPT_RK_OFFLINE
    lpt_t lpt;
    #endif
    #if SERGHEI_RE_MODEL
    re_t re;
    #endif
    #if SERGHEI_SCALAR_TRANSPORT
    ade_t ade;
    #endif
    #if SERGHEI_SUSPENDED_SEDIMENT
    st_t st;
    #endif

    // Granular timers for crop, subsurface RT, surface RT modules
    #if CROP_GROWTH_MODEL
    CropTimers crop;
    #endif
    #if SERGHEI_SUBSURFACE_TRANSPORT
    RTGwTimers rtgw;
    #endif
    #if SERGHEI_SURFACE_TRANSPORT
    RTSwTimers rtsw;
    #endif

private:
    coupling_t* g_coupling = nullptr;
    phase_t* g_phase = nullptr;
    io_loop_t* g_io_loop = nullptr;
    swe_t* g_swe = nullptr;

    #if SERGHEI_LPT || SERGHEI_LPT_RK || SERGHEI_LPT_RK_OFFLINE
    lpt_t* g_lpt = nullptr;
    #endif

    #if SERGHEI_RE_MODEL
    re_t* g_re = nullptr;
    #endif

    #if SERGHEI_SUBSURFACE_TRANSPORT
    RTGwTimers* g_rtgw = nullptr;
    #endif
    #if SERGHEI_SURFACE_TRANSPORT
    RTSwTimers* g_rtsw = nullptr;
    #endif
    #if CROP_GROWTH_MODEL
    CropTimers* g_crop = nullptr;
    #endif
    double* g_total = nullptr;

public:
    SergheiTimers() = default;

    void closure() {
        swe.bc.total += swe.bc.integrate;
        swe.other = swe.solve - (swe.flux.total + swe.update.total + swe.wetdrycorr.total + swe.dt.total + 
                  swe.bc.total + swe.halo.total + swe.ss.raininf);
        swe.total = swe.solve + swe.init.total + swe.io.in + swe.io.out + 
                    swe.integrate.total;
        swe.mpi = swe.io.mpi + swe.flux.mpi + swe.halo.mpi + swe.dt.mpi + swe.bc.mpi + swe.integrate.mpi;
#if SERGHEI_SUBSURFACE_TRANSPORT
        rtgw.other = rtgw.total - (rtgw.prepare + rtgw.bc + rtgw.boundaryAudit + rtgw.comm +
                     rtgw.disp + rtgw.linsys + rtgw.linsol + rtgw.reactKernel +
                     rtgw.reactAudit + rtgw.solidUpdate + rtgw.update + rtgw.convergence +
                     rtgw.massAudit + rtgw.postprocess);
#endif
#if SERGHEI_SURFACE_TRANSPORT
        // transportLaunch/Fence/Audit are a diagnostic breakdown of the
        // legacy transportKernel total and must not be subtracted twice.
        // Likewise commPack/commHostCopy/commMPIWait/commUnpack accumulate
        // stages inside exchangeMPI and are already contained in rtsw.comm.
        rtsw.other = rtsw.total - (rtsw.copyExch + rtsw.gwExch + rtsw.transportKernel +
                     rtsw.boundaryAudit + rtsw.reactKernel + rtsw.reactAudit +
                     rtsw.fluxAudit + rtsw.comm);
#endif
        phase.other = total - (phase.init + phase.compute + phase.finalise);
    }

    void gather(Parallel const &par) {
        if(par.masterproc){
            g_total = new double[par.nranks];
            g_phase = new phase_t[par.nranks];
            g_coupling = new coupling_t[par.nranks];
            g_io_loop = new io_loop_t[par.nranks];
            g_swe = new swe_t[par.nranks];
        }

        // Gather serghei.total
        SERGHEI_MPI_CHECK(MPI_Gather(&total, 1, MPI_DOUBLE,
                   par.masterproc ? g_total : nullptr, 1, MPI_DOUBLE,
                   0, MPI_COMM_WORLD));

        SERGHEI_MPI_CHECK(MPI_Gather(reinterpret_cast<double*>(&phase), phase_t::n, MPI_DOUBLE,
                   par.masterproc ? reinterpret_cast<double*>(g_phase) : nullptr,
                   phase_t::n, MPI_DOUBLE, 0, MPI_COMM_WORLD));

        // Gather coupling timers
        SERGHEI_MPI_CHECK(MPI_Gather(reinterpret_cast<double*>(&coupling), coupling_t::n, MPI_DOUBLE,
                   par.masterproc ? reinterpret_cast<double*>(g_coupling) : nullptr,
                   coupling_t::n, MPI_DOUBLE, 0, MPI_COMM_WORLD));

        // Gather io_loop timers
        SERGHEI_MPI_CHECK(MPI_Gather(reinterpret_cast<double*>(&io_loop), io_loop_t::n, MPI_DOUBLE,
                   par.masterproc ? reinterpret_cast<double*>(g_io_loop) : nullptr,
                   io_loop_t::n, MPI_DOUBLE, 0, MPI_COMM_WORLD));

        // Gather all timer values as doubles
        SERGHEI_MPI_CHECK(MPI_Gather(
            reinterpret_cast<double*>(&swe),
            swe_t::n,
            MPI_DOUBLE,
            par.masterproc ? reinterpret_cast<double*>(g_swe) : nullptr,
            swe_t::n,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
        ));

        #if SERGHEI_LPT || SERGHEI_LPT_RK || SERGHEI_LPT_RK_OFFLINE
        if(par.masterproc) g_lpt = new lpt_t[par.nranks];
        SERGHEI_MPI_CHECK(MPI_Gather(&lpt, lpt_t::n, MPI_DOUBLE,
                   par.masterproc ? g_lpt : nullptr, lpt_t::n, MPI_DOUBLE,
                   0, MPI_COMM_WORLD));
        #endif

        #if SERGHEI_RE_MODEL
        if(par.masterproc) g_re = new re_t[par.nranks];
        SERGHEI_MPI_CHECK(MPI_Gather(&re, re_t::n, MPI_DOUBLE,
                   par.masterproc ? g_re : nullptr, re_t::n, MPI_DOUBLE,
                   0, MPI_COMM_WORLD));
        #endif

        #if SERGHEI_SUBSURFACE_TRANSPORT
        if(par.masterproc) g_rtgw = new RTGwTimers[par.nranks];
        SERGHEI_MPI_CHECK(MPI_Gather(&rtgw, RTGwTimers::n, MPI_DOUBLE,
                   par.masterproc ? g_rtgw : nullptr, RTGwTimers::n, MPI_DOUBLE,
                   0, MPI_COMM_WORLD));
        #endif

        #if SERGHEI_SURFACE_TRANSPORT
        if(par.masterproc) g_rtsw = new RTSwTimers[par.nranks];
        SERGHEI_MPI_CHECK(MPI_Gather(&rtsw, RTSwTimers::n, MPI_DOUBLE,
                   par.masterproc ? g_rtsw : nullptr, RTSwTimers::n, MPI_DOUBLE,
                   0, MPI_COMM_WORLD));
        #endif

        #if CROP_GROWTH_MODEL
        if(par.masterproc) g_crop = new CropTimers[par.nranks];
        SERGHEI_MPI_CHECK(MPI_Gather(&crop, CropTimers::n, MPI_DOUBLE,
                   par.masterproc ? g_crop : nullptr, CropTimers::n, MPI_DOUBLE,
                   0, MPI_COMM_WORLD));
        #endif
    }

    void computeRelative(const SergheiTimers &base){
        total /= base.total;
        phase.init /= base.total;
        phase.compute /= base.total;
        phase.finalise /= base.total;
        phase.other /= base.total;
        coupling.total /= base.total;
        io_loop.fields /= base.total;
        io_loop.timeseries /= base.total;
        swe.total /= base.total;
        swe.solve /=  base.total;
        swe.init.total /= base.total;
        swe.io.in /= base.total;
        swe.io.out /= base.total; 
        swe.io.mpi /= base.total; 
        swe.flux.total /= base.total;
        swe.flux.dt /= base.total;
        swe.flux.mpi /= base.total;
        swe.update.total /= base.total;
        swe.wetdrycorr.total /= base.total;
        swe.halo.total /= base.total;
        swe.halo.mpi /= base.total;
        swe.dt.total /= base.total;
        swe.dt.mpi /= base.total;
        swe.bc.total /= base.total;
        swe.bc.integrate /= base.total;
        swe.bc.mpi /= base.total;
        swe.ss.raininf /= base.total;
        swe.other /= base.total;
        swe.integrate.total /= base.total;
        swe.integrate.mpi /= base.total; 
        swe.mpi /= base.total;

        #if SERGHEI_LPT || SERGHEI_LPT_RK || SERGHEI_LPT_RK_OFFLINE
        lpt.update /= base.total;
        lpt.init /= base.total;
        lpt.out /= base.total;
        #endif

        #if SERGHEI_RE_MODEL
        re.gw /= base.total;
        re.gwBC /= base.total;
        re.gwlinsys /= base.total;
        re.gwlinsol /= base.total;
        re.gwUpdateK /= base.total;
        re.gwUpdateQ /= base.total;
        re.gwUpdateWC /= base.total;
        re.gwMPI /= base.total;
        re.gwIntegrate /= base.total;
        re.out /= base.total;
        #endif

        #if SERGHEI_SUBSURFACE_TRANSPORT
        rtgw.total /= base.total;
        rtgw.prepare /= base.total;
        rtgw.bc /= base.total;
        rtgw.boundaryAudit /= base.total;
        rtgw.comm /= base.total;
        rtgw.disp /= base.total;
        rtgw.linsys /= base.total;
        rtgw.linsol /= base.total;
        rtgw.reactKernel /= base.total;
        rtgw.reactAudit /= base.total;
        rtgw.solidUpdate /= base.total;
        rtgw.update /= base.total;
        rtgw.convergence /= base.total;
        rtgw.massAudit /= base.total;
        rtgw.postprocess /= base.total;
        rtgw.other /= base.total;
        rtgw.timeseriesIntegrate /= base.total;
        #endif

        #if SERGHEI_SURFACE_TRANSPORT
        rtsw.total /= base.total;
        rtsw.copyExch /= base.total;
        rtsw.gwExch /= base.total;
        rtsw.transportKernel /= base.total;
        rtsw.transportLaunch /= base.total;
        rtsw.transportFence /= base.total;
        rtsw.transportAudit /= base.total;
        rtsw.boundaryAudit /= base.total;
        rtsw.reactKernel /= base.total;
        rtsw.reactAudit /= base.total;
        rtsw.fluxAudit /= base.total;
        rtsw.comm /= base.total;
        rtsw.commPack /= base.total;
        rtsw.commHostCopy /= base.total;
        rtsw.commMPIWait /= base.total;
        rtsw.commUnpack /= base.total;
        rtsw.other /= base.total;
        rtsw.timeseriesIntegrate /= base.total;
        #endif

        #if CROP_GROWTH_MODEL
        crop.total /= base.total;
        crop.rzm /= base.total;
        crop.calc /= base.total;
        crop.sync /= base.total;
        #endif
    }

    inline std::string indent(std::size_t n){
        return std::string(n, ' ');
    }
    inline std::string level(std::size_t n){
        return indent(2*n);
    }   

    void report(std::ofstream &file, int nrank){

        file << "serghei.total";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_total[r];
        file << std::endl;

        file << level(1) << "phase.init";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_phase[r].init;
        file << std::endl;
        file << level(1) << "phase.compute";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_phase[r].compute;
        file << std::endl;
        file << level(1) << "phase.finalise";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_phase[r].finalise;
        file << std::endl;
        file << level(1) << "phase.other";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_phase[r].other;
        file << std::endl;

        file << level(1) << "coupling.total";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_coupling[r].total;
        file << std::endl;

        file << level(1) << "io_loop.fields";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_io_loop[r].fields;
        file << std::endl;
        file << level(1) << "io_loop.timeseries";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_io_loop[r].timeseries;
        file << std::endl;

        file << level(1) << "swe.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].total;
        file << std::endl;

        file << level(2) << "swe.init.total"; 
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].init.total;
        file << std::endl;

        file << level(2) <<"swe.io.in";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].io.in;
        file << std::endl;

        file << level(2) << "swe.io.out";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].io.out;
        file << std::endl;
        
        file << level(2) << "swe.io.mpi";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].io.mpi;
        file << std::endl;
        
        file << level(2) << "swe.solve";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].solve;
        file << std::endl;

        file << level(3) << "swe.flux.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].flux.total;
        file << std::endl;
        
        file << level(4) << "swe.flux.dt";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].flux.dt;
        file << std::endl;

        file << level(4) << "swe.flux.mpi";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].flux.mpi;
        file << std::endl;

        file << level(3) << "swe.update.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].update.total;
        file << std::endl;

        file << level(3) << "swe.wetdrycorr.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].wetdrycorr.total;
        file << std::endl;

        file << level(3) << "swe.halo.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].halo.total;
        file << std::endl;

        file << level(4) << "swe.halo.mpi";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].halo.mpi;
        file << std::endl;

        file << level(3) << "swe.dt.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].dt.total;
        file << std::endl;

        file << level(4) << "swe.dt.mpi";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].dt.mpi;
        file << std::endl;

        file << level(3) << "swe.bc.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].bc.total;
        file << std::endl;

        file << level(4) << "swe.bc.integrate";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].bc.integrate;
        file << std::endl;

        file << level(4) << "swe.bc.mpi";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].bc.mpi;
        file << std::endl;

        file << level(3) << "swe.ss.raininf";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].ss.raininf;
        file << std::endl;

        file << level(3) << "swe.other";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].other;
        file << std::endl;

        file << level(2) << "swe.integrate.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].integrate.total;
        file << std::endl;

        file << level(3) << "swe.integrate.mpi";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].integrate.mpi;
        file << std::endl;
        
        file << level(2) << "(swe.mpi)";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].mpi;
        file << std::endl;


        #if SERGHEI_LPT || SERGHEI_LPT_RK || SERGHEI_LPT_RK_OFFLINE
        // Report LPT timers
        file << level(1) << "lpt.update";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_lpt[r].update;
        file << std::endl;

        file << level(1) << "lpt.init";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_lpt[r].init;
        file << std::endl;

        file << level(1) << "lpt.out";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_lpt[r].out;
        file << std::endl;
        #endif

        #if SERGHEI_RE_MODEL
        file << level(1) << "re.gw";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gw;
        file << std::endl;

        file << level(2) << "re.gwBC";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gwBC;
        file << std::endl;

        file << level(2) << "re.gwlinsys";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gwlinsys;
        file << std::endl;

        file << level(2) << "re.gwlinsol";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gwlinsol;
        file << std::endl;

        file << level(2) << "re.gwUpdateK";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gwUpdateK;
        file << std::endl;

        file << level(2) << "re.gwUpdateQ";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gwUpdateQ;
        file << std::endl;

        file << level(2) << "re.gwUpdateWC";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gwUpdateWC;
        file << std::endl;

        file << level(2) << "re.gwMPI";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gwMPI;
        file << std::endl;

        file << level(2) << "re.gwIntegrate";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gwIntegrate;
        file << std::endl;

        file << level(1) << "re.out";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].out;
        file << std::endl;
        #endif

        #if SERGHEI_SUBSURFACE_TRANSPORT
        file << level(1) << "rtgw.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].total;
        file << std::endl;

        file << level(2) << "rtgw.prepare";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].prepare;
        file << std::endl;
        file << level(2) << "rtgw.bc";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].bc;
        file << std::endl;
        file << level(2) << "rtgw.boundaryAudit";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].boundaryAudit;
        file << std::endl;

        file << level(2) << "rtgw.comm";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].comm;
        file << std::endl;

        file << level(2) << "rtgw.disp";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].disp;
        file << std::endl;

        file << level(2) << "rtgw.linsys";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].linsys;
        file << std::endl;

        file << level(2) << "rtgw.linsol";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].linsol;
        file << std::endl;

        file << level(2) << "rtgw.reactKernel";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].reactKernel;
        file << std::endl;
        file << level(2) << "rtgw.reactAudit";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].reactAudit;
        file << std::endl;
        file << level(2) << "rtgw.solidUpdate";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].solidUpdate;
        file << std::endl;

        file << level(2) << "rtgw.update";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].update;
        file << std::endl;

        file << level(2) << "rtgw.convergence";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].convergence;
        file << std::endl;
        file << level(2) << "rtgw.massAudit";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].massAudit;
        file << std::endl;
        file << level(2) << "rtgw.postprocess";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].postprocess;
        file << std::endl;
        file << level(2) << "rtgw.other";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].other;
        file << std::endl;
        file << level(1) << "rtgw.timeseriesIntegrate";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].timeseriesIntegrate;
        file << std::endl;
        #endif

        #if SERGHEI_SURFACE_TRANSPORT
        file << level(1) << "rtsw.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].total;
        file << std::endl;

        file << level(2) << "rtsw.copyExch";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].copyExch;
        file << std::endl;
        file << level(2) << "rtsw.gwExch";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].gwExch;
        file << std::endl;
        file << level(2) << "rtsw.transportKernel";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].transportKernel;
        file << std::endl;
        file << level(2) << "rtsw.transportLaunch";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].transportLaunch;
        file << std::endl;
        file << level(2) << "rtsw.transportFence";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].transportFence;
        file << std::endl;
        file << level(2) << "rtsw.transportAudit";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].transportAudit;
        file << std::endl;
        file << level(2) << "rtsw.boundaryAudit";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].boundaryAudit;
        file << std::endl;
        file << level(2) << "rtsw.reactKernel";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].reactKernel;
        file << std::endl;
        file << level(2) << "rtsw.reactAudit";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].reactAudit;
        file << std::endl;
        file << level(2) << "rtsw.fluxAudit";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].fluxAudit;
        file << std::endl;

        file << level(2) << "rtsw.comm";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].comm;
        file << std::endl;
        file << level(3) << "rtsw.commPack";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].commPack;
        file << std::endl;
        file << level(3) << "rtsw.commHostCopy";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].commHostCopy;
        file << std::endl;
        file << level(3) << "rtsw.commMPIWait";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].commMPIWait;
        file << std::endl;
        file << level(3) << "rtsw.commUnpack";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].commUnpack;
        file << std::endl;
        file << level(3) << "rtsw.exchangeCount";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].exchangeCount;
        file << std::endl;
        file << level(3) << "rtsw.exchangeBytes";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].exchangeBytes;
        file << std::endl;
        file << level(2) << "rtsw.other";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].other;
        file << std::endl;
        file << level(1) << "rtsw.timeseriesIntegrate";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].timeseriesIntegrate;
        file << std::endl;
        #endif

        #if CROP_GROWTH_MODEL
        file << level(1) << "crop.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_crop[r].total;
        file << std::endl;

        file << level(2) << "crop.rzm";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_crop[r].rzm;
        file << std::endl;

        file << level(2) << "crop.calc";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_crop[r].calc;
        file << std::endl;

        file << level(2) << "crop.sync";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_crop[r].sync;
        file << std::endl;
        #endif
    }

    void reportClosure(std::ofstream &file, int nrank) {
        auto writeClosure = [&](const std::string &name, auto totalAt, auto otherAt) {
            file << name;
            for (int r = 0; r < nrank; ++r) {
                const double parent = totalAt(r);
                const double other = otherAt(r);
                const double rawRatio = parent > 0.0 ? 1.0 - std::abs(other) / parent : 1.0;
                const double ratio = rawRatio < 0.0 ? 0.0 : (rawRatio > 1.0 ? 1.0 : rawRatio);
                file << "\t" << ratio;
            }
            file << std::endl;
        };

        writeClosure("phase.closureRatio",
                     [this](int r) { return g_total[r]; },
                     [this](int r) { return g_phase[r].other; });
#if SERGHEI_SUBSURFACE_TRANSPORT
        writeClosure("rtgw.closureRatio",
                     [this](int r) { return g_rtgw[r].total; },
                     [this](int r) { return g_rtgw[r].other; });
#endif
#if SERGHEI_SURFACE_TRANSPORT
        writeClosure("rtsw.closureRatio",
                     [this](int r) { return g_rtsw[r].total; },
                     [this](int r) { return g_rtsw[r].other; });
#endif
    }

    double inline getAverageTimer(double mytimer, int nrank){
        double avg = 0.0;
        for(int r=0; r<nrank; ++r) avg += mytimer;
        if(avg == 0.0) return 1.0;
        return avg / static_cast<double>(nrank);
    }

    double inline getAverageTimer(std::function<double(int)> getTimer, int nrank){
        double avg = 0.0;
        for(int r = 0; r < nrank; ++r) avg += getTimer(r);
        if(avg == 0.0) return 1.0;
        return avg / static_cast<double>(nrank);
    }

    void reportLoadBalance(std::ofstream &file, int nrank){
        double avg = 0.0;

        avg = getAverageTimer([this](int r){ return g_total[r]; }, nrank);
        file << "serghei.total";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_total[r] / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_phase[r].init; }, nrank);
        file << level(1) << "phase.init";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_phase[r].init / avg;
        file << std::endl;
        avg = getAverageTimer([this](int r){ return g_phase[r].compute; }, nrank);
        file << level(1) << "phase.compute";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_phase[r].compute / avg;
        file << std::endl;
        avg = getAverageTimer([this](int r){ return g_phase[r].finalise; }, nrank);
        file << level(1) << "phase.finalise";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_phase[r].finalise / avg;
        file << std::endl;
        avg = getAverageTimer([this](int r){ return g_phase[r].other; }, nrank);
        file << level(1) << "phase.other";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_phase[r].other / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_coupling[r].total; }, nrank);
        file << level(1) << "coupling.total";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_coupling[r].total / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_io_loop[r].fields; }, nrank);
        file << level(1) << "io_loop.fields";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_io_loop[r].fields / avg;
        file << std::endl;
        avg = getAverageTimer([this](int r){ return g_io_loop[r].timeseries; }, nrank);
        file << level(1) << "io_loop.timeseries";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_io_loop[r].timeseries / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].total; }, nrank);
        file << level(1) << "swe.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].total / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].init.total; }, nrank);
        file << level(2) << "swe.init.total"; 
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].init.total / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].io.in; }, nrank);
        file << level(2) <<"swe.io.in";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].io.in / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].io.out; }, nrank);
        file << level(2) << "swe.io.out";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].io.out / avg;
        file << std::endl;
        
        avg = getAverageTimer([this](int r){ return g_swe[r].io.mpi; }, nrank);
        file << level(2) << "swe.io.mpi";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].io.mpi / avg;
        file << std::endl;
        
        avg = getAverageTimer([this](int r){ return g_swe[r].solve; }, nrank);
        file << level(2) << "swe.solve";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].solve / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].flux.total; }, nrank);
        file << level(3) << "swe.flux.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].flux.total / avg;
        file << std::endl;
        
        avg = getAverageTimer([this](int r){ return g_swe[r].flux.dt; }, nrank);
        file << level(4) << "swe.flux.dt";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].flux.dt / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].flux.mpi; }, nrank);
        file << level(4) << "swe.flux.mpi";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].flux.mpi / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].update.total; }, nrank);
        file << level(3) << "swe.update.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].update.total / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].wetdrycorr.total; }, nrank);
        file << level(3) << "swe.wetdrycorr.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].wetdrycorr.total / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].halo.total; }, nrank);
        file << level(3) << "swe.halo.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].halo.total / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].halo.mpi; }, nrank);
        file << level(4) << "swe.halo.mpi";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].halo.mpi / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].dt.total; }, nrank);
        file << level(3) << "swe.dt.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].dt.total / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].dt.mpi; }, nrank);
        file << level(4) << "swe.dt.mpi";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].dt.mpi / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].bc.total; }, nrank);
        file << level(3) << "swe.bc.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].bc.total / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].bc.integrate; }, nrank);
        file << level(4) << "swe.bc.integrate";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].bc.integrate / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].bc.mpi; }, nrank);
        file << level(4) << "swe.bc.mpi";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].bc.mpi / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].ss.raininf; }, nrank);
        file << level(3) << "swe.ss.raininf";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].ss.raininf / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].other; }, nrank);
        file << level(3) << "swe.other";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].other / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].integrate.total; }, nrank);
        file << level(2) << "swe.integrate.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].integrate.total / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_swe[r].integrate.mpi; }, nrank);
        file << level(3) << "swe.integrate.mpi";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].integrate.mpi / avg;
        file << std::endl;
        
        avg = getAverageTimer([this](int r){ return g_swe[r].mpi; }, nrank);
        file << level(2) << "(swe.mpi)";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_swe[r].mpi / avg;
        file << std::endl;

        #if SERGHEI_LPT || SERGHEI_LPT_RK || SERGHEI_LPT_RK_OFFLINE
        avg = getAverageTimer([this](int r){ return g_lpt[r].update; }, nrank);
        file << level(1) << "lpt.update";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_lpt[r].update / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_lpt[r].init; }, nrank);
        file << level(1) << "lpt.init";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_lpt[r].init / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_lpt[r].out; }, nrank);
        file << level(1) << "lpt.out";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_lpt[r].out / avg;
        file << std::endl;
        #endif

        #if SERGHEI_RE_MODEL
        avg = getAverageTimer([this](int r){ return g_re[r].gw; }, nrank);
        file << level(1) << "re.gw";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gw / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_re[r].gwBC; }, nrank);
        file << level(2) << "re.gwBC";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gwBC / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_re[r].gwlinsys; }, nrank);
        file << level(2) << "re.gwlinsys";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gwlinsys / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_re[r].gwlinsol; }, nrank);
        file << level(2) << "re.gwlinsol";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gwlinsol / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_re[r].gwUpdateK; }, nrank);
        file << level(2) << "re.gwUpdateK";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gwUpdateK / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_re[r].gwUpdateQ; }, nrank);
        file << level(2) << "re.gwUpdateQ";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gwUpdateQ / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_re[r].gwUpdateWC; }, nrank);
        file << level(2) << "re.gwUpdateWC";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gwUpdateWC / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_re[r].gwMPI; }, nrank);
        file << level(2) << "re.gwMPI";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gwMPI / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_re[r].gwIntegrate; }, nrank);
        file << level(2) << "re.gwIntegrate";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].gwIntegrate / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_re[r].out; }, nrank);
        file << level(1) << "re.out";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_re[r].out / avg;
        file << std::endl;
        #endif

        #if SERGHEI_SUBSURFACE_TRANSPORT
        avg = getAverageTimer([this](int r){ return g_rtgw[r].total; }, nrank);
        file << level(1) << "rtgw.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].total / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtgw[r].prepare; }, nrank);
        file << level(2) << "rtgw.prepare";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].prepare / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtgw[r].bc; }, nrank);
        file << level(2) << "rtgw.bc";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].bc / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtgw[r].boundaryAudit; }, nrank);
        file << level(2) << "rtgw.boundaryAudit";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].boundaryAudit / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtgw[r].comm; }, nrank);
        file << level(2) << "rtgw.comm";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].comm / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtgw[r].disp; }, nrank);
        file << level(2) << "rtgw.disp";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].disp / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtgw[r].linsys; }, nrank);
        file << level(2) << "rtgw.linsys";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].linsys / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtgw[r].linsol; }, nrank);
        file << level(2) << "rtgw.linsol";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].linsol / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtgw[r].reactKernel; }, nrank);
        file << level(2) << "rtgw.reactKernel";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].reactKernel / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtgw[r].reactAudit; }, nrank);
        file << level(2) << "rtgw.reactAudit";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].reactAudit / avg;
        file << std::endl;
        avg = getAverageTimer([this](int r){ return g_rtgw[r].solidUpdate; }, nrank);
        file << level(2) << "rtgw.solidUpdate";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].solidUpdate / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtgw[r].update; }, nrank);
        file << level(2) << "rtgw.update";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].update / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtgw[r].convergence; }, nrank);
        file << level(2) << "rtgw.convergence";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].convergence / avg;
        file << std::endl;
        avg = getAverageTimer([this](int r){ return g_rtgw[r].massAudit; }, nrank);
        file << level(2) << "rtgw.massAudit";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].massAudit / avg;
        file << std::endl;
        avg = getAverageTimer([this](int r){ return g_rtgw[r].postprocess; }, nrank);
        file << level(2) << "rtgw.postprocess";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].postprocess / avg;
        file << std::endl;
        avg = getAverageTimer([this](int r){ return g_rtgw[r].other; }, nrank);
        file << level(2) << "rtgw.other";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].other / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtgw[r].timeseriesIntegrate; }, nrank);
        file << level(1) << "rtgw.timeseriesIntegrate";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].timeseriesIntegrate / avg;
        file << std::endl;
        #endif

        #if SERGHEI_SURFACE_TRANSPORT
        avg = getAverageTimer([this](int r){ return g_rtsw[r].total; }, nrank);
        file << level(1) << "rtsw.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].total / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtsw[r].copyExch; }, nrank);
        file << level(2) << "rtsw.copyExch";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].copyExch / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtsw[r].transportKernel; }, nrank);
        file << level(2) << "rtsw.transportKernel";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].transportKernel / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtsw[r].boundaryAudit; }, nrank);
        file << level(2) << "rtsw.boundaryAudit";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].boundaryAudit / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtsw[r].reactKernel; }, nrank);
        file << level(2) << "rtsw.reactKernel";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].reactKernel / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtsw[r].reactAudit; }, nrank);
        file << level(2) << "rtsw.reactAudit";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].reactAudit / avg;
        file << std::endl;
        avg = getAverageTimer([this](int r){ return g_rtsw[r].fluxAudit; }, nrank);
        file << level(2) << "rtsw.fluxAudit";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].fluxAudit / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtsw[r].timeseriesIntegrate; }, nrank);
        file << level(1) << "rtsw.timeseriesIntegrate";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].timeseriesIntegrate / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtsw[r].gwExch; }, nrank);
        file << level(2) << "rtsw.gwExch";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].gwExch / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtsw[r].comm; }, nrank);
        file << level(2) << "rtsw.comm";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].comm / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtsw[r].other; }, nrank);
        file << level(2) << "rtsw.other";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].other / avg;
        file << std::endl;
        #endif

        #if CROP_GROWTH_MODEL
        avg = getAverageTimer([this](int r){ return g_crop[r].total; }, nrank);
        file << level(1) << "crop.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_crop[r].total / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_crop[r].rzm; }, nrank);
        file << level(2) << "crop.rzm";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_crop[r].rzm / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_crop[r].calc; }, nrank);
        file << level(2) << "crop.calc";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_crop[r].calc / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_crop[r].sync; }, nrank);
        file << level(2) << "crop.sync";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_crop[r].sync / avg;
        file << std::endl;
        #endif
    }
};
