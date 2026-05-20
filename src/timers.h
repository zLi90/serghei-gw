#pragma once

#include <fstream>
#include <iostream>
#include "Parallel.h"
#include <mpi.h>

// ============================================
// [FROM CODE1] Granular timer sub-structs for crop, RTGW, RTSW
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
    static constexpr unsigned short n = 9;
    double total = 0;
    double bc = 0;       // boundary condition application
    double comm = 0;     // MPI communication
    double disp = 0;     // dispersion tensor computation
    double linsys = 0;   // linear system construction
    double linsol = 0;   // linear solver
    double react = 0;    // chemical reaction integration
    double update = 0;   // concentration & solid phase update
    double integrate = 0; // integrator total time
};
#endif

#if SERGHEI_SURFACE_TRANSPORT
struct RTSwTimers {
    static constexpr unsigned short n = 6;
    double total = 0;
    double transport = 0;  // physical transport
    double react = 0;      // chemical reaction integration
    double integrate = 0;  // integrator
    double gwExch = 0;     // surface-subsurface solute exchange time
    double comm = 0;       // MPI communication time
};
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
        static constexpr unsigned short n = 1;
        double output = 0.0;
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

    // [FROM CODE1] Granular timers for crop, subsurface RT, surface RT modules
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
    }

    void gather(Parallel const &par) {
        if(par.masterproc){
            g_total = new double[par.nranks];
            g_coupling = new coupling_t[par.nranks];
            g_io_loop = new io_loop_t[par.nranks];
            g_swe = new swe_t[par.nranks];
        }

        // Gather serghei.total
        MPI_Gather(&total, 1, MPI_DOUBLE,
                   par.masterproc ? g_total : nullptr, 1, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);

        // Gather coupling timers
        MPI_Gather(reinterpret_cast<double*>(&coupling), coupling_t::n, MPI_DOUBLE,
                   par.masterproc ? reinterpret_cast<double*>(g_coupling) : nullptr,
                   coupling_t::n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Gather io_loop timers
        MPI_Gather(reinterpret_cast<double*>(&io_loop), io_loop_t::n, MPI_DOUBLE,
                   par.masterproc ? reinterpret_cast<double*>(g_io_loop) : nullptr,
                   io_loop_t::n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Gather all timer values as doubles
        MPI_Gather(
            reinterpret_cast<double*>(&swe),
            swe_t::n,
            MPI_DOUBLE,
            par.masterproc ? reinterpret_cast<double*>(g_swe) : nullptr,
            swe_t::n,
            MPI_DOUBLE,
            0,
            MPI_COMM_WORLD
        );

        #if SERGHEI_LPT || SERGHEI_LPT_RK || SERGHEI_LPT_RK_OFFLINE
        if(par.masterproc) g_lpt = new lpt_t[par.nranks];
        MPI_Gather(&lpt, lpt_t::n, MPI_DOUBLE,
                   par.masterproc ? g_lpt : nullptr, lpt_t::n, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        #endif

        #if SERGHEI_RE_MODEL
        if(par.masterproc) g_re = new re_t[par.nranks];
        MPI_Gather(&re, re_t::n, MPI_DOUBLE,
                   par.masterproc ? g_re : nullptr, re_t::n, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        #endif

        #if SERGHEI_SUBSURFACE_TRANSPORT
        if(par.masterproc) g_rtgw = new RTGwTimers[par.nranks];
        MPI_Gather(&rtgw, RTGwTimers::n, MPI_DOUBLE,
                   par.masterproc ? g_rtgw : nullptr, RTGwTimers::n, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        #endif

        #if SERGHEI_SURFACE_TRANSPORT
        if(par.masterproc) g_rtsw = new RTSwTimers[par.nranks];
        MPI_Gather(&rtsw, RTSwTimers::n, MPI_DOUBLE,
                   par.masterproc ? g_rtsw : nullptr, RTSwTimers::n, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        #endif

        #if CROP_GROWTH_MODEL
        if(par.masterproc) g_crop = new CropTimers[par.nranks];
        MPI_Gather(&crop, CropTimers::n, MPI_DOUBLE,
                   par.masterproc ? g_crop : nullptr, CropTimers::n, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);
        #endif
    }

    void computeRelative(const SergheiTimers &base){
        total /= base.total;
        coupling.total /= base.total;
        io_loop.output /= base.total;
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
        rtgw.bc /= base.total;
        rtgw.comm /= base.total;
        rtgw.disp /= base.total;
        rtgw.linsys /= base.total;
        rtgw.linsol /= base.total;
        rtgw.react /= base.total;
        rtgw.update /= base.total;
        rtgw.integrate /= base.total;
        #endif

        #if SERGHEI_SURFACE_TRANSPORT
        rtsw.total /= base.total;
        rtsw.transport /= base.total;
        rtsw.react /= base.total;
        rtsw.integrate /= base.total;
        rtsw.gwExch /= base.total;
        rtsw.comm /= base.total;
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
        for(int r = 0; r < nrank; ++r) file << "\t" << total;
        file << std::endl;

        file << level(1) << "coupling.total";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_coupling[r].total;
        file << std::endl;

        file << level(1) << "io_loop.output";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_io_loop[r].output;
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

        file << level(2) << "rtgw.bc";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].bc;
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

        file << level(2) << "rtgw.react";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].react;
        file << std::endl;

        file << level(2) << "rtgw.update";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].update;
        file << std::endl;

        file << level(2) << "rtgw.integrate";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].integrate;
        file << std::endl;
        #endif

        #if SERGHEI_SURFACE_TRANSPORT
        file << level(1) << "rtsw.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].total;
        file << std::endl;

        file << level(2) << "rtsw.transport";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].transport;
        file << std::endl;

        file << level(2) << "rtsw.react";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].react;
        file << std::endl;

        file << level(2) << "rtsw.integrate";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].integrate;
        file << std::endl;

        file << level(2) << "rtsw.gwExch";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].gwExch;
        file << std::endl;

        file << level(2) << "rtsw.comm";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].comm;
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

        avg = getAverageTimer([this](int r){ return g_coupling[r].total; }, nrank);
        file << level(1) << "coupling.total";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_coupling[r].total / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_io_loop[r].output; }, nrank);
        file << level(1) << "io_loop.output";
        for(int r = 0; r < nrank; ++r) file << "\t" << g_io_loop[r].output / avg;
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

        avg = getAverageTimer([this](int r){ return g_rtgw[r].bc; }, nrank);
        file << level(2) << "rtgw.bc";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].bc / avg;
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

        avg = getAverageTimer([this](int r){ return g_rtgw[r].react; }, nrank);
        file << level(2) << "rtgw.react";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].react / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtgw[r].update; }, nrank);
        file << level(2) << "rtgw.update";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].update / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtgw[r].integrate; }, nrank);
        file << level(2) << "rtgw.integrate";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtgw[r].integrate / avg;
        file << std::endl;
        #endif

        #if SERGHEI_SURFACE_TRANSPORT
        avg = getAverageTimer([this](int r){ return g_rtsw[r].total; }, nrank);
        file << level(1) << "rtsw.total";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].total / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtsw[r].transport; }, nrank);
        file << level(2) << "rtsw.transport";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].transport / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtsw[r].react; }, nrank);
        file << level(2) << "rtsw.react";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].react / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtsw[r].integrate; }, nrank);
        file << level(2) << "rtsw.integrate";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].integrate / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtsw[r].gwExch; }, nrank);
        file << level(2) << "rtsw.gwExch";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].gwExch / avg;
        file << std::endl;

        avg = getAverageTimer([this](int r){ return g_rtsw[r].comm; }, nrank);
        file << level(2) << "rtsw.comm";
        for (int r = 0; r < nrank; ++r) file << "\t" << g_rtsw[r].comm / avg;
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
