#ifndef _RT_FUNCTION_H_
#define _RT_FUNCTION_H_

#include "const.h"
#include "define.h"
#include "SArray.h"
#include "Indexing.h"
#include "GwDomain.h"
#include "GwMPI.h"
#include "GwState.h"
#include "State.h"
#include <set>
#include <math.h>


class RTFunction   {

private:
	Kokkos::Timer timer;

public:
    /* --------------------------------------------------
        Top-level RT solver
    -------------------------------------------------- */
    inline void rt_solve(RTState &rt, GwState &gw, GwDomain &gdom, GwMPI &gmpi, Parallel &par)  {
        real dt_tmp;
		
		/*
			APPLY TRANSPORT BOUNDARY CONDITIONS
		*/
		
		/*
			CALCULATE DISPERSION TENSOR (rt.d)
		*/
		
		/*
			SOLVE THE TRANSPORT EQUATION 
			e.g., c = c_old + dt*(advection+dispersion)
		*/
        gmpi.mpi_sendrecv(rt.c, gdom, par);
		/*
			UPDATE TIME STEP, dt
			--> result named dt_tmp
		*/
        ierr = MPI_Allreduce(&dt_tmp, &rt.dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        Kokkos::parallel_for(gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {rt.c(iGlob,0) = rt.c(iGlob,1);});
    }

    
};

#endif
