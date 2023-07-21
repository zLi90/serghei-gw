#ifndef _EXCHANGE_H_
#define _EXCHANGE_H_

#include "define.h"
#include "mpi.h"
#include "Indexing.h"

class Exchange {

protected:

  int const maxPack = 3; //the number of sw variables

  MPI_Request sReq[2];
  MPI_Request rReq[2];

  MPI_Status  sStat[2];
  MPI_Status  rStat[2];

  int nPack;
  int nUnpack;

  realArr haloSendBufS;
  realArr haloSendBufN;
  realArr haloSendBufW;
  realArr haloSendBufE;
  realArr haloRecvBufS;
  realArr haloRecvBufN;
  realArr haloRecvBufW;
  realArr haloRecvBufE;

public:
  real MPItime=0;
  real exchangeTime=0;

  inline void allocate(Domain &dom) {
	// haloSendBuf* and haloRecvBuf* arrays are 1D containers (arrays) which should be able to store all state variables, for the exchanging band
    haloSendBufS = realArr("haloSendBufS",maxPack*haloc*dom.nx);
    haloSendBufN = realArr("haloSendBufN",maxPack*haloc*dom.nx);
    haloSendBufW = realArr("haloSendBufW",maxPack*dom.ny*haloc);
    haloSendBufE = realArr("haloSendBufE",maxPack*dom.ny*haloc);
    haloRecvBufS = realArr("haloRecvBufS",maxPack*haloc*dom.nx);
    haloRecvBufN = realArr("haloRecvBufN",maxPack*haloc*dom.nx);
    haloRecvBufW = realArr("haloRecvBufW",maxPack*dom.ny*haloc);
    haloRecvBufE = realArr("haloRecvBufE",maxPack*dom.ny*haloc);
  }


  inline void haloInit() {
    nPack   = 0;
    nUnpack = 0;
  }

	// packs state variable a into haloSendBuf arrays
  inline void haloPack_x(Domain &dom, realArr &a) {
    haloPack_x_ext(dom, a, haloSendBufW, haloSendBufE, nPack);
    nPack = nPack + 1;		// to keep track of how many variable are getting packed
  }
  inline void haloPack_x_ext(Domain &dom, realArr &a, realArr &haloSendBufW, realArr &haloSendBufE, int const nPack) {
	// span the x-halo columns
    Kokkos::parallel_for("haloPack_x_span", dom.ny*haloc , KOKKOS_LAMBDA (int iGlob) {
		int rx,ry;
		int nGlob = dom.ny*haloc;
	 	unpackIndices(iGlob,dom.ny,haloc,ry,rx);
	   int ii1=(haloc+ry)*(dom.nx+2*haloc)+haloc+rx;		// west
	   int ii2=(haloc+ry)*(dom.nx+2*haloc)+haloc+dom.nx-1+rx;	// east
      haloSendBufW(nPack*nGlob+iGlob) = a(ii1);
      haloSendBufE(nPack*nGlob+iGlob) = a(ii2);
    });
  }


  inline void haloPack_y(Domain &dom, realArr &a) {
    haloPack_y_ext(dom, a, haloSendBufS, haloSendBufN, nPack);
    nPack = nPack + 1;       // to keep track of how many variable are getting packed
  }
  inline void haloPack_y_ext(Domain &dom, realArr &a, realArr &haloSendBufS, realArr &haloSendBufN, int const nPack) {
	  	// span the y-halo rows
    Kokkos::parallel_for("haloPack_y_span", haloc*dom.nx , KOKKOS_LAMBDA (int iGlob) {
	 	int rx,ry;
		int nGlob = haloc*dom.nx;
		unpackIndices(iGlob,haloc,dom.nx,ry,rx);
	 	int ii1=(haloc+ry)*(dom.nx+2*haloc)+haloc+rx;			// north
	 	int ii2=(haloc+dom.ny-1+ry)*(dom.nx+2*haloc)+haloc+rx;	// south
		haloSendBufN(nPack*nGlob+iGlob) = a(ii1);
		haloSendBufS(nPack*nGlob+iGlob) = a(ii2);

    });
  }

	// unpacks haloRecvBuf arrays into state variable a
  inline void haloUnpack_x(Domain &dom, realArr &a) {
    haloUnpack_x_ext(dom, a, haloRecvBufW, haloRecvBufE, nUnpack);
    nUnpack = nUnpack + 1;
  }
  inline void haloUnpack_x_ext(Domain &dom, realArr &a, realArr &haloRecvBufW, realArr &haloRecvBufE, int const nUnpack) {
    Kokkos::parallel_for("haloUnpack_x_span", dom.ny*haloc , KOKKOS_LAMBDA (int iGlob) {
	 	int rx,ry;
	 	int nGlob = dom.ny*haloc;
		unpackIndices(iGlob,dom.ny,haloc,ry,rx);
	   int ii1=(haloc+ry)*(dom.nx+2*haloc)+haloc-rx-1;
	   int ii2=(haloc+ry)*(dom.nx+2*haloc)+haloc+dom.nx+rx;
      a(ii1) = haloRecvBufW(nUnpack*nGlob+iGlob);
      a(ii2) = haloRecvBufE(nUnpack*nGlob+iGlob);

    });
  }


  inline void haloUnpack_y(Domain &dom, realArr &a) {
    haloUnpack_y_ext(dom, a, haloRecvBufS, haloRecvBufN, nUnpack);
    nUnpack = nUnpack + 1;
  }
  inline void haloUnpack_y_ext(Domain &dom, realArr &a, realArr &haloRecvBufS, realArr &haloRecvBufN, int const nUnpack) {
    Kokkos::parallel_for( "haloUnpack_y_span",haloc*dom.nx , KOKKOS_LAMBDA (int iGlob) {
		int rx,ry;
		int nGlob = haloc*dom.nx;
		unpackIndices(iGlob,haloc,dom.nx,ry,rx);
	 	int ii1=(haloc-ry-1)*(dom.nx+2*haloc)+haloc+rx;
	 	int ii2=(haloc+dom.ny+ry)*(dom.nx+2*haloc)+haloc+rx;
      a(ii1) = haloRecvBufN(nUnpack*nGlob+iGlob);
      a(ii2) = haloRecvBufS(nUnpack*nGlob+iGlob);

    });
  }

	// MPI wrapper to exchange east/west halo regions
	inline int haloExchange_x(Domain &dom, Parallel &par) {
	#if SERGHEI_DEBUG_MPI
	std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
	#endif
    	int ierr;
		ierr=1;


    	if (par.nproc_x > 1) {
      	Kokkos::fence();	// ensure no kernels are running, we need everything available in host memory
      	//Pre-post the receives
      	ierr = MPI_Irecv( haloRecvBufW.data() , nPack*dom.ny*haloc , SERGHEI_MPI_REAL , par.neigh(1,0) , 0 , MPI_COMM_WORLD , &rReq[0] );
      	ierr = MPI_Irecv( haloRecvBufE.data() , nPack*dom.ny*haloc , SERGHEI_MPI_REAL , par.neigh(1,2) , 1 , MPI_COMM_WORLD , &rReq[1] );

     	//Send the data
      	ierr = MPI_Isend( haloSendBufW.data() , nPack*dom.ny*haloc , SERGHEI_MPI_REAL , par.neigh(1,0) , 1 , MPI_COMM_WORLD , &sReq[0] );
      	ierr = MPI_Isend( haloSendBufE.data() , nPack*dom.ny*haloc , SERGHEI_MPI_REAL , par.neigh(1,2) , 0 , MPI_COMM_WORLD , &sReq[1] );

      	//Wait for the sends and receives to finish
      	ierr = MPI_Waitall(2, sReq, sStat);
      	ierr = MPI_Waitall(2, rReq, rStat);

    	}

		// the outer boundaries are assumed to be periodic by default in the previous exchange
		// we now need to manage boundary conditions (halo) for outer domain

		if(dom.BCtype==BC_TRANSMISSIVE){
			if(dom.iW){
				haloTransmissive(nPack*dom.ny*haloc, haloSendBufW, haloRecvBufW);
			}
			if(dom.iE){
				haloTransmissive(nPack*dom.ny*haloc, haloSendBufE, haloRecvBufE);
			}
		}
		if(dom.BCtype==BC_REFLECTIVE){
			if(dom.iW){
				haloReflective(nPack*dom.ny*haloc, haloRecvBufW);
			}
			if(dom.iE){
				haloReflective(nPack*dom.ny*haloc, haloRecvBufE);
			}
		}


		return ierr;

	}


	inline int haloExchange_y(Domain &dom, Parallel &par) {
	#if SERGHEI_DEBUG_MPI
	std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
	#endif
   	int ierr;
		ierr=1;

    	if (par.nproc_y > 1) {
      	Kokkos::fence();

			//Pre-post the receives
			ierr = MPI_Irecv( haloRecvBufN.data() , nPack*haloc*dom.nx , SERGHEI_MPI_REAL , par.neigh(0,1) , 0 , MPI_COMM_WORLD , &rReq[0] );
			ierr = MPI_Irecv( haloRecvBufS.data() , nPack*haloc*dom.nx , SERGHEI_MPI_REAL , par.neigh(2,1) , 1 , MPI_COMM_WORLD , &rReq[1] );

			//Send the data
			ierr = MPI_Isend( haloSendBufN.data() , nPack*haloc*dom.nx , SERGHEI_MPI_REAL , par.neigh(0,1) , 1 , MPI_COMM_WORLD , &sReq[0] );
			ierr = MPI_Isend( haloSendBufS.data() , nPack*haloc*dom.nx , SERGHEI_MPI_REAL , par.neigh(2,1) , 0 , MPI_COMM_WORLD , &sReq[1] );

			//Wait for the sends and receives to finish
			ierr = MPI_Waitall(2, sReq, sStat);
			ierr = MPI_Waitall(2, rReq, rStat);

    	}

		// the outer boundaries are assumed to be periodic by default in the previous exchange
		// we now need to manage boundary conditions (halo) for outer domain

		if(dom.BCtype==BC_TRANSMISSIVE){
			if(dom.iN){
				haloTransmissive(nPack*haloc*dom.nx, haloSendBufN, haloRecvBufN);
			}
			if(dom.iS){
				haloTransmissive(nPack*haloc*dom.nx, haloSendBufS, haloRecvBufS);
			}
		}
		if(dom.BCtype==BC_REFLECTIVE){
			if(dom.iN){
				haloReflective(nPack*haloc*dom.nx, haloRecvBufN);
			}
			if(dom.iS){
				haloReflective(nPack*haloc*dom.nx, haloRecvBufS);
			}
		}

		return ierr;
	}


  /*inline void haloPeriodic(const int ncells, realArr const &haloSendBuf1, realArr const &haloSendBuf2, realArr &haloRecvBuf1, realArr &haloRecvBuf2) {
    Kokkos::parallel_for("haloPeriodic", ncells , KOKKOS_LAMBDA (int iGlob) {
		haloRecvBuf1(iGlob) = haloSendBuf2(iGlob);
      haloRecvBuf2(iGlob) = haloSendBuf1(iGlob);
    });
  }*/

  inline void haloReflective(const int ncells, realArr &haloRecv) {
    Kokkos::parallel_for("haloReflective", ncells , KOKKOS_LAMBDA (int iGlob) {
      haloRecv(iGlob) = 0.0;
    });
  }
  inline void haloTransmissive(const int ncells, realArr const &haloSend, realArr &haloRecv) {
    Kokkos::parallel_for("haloTransmissive", ncells , KOKKOS_LAMBDA (int iGlob) {
		haloRecv(iGlob) = haloSend(iGlob);
    });
  }

	// high level wrapper to do the initial exchange of surface parameters
	inline void exchangeIniMPI(State &state, Domain &dom, Exchange &exch, Parallel &par){
	#if SERGHEI_DEBUG_MPI
	std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
	#endif
	    // Allocate the MPI exchange buffers
    	exch.allocate(dom);

		 // Exchange surface elevation in x-direction
		 exch.haloInit      ();
		 exch.haloPack_x   (dom, state.z);
		 exch.haloPack_x   (dom, state.roughness);
		 exch.haloExchange_x(dom, par);
		 exch.haloUnpack_x (dom, state.z);
		 exch.haloUnpack_x (dom, state.roughness);

		 // Exchange surface elevation in y-direction
		 exch.haloInit      ();
		 exch.haloPack_y   (dom, state.z);
		 exch.haloPack_y   (dom, state.roughness);
		 exch.haloExchange_y(dom, par);
		 exch.haloUnpack_y (dom, state.z);
		 exch.haloUnpack_y (dom, state.roughness);

		exch.exchangeMPIh(state,dom,exch,par);
		exch.exchangeMPIhuhv(state,dom,exch,par);


	}

	// high level wrapper for water depth
	inline void exchangeMPIh(State &state, Domain &dom, Exchange &exch, Parallel &par){
	#if SERGHEI_DEBUG_MPI
	std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
	#endif
		Kokkos::Timer timer;	// only to keep track of time

		 // Exchange depth in x-direction
		 exch.haloInit      ();
		 exch.haloPack_x   (dom, state.h);		// re-orders the entries in state.h array which are on the east/west halo regions into data packs which will be sent east and west
		 exch.haloExchange_x(dom, par);			// MPI send/receives the data packs, and handles boundaries
		 exch.haloUnpack_x (dom, state.h);		// re-order the updated data packs back into the halo regions in state.h array

		 // Exchange depth in y-direction
		 exch.haloInit      ();
		 exch.haloPack_y   (dom, state.h);		// re-orders the entries in state.h array which are on the north/south halo regions into data packs which will be sent north and south
		 exch.haloExchange_y(dom, par);			// MPI send/receives the data packs, and handles boundaries
		 exch.haloUnpack_y (dom, state.h);		// re-order the updated data packs back into the state.h array
    dom.timers.exchange += timer.seconds();
	}

	// high level wrapper for momentum
	inline void exchangeMPIhuhv(State &state, Domain &dom, Exchange &exch, Parallel &par){
	#if SERGHEI_DEBUG_MPI
	std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
	#endif

		Kokkos::Timer timer;

		 // Exchange momentum in x-direction
		 exch.haloInit      ();
		 exch.haloPack_x   (dom, state.hu);
		 exch.haloPack_x   (dom, state.hv);
		 exch.haloExchange_x(dom, par);
		 exch.haloUnpack_x (dom, state.hu);
		 exch.haloUnpack_x (dom, state.hv);

		 // Exchange momentum in y-direction
		 exch.haloInit      ();
		 exch.haloPack_y   (dom, state.hu);
		 exch.haloPack_y   (dom, state.hv);
		 exch.haloExchange_y(dom, par);
		 exch.haloUnpack_y (dom, state.hu);
		 exch.haloUnpack_y (dom, state.hv);

		 dom.timers.exchange += timer.seconds();
	}





};

#endif
