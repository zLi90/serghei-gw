#include "../serghei.h"
class damBreakConfig{
	public:
	geometry::point p; //point defining the dam break
	real size; //distance defining a dam break
	real hLow, hHigh; // dam break depths
};

void createDamBreakX(Domain &dom, damBreakConfig &db, int nx_override){
	// domain definition
	dom.nx_glob = 100;
  if(nx_override > 0) dom.nx_glob = nx_override;
	dom.ny_glob = 4;
	dom.xll = dom.yll = 0.;
	dom.dxConst = 0.5;
	dom.BCtype = BC_REFLECTIVE;
	// minimal parameters
	dom.simLength= 8;
	// initial conditions
	db.p(_X)= dom.nx_glob * dom.dxConst * 0.5;	
	db.hHigh = 4.;	
	db.hLow = 1.;
}

void damBreakX(const Domain &dom, State &state, const damBreakConfig &db){
	Kokkos::parallel_for("damBreakCreator",dom.nCell,KOKKOS_LAMBDA (int iGlob){
 	 	int ii = dom.getIndex(iGlob);
		state.h(ii) = db.hLow;
		if(dom.getCellCenter(iGlob)(_X) <= db.p(_X)) state.h(ii) = db.hHigh;
	});
}

void createDamBreakCircle(Domain &dom, damBreakConfig &db, int nx_override){
	// domain definition
	dom.nx_glob = 100;
	dom.ny_glob = 100;
  dom.xll = dom.yll = 0.;
  if(nx_override > 0){
    dom.nx_glob = nx_override;
    dom.ny_glob = nx_override;
  }
	dom.dxConst = 0.5;
	dom.BCtype = BC_REFLECTIVE;
	// minimal parameters
	dom.simLength= 8;
	// initial conditions
	db.p(_X)= dom.nx_glob * dom.dxConst * 0.5;	
	db.p(_Y)= dom.ny_glob * dom.dxConst * 0.5;	
	db.hHigh = 4.;	
	db.hLow = 1.;
	db.size = dom.nx_glob*0.5 / 5.; // radius 
}

void damBreakCircle(const Domain &dom, State &state, const damBreakConfig &db){
	Kokkos::parallel_for("damBreakCreator",dom.nCell,KOKKOS_LAMBDA (int iGlob){
 	 	int ii =  dom.getIndex(iGlob);
		state.h(ii) = db.hLow;
    real d = geometry::distance(dom.getCellCenter(iGlob),db.p);
		if(d <= db.size) state.h(ii) = db.hHigh;
	});
}

int main(int argc, char** argv) {
  int minArgs=3;
  int maxArgs=6;
  int mmaxArgs=8;
  int nx_override=-1;
	
	
	#if SERGHEI_SWE_MODEL==0
		std::cout << RERROR << "SERGHEI was compiled with SERGHEI_SWE_MODEL=OFF, therefore this test cannot be executed" << std::endl;
		return 0;
	#endif

	// Read command line arguments 
	 if (argc < minArgs || argc > mmaxArgs){
	  std::cerr << RERROR "The program is run as: ./nprogram caseName [damBreakX, damBreakCircle] Nthreads [[Nx,parNx,parNy,simLength,outFreq]]" << std::endl;
		  return 1;
	}
   if (argc > minArgs){
    if(argc < maxArgs){
	  std::cerr << RERROR "If the optional domain arguments Nx, parNx and parNy, are included, then all of them must be given" << std::endl;
		  return 1;  
    }else if(argc > maxArgs && argc < mmaxArgs){
      std::cerr << RERROR "Optional arguments simLength and outFreq expected but only one provided" << std::endl;
      return 1;
    }
   }

	{
	int found=0;
	SERGHEI serghei;
  
  // domain decomposition defaults
	serghei.par.nproc_x = 1;
	serghei.par.nproc_y = 1;

  real simLength = -1;
  real outFreq = -1;

  // override defaults
  if(argc > minArgs){
    nx_override = atoi(argv[3]);
    serghei.par.nproc_x = atoi(argv[4]);
    serghei.par.nproc_y = atoi(argv[5]);
    std::cout << GOK "Overriding of default problem settings with nx = " << nx_override << ", nproc_x = " << serghei.par.nproc_x <<  ", nproc_y = " << serghei.par.nproc_y << std::endl;
    if(argc > maxArgs){
        simLength = atof(argv[6]);
        outFreq = atof(argv[7]);
        std::cout << GOK "Overriding of default problem settings with simLength = " << simLength << " and outFreq = " << outFreq << std::endl;
    }
  }
	
  std::string testCase = argv[1];


  std::cout << BOLD "INTEGRATION TEST - " << testCase << RESET << std::endl;

	serghei.outFolder = "./" + testCase + "/";

	serghei.par.nthreads = atoi(argv[2]);

	serghei.init.read = 0;	// override reading of files

	damBreakConfig db;



	if(!testCase.compare("damBreakX")){
		createDamBreakX(serghei.dom, db, nx_override);
		found++;
	}
	if(!testCase.compare("damBreakCircle")){
		createDamBreakCircle(serghei.dom,db, nx_override);
		found++;
	}
    if(simLength >= 0) serghei.dom.simLength = simLength;
	
    // minimal parameters
	serghei.dom.cfl = 0.5;
	serghei.io.outFreq =  floor(serghei.dom.simLength/8);
	serghei.io.outFormat = OUT_NETCDF;
	serghei.io.allowIni = 0;
    
  if(outFreq > 0) serghei.io.outFreq = outFreq;
  if(outFreq == 0){
    serghei.io.outFreq = simLength;
    serghei.io.outFormat = -100;
  }
    

	if(!found){
		std::cerr << RERROR << "Invalid test case " << std::endl;
		return 1;
	}

	
	if(!serghei.start(argc, argv)) return 1;
	
	if(serghei.par.masterproc){
    std::cout << GOK << "Number of cells: "<< serghei.dom.nCellGlobal << std::endl;
	  std::cout << GOK << "Domain size: " << serghei.dom.nx_glob << "x" << serghei.dom.ny_glob << std::endl;
	  std::cout << GOK << "Number of MPI ranks: " << serghei.par.nranks <<std::endl;
	  std::cout << BDASH << "Creating initial state" << std::endl;
  }

	if(!testCase.compare("damBreakX")) damBreakX(serghei.dom, serghei.state, db);
	if(!testCase.compare("damBreakCircle")) damBreakCircle(serghei.dom, serghei.state, db);
  
  serghei.exch.exchangeIniMPI(serghei.state,serghei.dom,serghei.exch,serghei.par);
  serghei.init.boundaryIni(serghei.state,serghei.dom,serghei.par,serghei.ebc.extbc);

  serghei.io.outputIni(serghei.state, serghei.dom, serghei.ss.swss, serghei.par,serghei.outFolder);
	
	if(!serghei.compute()) return 1;
	if(!serghei.finalise()) return 1;
  
	} // scope guard required to ensure serghei destructor is called

	Kokkos::finalize();
  MPI_Finalize();

	return 0;
}