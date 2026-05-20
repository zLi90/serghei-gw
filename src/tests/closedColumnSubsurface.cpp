#include "../serghei.h"
#if SERGHEI_SUBSURFACE_MODEL

void createTest(GwState &gw, Domain &dom, GwDomain &gdom)	{
	// domain definition
	dom.nx_glob = 1;
	dom.ny_glob = 1;
	dom.dxConst = 1.0;
	dom.BCtype = BC_REFLECTIVE;
	// parameters
	dom.simLength = 10;
	// subsurface
	gdom.thickH = 1.0;
	gdom.dz_multiplier = 1.0;
	gdom.dz_base = 0.05;
	gdom.nz_glob = 20;
	gdom.dt_init = 0.01;
	gdom.dt_max = 10.0;
	gdom.nSoilID = 1;
	gdom.gw_scheme = 2;
	gdom.aev = 0.0;
	gdom.async = 0;
	// soil parameters
	gdom.wcs = 0.3;	gdom.wcr = 0.0;
	gdom.alpha = 1.5; gdom.n = 2.0;
	gdom.Ks = 1e-3;

	// soil parameters
	gw.nVGparam = 7;

	// initial condition: prescribed water content
	gw.initialMode = 3;
	gdom.wc_ic = 0.2;
}

#endif

int main(int argc, char** argv) {
#if SERGHEI_SUBSURFACE_MODEL==0
		std::cout << RERROR << "SERGHEI was compiled with SERGHEI_SUBSURFACE_MODEL=OFF, therefore this test cannot be executed" << std::endl;
		return 0;
#else
	// Read command line arguments
	 if (argc != 3){
	  std::cerr << RERROR "The program is run as: ./nprogram caseName [column] Nthreads" << std::endl;
		  return 1;
	}

	{
	int found=0;
	SERGHEI serghei;

	std::string testCase = argv[1];
	std::cout << BOLD "INTEGRATION TEST - " << testCase << RESET << std::endl;

	serghei.outFolder = "./" + testCase + "/";

	serghei.par.nthreads = atoi(argv[2]);

	serghei.init.read = 0;	// override reading of files
	serghei.ginit.read = 0;

	// domain decomposition
	serghei.par.nproc_x = 1;
	serghei.par.nproc_y = 1;
	serghei.par.nranks = 1;
	// minimal parameters
	serghei.dom.cfl = 0.1;
	serghei.io.outFreq =  1;
	serghei.io.obsFreq =  1;
	serghei.io.outFormat = OUT_NETCDF;
	serghei.io.nScreen=10;

	if (!testCase.compare("column"))	{
		createTest(serghei.gw, serghei.dom, serghei.gdom);
		found++;
	}
	if(!found){
		std::cerr << RERROR << "Invalid test case " << std::endl;
		return 1;
	}

	std::cout << GOK << "Domain size: " << serghei.dom.nx_glob << "x" << serghei.dom.ny_glob << std::endl;
	std::cout << GOK << "Number of MPI ranks: " << serghei.par.nranks <<std::endl;

	if(!serghei.start(argc, argv)) return 1;

	std::cout << GOK << "Number of cells: "<< serghei.dom.nCell << std::endl;
	std::cout << BDASH << "Creating initial state" << std::endl;

	if(!serghei.compute()) return 1;
	if(!serghei.finalise()) return 1;

	} // scope guard required to ensure serghei destructor is called

	Kokkos::finalize();
	MPI_Finalize();

	return 0;
#endif
}
