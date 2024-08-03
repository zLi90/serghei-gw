#include "../serghei.h"
class parabolicBowlConfig{
	public:
	real h0 = 0.10;
	real L = 4.0;
	real a = 1.0;
	real eta = 0.5;
	real omega = sqrt(2.0 * GRAV * h0) / a;
	real period;
	real r0 = 0.80;
	real A = (a*a-r0*r0)/(a*a+r0*r0);
};

void createParabolicBowl(Domain &dom, parabolicBowlConfig &pb){
	// domain definition
	dom.nx_glob = 200;
	dom.ny_glob = 200;
	dom.dxConst = pb.L/dom.nx_glob;
	dom.BCtype = BC_REFLECTIVE;
	// minimal parameters
	pb.period = 2.*M_PI*pb.omega;
	dom.simLength= 2*pb.period;
}

KOKKOS_INLINE_FUNCTION double radius(double &x, double &y, const parabolicBowlConfig &pb){
  return sqrt((x-pb.L*0.5)*(x-pb.L*0.5) + (y-pb.L*0.5)*(y-pb.L*0.5));

}
KOKKOS_INLINE_FUNCTION double paraboloid(double &x, double &y, const parabolicBowlConfig &pb){
	double r = radius(x,y,pb);
  return  -pb.h0 * (1. - r*r/pb.a/pb.a);
}

void parabolicBowlPlanar(const Domain &dom, State &state, parabolicBowlConfig &pb){
	Kokkos::parallel_for("parabolicBowlPlanarCreator",dom.nCell,KOKKOS_LAMBDA (int iGlob){
 	 	int ii = dom.getIndex(iGlob);
    double x = dom.getCellCenter(iGlob)(_X);
    double y = dom.getCellCenter(iGlob)(_Y);
    double r = radius(x,y,pb);
    state.z(ii) = paraboloid(x,y,pb);
		state.h(ii) = (pb.eta * pb.h0 / (pb.a * pb.a)) * (2.0 * (x - 0.5*pb.L) * cos(pb.omega * dom.etime) - 2.0 * (y - 0.5*pb.L) * sin(pb.omega * dom.etime) - pb.eta) - state.z(ii);
		if(state.h(ii) < 0) state.h(ii) = 0.;
    state.hu(ii) = -pb.eta * pb.omega * sin(pb.omega * dom.etime) * state.h(ii);
    state.hv(ii) = pb.eta * pb.omega * cos(pb.omega * dom.etime) * state.h(ii);
	});
}

void parabolicBowlRadial(const Domain &dom, State &state, parabolicBowlConfig &pb){
	Kokkos::parallel_for("parabolicBowlRadialCreator",dom.nCell,KOKKOS_LAMBDA (int iGlob){
 	 	int ii = dom.getIndex(iGlob);
    double x = dom.getCellCenter(iGlob)(_X);
    double y = dom.getCellCenter(iGlob)(_Y);
    double r = radius(x,y,pb);
    state.z(ii) = paraboloid(x,y,pb);

		real Aa = 1. - pb.A*pb.A;
		real cAa = 1.0 - pb.A*cos(pb.omega*dom.etime);
		state.h(ii) = pb.h0*((sqrt(Aa)/cAa)-1.-r*(Aa/(cAa*cAa)-1.)) - state.z(ii);
		if(state.h(ii) < 0) state.h(ii) = 0.;
		cAa = 1.-pb.A*cos(pb.omega*dom.etime);
    state.hu(ii) = 1./cAa * 0.5 * pb.omega * (x - 0.5*pb.L) * pb.A * sin(pb.omega * dom.etime)*state.h(ii);
    state.hv(ii) = 1./cAa * 0.5 * pb.omega * (y - 0.5*pb.L) * pb.A * sin(pb.omega * dom.etime)*state.h(ii);
	});
}

int main(int argc, char** argv) {
	
	#if SERGHEI_SWE_MODEL==0
		std::cout << RERROR << "SERGHEI was compiled with SERGHEI_SWE_MODEL=OFF, therefore this test cannot be executed" << std::endl;
		return 0;
	#endif

	// Read command line arguments 
	 if (argc != 3){
	  std::cerr << RERROR "The program is run as: ./nprogram caseName [planar,radialSymmetric] Nthreads" << std::endl;
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

	parabolicBowlConfig pb;

	// domain decomposition
	serghei.par.nproc_x = 1;
	serghei.par.nproc_y = 1;
	serghei.par.nranks = 1;
	// minimal parameters
	serghei.dom.cfl = 0.5;
	serghei.io.nScreen = 200;
	serghei.io.outFormat = OUT_NETCDF;
	serghei.io.allowIni = 0;

	if(!testCase.compare("planar")){
		serghei.state.hmin = 1E-4;
		createParabolicBowl(serghei.dom, pb);
		found++;
	}
	if(!testCase.compare("radialSymmetric")){
		pb.omega = sqrt(8.0 * GRAV * pb.h0) / pb.a;
		serghei.state.hmin = 1E-10;
		createParabolicBowl(serghei.dom, pb);
		found++;
	}
	
	serghei.io.outFreq =  pb.period*0.25;

	if(!found){
		std::cerr << RERROR << "Invalid test case " << std::endl;
		return 1;
	}
	

	std::cout << GOK << "Domain size: " << serghei.dom.nx_glob << "x" << serghei.dom.ny_glob << std::endl;
	std::cout << GOK << "Number of MPI ranks: " << serghei.par.nranks <<std::endl;
	

	std::cout << GOK << "Oscillation period : " << pb.period << " s" << std::endl;
	std::cout << GOK << "Output frequency: " << serghei.io.outFreq << " s" << std::endl;

	if(!serghei.start(argc, argv)) return 1;
	
	std::cout << GOK << "Number of cells: "<< serghei.dom.nCell << std::endl;
	std::cout << BDASH << "Creating initial state" << std::endl;

	if(!testCase.compare("planar")) parabolicBowlPlanar(serghei.dom, serghei.state, pb);
	if(!testCase.compare("radialSymmetric")) parabolicBowlRadial(serghei.dom, serghei.state, pb);


  serghei.io.outputIni(serghei.state, serghei.dom, serghei.ss.swss, serghei.par,serghei.outFolder);
	
	if(!serghei.compute()) return 1;
	if(!serghei.finalise()) return 1;
  
	} // scope guard required to ensure serghei destructor is called

	Kokkos::finalize();
  MPI_Finalize();

	return 0;
}
