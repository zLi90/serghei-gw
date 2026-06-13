#ifndef _SOURCESINK_H_
#define _SOURCESINK_H_

#include "define.h"
#include "SArray.h"
#include "State.h"

#include "GwDomain.h"
#include "GwState.h"


#define INF_NONE 0
#define INF_CONSTANT 1
#define INF_HORTON 2
#define INF_GREENAMPT 3
#define INF_SINK 4

namespace {
/** log for device/host; double intermediate for float stability */
KOKKOS_INLINE_FUNCTION real ga_log_pos(real x)
{
    return (real)log((double)x);
}

/**
 * One-step cumulative infiltration depth ΔF over ts (SWMM grnampt_getF2), then min(ΔF, ia_cap*ts).
 * - c1==0: ΔF = Ks*ts.
 * - ts<10 and F1 > 0.01*c1: explicit F2 = F1 + Ks(1+c1/F1)ts, F2 = max(F2, F1+Ks*ts).
 * - else: Newton on integrated GA F2 - c1*ln(F2+c1) = F1 - c1*ln(F1+c1) + Ks*ts.
 * F1>0 uses physical wetting depth as denominator; only F1<=0 uses Ffloor (avoids Fs<<infDry blow-up).
 */
KOKKOS_INLINE_FUNCTION real greenAmptStep_dF(real ks, real c1, real F1, real Ffloor, real ts, real ia_cap)
{
    if (ts <= (real)TOL_MACHINE_ACCURACY) { return (real)0.0; }
    real cap = ia_cap * ts;
    real f2min = F1 + ks * ts;
    if (c1 <= (real)0.0) {
        real dF = ks * ts;
        if (dF > cap) { dF = cap; }
        return dF;
    }
    real f2;
    if (ts < (real)10.0 && F1 > (real)0.01 * c1) {
        real Fd = (F1 > (real)0.0) ? F1 : Ffloor;
        f2 = F1 + ks * ((real)1.0 + c1 / Fd) * ts;
        if (f2 < f2min) { f2 = f2min; }
    } else {
        real c2 = c1 * ga_log_pos(F1 + c1) - ks * ts;
        if (F1 > (real)0.0) {
            f2 = F1;
        } else {
            f2 = max(f2min, (real)0.01 * c1);
        }
        for (int i = 0; i < 20; ++i) {
            real fpc = f2 + c1;
            if (fpc <= (real)0.0) {
                f2 = f2min;
                break;
            }
            real denom = (real)1.0 - c1 / fpc;
            if (myfabs(denom) < (real)1e-14) { break; }
            real df2 = (f2 - F1 - c1 * ga_log_pos(fpc) + c2) / denom;
            if (myfabs(df2) < (real)0.00001) { break; }
            f2 -= df2;
        }
        if (f2 < f2min) { f2 = f2min; }
    }
    real dF = f2 - F1;
    if (dF > cap) { dF = cap; }
    return dF;
}


KOKKOS_INLINE_FUNCTION real greenAmptCellRate(real ks, real psi, real dth, real F, real pond, real r, real dt, real Ffloor)
{
    real ia = r + pond / dt;
    if (ia <= (real)TOL_MACHINE_ACCURACY) { return (real)0.0; }
    if (ia <= ks) { return ia; }
    real c1 = (psi + pond) * dth;
    real dks = ia - ks;
    real Fs = ks * c1 / dks;
    if (Fs < (real)0.0) { Fs = (real)0.0; }
    if (F + ia * dt < Fs) { return ia; }
    if (F >= Fs) { return greenAmptStep_dF(ks, c1, F, Ffloor, dt, ia) / dt; }
    real t1 = (Fs - F) / ia;
    if (t1 >= dt - (real)TOL_MACHINE_ACCURACY) { return ia; }
    real dF1 = ia * t1;
    real ts = dt - t1;
    real dF2 = greenAmptStep_dF(ks, c1, Fs, Ffloor, ts, ia);
    return (dF1 + dF2) / dt;
}
}

//class State;	// forward declaration

/*
  TimeSeries provides a construct/class to store time series.
 */
class TimeSeries{

public:

  int np;     // number of points in time
  int nc;       // number of grid cells with different time series values
  int nx = 1; // number of partitions in x direction
  int ny = 1; // number of partitions in y direction

  realArr time;
  realArr value;
  realArr2 values;
  int timeIndex = 0;

/*
  // WARNING valid only for piece-wise constant time data
  inline real interpolate (real const &t, int spaceIndex){
    if(t >= time (np - 1)){
      timeIndex = np - 1;
    }
    else{
	     if (t >= time (timeIndex + 1)) timeIndex++;
    }
    return (value (np * spaceIndex + timeIndex));
  }
*/
  void initialise(int n){
    np = n;
    time = realArr("time",np);
    value = realArr("value",np);
  };

};

KOKKOS_INLINE_FUNCTION void findTimeBlock (TimeSeries &ts, real const &t){
  if(t >= ts.time (ts.np - 1)){
    ts.timeIndex = ts.np - 1;
  }
  else{
	   if (t >= ts.time (ts.timeIndex + 1)) ts.timeIndex++;
  }
};

KOKKOS_INLINE_FUNCTION real interpolatePiecewise (TimeSeries const &ts, real const &t, int const spaceIndex ){
  return (ts.value (ts.np * spaceIndex + ts.timeIndex));
};

KOKKOS_INLINE_FUNCTION real interpolateLinear(TimeSeries &ts, real const &t){
  #if SERGHEI_DEBUG_WORKFLOW > 1
	std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET << std::endl;
  #endif
  int ii,jj;
  findTimeBlock(ts,t);
  ii = ts.timeIndex;
  jj=ii+1;
  if(ii == ts.np - 1) jj = ii;
  real v = ts.value(ii) + (ts.value(jj) - ts.value(ii))/(ts.time(jj)-ts.time(ii))*(t-ts.time(ii));
  return(v);
};

KOKKOS_INLINE_FUNCTION real interpolateValues(TimeSeries &ts, real const &t, int icol){
  int ii,jj;
  findTimeBlock(ts,t);
  ii = ts.timeIndex;
  jj=ii+1;
  if(ii == ts.np - 1) jj = ii;
  real v = ts.values(ii,icol) + (ts.values(jj,icol) - ts.values(ii,icol))/(ts.time(jj)-ts.time(ii))*(t-ts.time(ii));
  return(v);
};



class ConstantInfiltration{

private:

  real _constCap;

public:

  ConstantInfiltration(real constCap) :  _constCap(constCap) {} 	 // constructor
  real operator()()
  {
    return _constCap;
  }

};


class InfiltrationModel{

private:
/*
  KOKKOS_INLINE_FUNCTION real horton(const int ii,const real dt) const
  {
    real t = infTime(ii) + dt;
    infTime(ii) = t;
    real infCap = fc + (f0-fc)*exp(-k * t);
    return(infCap);
  }
*/
  realArr infTime;

public:
    int model = -999;
    int nLabels = 0;
    realArr constCap ;
    // Horton
    realArr k;
    realArr fc;
    realArr f0;
    // Green-Ampt (per infiltration class id, same indexing as k/fc/f0)
    realArr ks;
    realArr psi;
    realArr dtheta;
    // Sink
    real sink_area = -999;
    real sink_capacity = -999;
    realArr storage;

    real infDry = 1E-8;  // [L] threshold to consider dry for infiltration purposes

    // This is not a state variable
    // which is why it is here and not in class State.
    // It is necessary for output
    // and because it is a variable in the GreenAmpt model
    realArr infVol; // accumulated infiltration volume
    realArr rate;		// infiltration rate
    /** Green–Ampt: cumulative infiltrated depth F [L] per cell (after supply limiting). Public: updated in TimeIntegrator. */
    realArr gaF;
    intArr infLabel; // labels for heterogeneous infiltration

    // this is a function pointer which allows to redirect
    // to the specific infiltration capacity function.
    // the goal is to avoid evaluating which model to use every time step
    // real (InfiltrationModel::*capacity)(const int ii, const real dt) const;


    // Define infiltration capacity models
    /*
    real constant(const int ii, const real t) const{
    	return(constCap);
    }
    */


    void allocate(const Domain &dom){
        if(model){
          rate = realArr("rate",dom.nCellMem);
          infVol = realArr("infVol",dom.nCellMem);
          if(model == INF_HORTON) infTime = realArr("infTime",dom.nCellMem);
          if(model == INF_GREENAMPT) {
            gaF = realArr("gaF", dom.nCellMem);
            Kokkos::deep_copy(gaF, 0.0);
          }
        }
    }



  int assignModel(Parallel &par){
    int error = 0;

	#if SERGHEI_DEBUG_INFILTRATION
	std::cout << GGD << GRAY << __PRETTY_FUNCTION__ << RESET <<  "model: " << model << std::endl;
  	#endif
    switch(model){
      case INF_NONE:
		 if(par.masterproc){
                	std::cerr << BDASH << "No infiltration capacity" << std::endl;
		 }
         break;
      case INF_CONSTANT:
      	// capacity = &InfiltrationModel::constant;
        std::cerr << BDASH << "Constant infiltration capacity" << std::endl;
        for(int id=1; id<nLabels; id++){
        	if(constCap(id) < 0){
		  		if(par.masterproc){
             		std::cerr << RERROR << "Infiltration rate not found for constant infiltration model" << std::endl;
		  		}
                error++;
            }
        }
        break;
      case INF_HORTON:
      	// capacity = &InfiltrationModel::horton;
		if(par.masterproc){
        	std::cerr << BDASH << "Horton infiltration capacity" << std::endl;
		}
        for(int ii=0; ii<nLabels; ii++){
        	//std::cerr << ii << "\t" << k(ii) << "\t" << fc(ii) << "\t" << f0(ii) << std::endl;
            if(k(ii) < 0){
				if(par.masterproc){
                	std::cerr << RERROR << "Shape factor not found for Horton infiltration model" << std::endl;
				}
                error++;
            }
            if(f0(ii) < 0){
				if(par.masterproc){
               		std::cerr << RERROR << "Initial infiltration capacity not found for Horton infiltration model" << std::endl;
				}
                error++;
            }
            if(fc(ii) < 0){
				if(par.masterproc){
        	    	std::cerr << RERROR << "Asymptotic infiltration capacity not found for Horton infiltration model" << std::endl;
				}
                error++;
            }
        }
        break;
      case INF_GREENAMPT:
		 if(par.masterproc){
         	std::cerr << BDASH << "Green-Ampt infiltration capacity" << std::endl;
		 }
        for(int ii=0; ii<nLabels; ii++){
            if(ks(ii) < 0){
		 	if(par.masterproc){
            	std::cerr << RERROR << "Saturated hydraulic conductivity not found for Green-Ampt infiltration model (class " << ii << ")" << std::endl;
			}
            error++;
            }
            if(psi(ii) < 0){
		 	if(par.masterproc){
          		std::cerr << RERROR << "Average suction head not found for Green-Ampt infiltration model (class " << ii << ")" << std::endl;
			}
            error++;
            }
            if(dtheta(ii) < 0){
			if(par.masterproc){
               	std::cerr << RERROR << "Water content difference not found for Green-Ampt infiltration model (class " << ii << ")" << std::endl;
			}
         	error++;
            }
        }
         break;
         case INF_SINK:
             if(sink_area < 0 || sink_capacity < 0){
                if(par.masterproc){
                    std::cerr << RERROR << "Sink area and capacity must be > 0" << std::endl;
                }
                error++;
             }
             break;
      default:
		if(par.masterproc){
        	std::cerr << RERROR << "Error processing data in infiltration.input using infiltration model " << model << "." << std::endl;
		}
        error++;
        break;
      }
	  if(error > 0) return 0;
      return 1;
    }



    /** Constant / Horton only. Green-Ampt is handled in SourceSinkData::computeInfiltration (needs ss.rainRate). */
    inline void ComputeInfiltrationCapacity(const State &state, const Domain &dom){
      #if SERGHEI_DEBUG_WORKFLOW
        std::cerr << GGD << __PRETTY_FUNCTION__ << std::endl;
      #endif
        if(model){
                if (model == INF_GREENAMPT) {
                    return;
                }
                realArr &inf_p = rate;
                intArr infLabel = this->infLabel;
                realArr constCap = this->constCap;

                switch(model){
                    case INF_CONSTANT:
                        Kokkos::parallel_for("inf_constant", dom.nCell, KOKKOS_LAMBDA (int iGlob){
                            int ii = dom.getIndex(iGlob);
                            int id = infLabel(ii);
                            inf_p(ii) = constCap(id);
                        });
                        break;
                    case INF_HORTON:
                    {
                        realArr fc = this->fc;
                        realArr f0 = this->f0;
                        realArr k = this->k;
                        realArr &infTime_p = infTime;
                        Kokkos::parallel_for("inf_horton", dom.nCell, KOKKOS_LAMBDA (int iGlob)
                        {
                            int ii = dom.getIndex(iGlob);
                            int id = infLabel(ii);
                            real t = infTime_p(ii) + dom.dt;
                            infTime_p(ii) = t;
                            inf_p(ii) = fc(id) + (f0(id)-fc(id))*exp(-k(id) * t);
                        });
                    }
                        break;
                }
            }
    }
    
    // simple drainge model
    inline void ComputeDrain(const State &state, const Domain &dom) {
        if (model == INF_SINK) {
            realArr &inf_p = rate;
            intArr infLabel = this->infLabel;
            realArr constCap = this->constCap;
            Kokkos::parallel_for("inf_sink", dom.nCell, KOKKOS_LAMBDA (int iGlob){
                int ii = dom.getIndex(iGlob);
                int id = infLabel(ii);
                real vc, dh, flux = 0;
                if (id > 0) {
                    // copute drainge if water flows towards the drainage outlet
                    //vc = mysqrt(GRAV*state.h(ii));
                    //flux += vc * state.h(ii) * dom.dx();

                    if (state.hu(ii) < 0)   {
                        // assume flow is critical
                        vc = mysqrt(GRAV*state.h(ii+1));
                        flux += vc * state.h(ii+1) * dom.dx();
                    }
                    if (state.hu(ii-1) > 0)   {
                        vc = mysqrt(GRAV*state.h(ii-1));
                        flux += vc * state.h(ii-1) * dom.dx();
                    }
                    if (state.hv(ii) < 0)   {
                        vc = mysqrt(GRAV*state.h(ii+dom.nx+2*hc));
                        flux += vc * state.h(ii+dom.nx+2*hc) * dom.dx();
                    }
                    if (state.hv(ii-dom.nx-2*hc) > 0)   {
                        vc = mysqrt(GRAV*state.h(ii-dom.nx+2*hc));
                        flux += vc * state.h(ii-dom.nx+2*hc) * dom.dx();
                    }
                    // get the equivalent drainage depth
                    dh = flux * dom.dt * sink_area / (dom.dx() * dom.dx());
                    if (storage(ii) < sink_capacity)   {
                        inf_p(ii) = dh;
                        storage(ii) += dh;
                    }
                    else {
                        inf_p(ii) = 0.0;
                    }
                }
                else {
                    inf_p(ii) = 0.0;
                }
            });
        }
    }
};


class SourceSinkData{

public:

    TimeSeries rain, evap, wind;
    InfiltrationModel inf;
    realArr rainRate, evapRate, windspd, winddir;

    void allocateSW (Domain const &dom){
        if (dom.isRain) {rainRate  = realArr ("rainRate", dom.nCellMem);}
        if (dom.isEvap) {evapRate  = realArr ("evapRate", dom.nCellMem);}
        if (dom.isWind) {
            windspd = realArr("windspd", dom.nCellMem);
            winddir = realArr("winddir", dom.nCellMem);
        }
        if (inf.model)  {inf.allocate(dom);}
    }

  /**
   * Surface infiltration capacity into inf.rate: delegates to inf for constant/Horton;
   * Green-Ampt is implemented here so ia can use this->rainRate when dom.isRain (same as TimeIntegrator).
   */
  inline void computeInfiltration(const State &state, const Domain &dom) {
        if (!inf.model) { return; }
        if (inf.model == INF_GREENAMPT) {
            const bool read_rain = dom.isRain;
            realArr rain_p;
            if (read_rain) { rain_p = rainRate; }
            realArr &inf_p = inf.rate;
            intArr infLabel = inf.infLabel;
            realArr ks_p = inf.ks;
            realArr psi_p = inf.psi;
            realArr dth_p = inf.dtheta;
            real Ffloor = inf.infDry;
            realArr &gaF_p = inf.gaF;
            real dt = dom.dt;
            Kokkos::parallel_for("inf_greenampt", dom.nCell, KOKKOS_LAMBDA (int iGlob)
            {
                int ii = dom.getIndex(iGlob);
                int id = infLabel(ii);
                if (id <= 0) {
                    inf_p(ii) = 0.0;
                    return;
                }
                real r = (real)0.0;
                if (read_rain) { r = rain_p(ii); }
                inf_p(ii) = greenAmptCellRate(
                    ks_p(id), psi_p(id), dth_p(id), gaF_p(ii), state.h(ii), r, dt, Ffloor);
            });
            return;
        }
        inf.ComputeInfiltrationCapacity(state, dom);
  }


  inline void ComputeRain (const Domain &dom){
    if(dom.isRain){

	// ----------------------------------------------------------------------
	// get global values to map to the correct rain subdomain
	// ----------------------------------------------------------------------
	//int nx = dom.nx_glob; // computational cell number in x direction
  //int ny = dom.ny_glob; // computational cell number in y direction

	int rainx = rain.nx;  // rain subdomain number in x direction
	int rainy = rain.ny;  // rain subdomain number in y direction

	int intervalx = dom.nx / rainx; // approximate number of cells in
				    // a subdomain in x direction
	int intervaly = dom.ny / rainy; // approximate number of cells in
				    // a subdomain in y direction
	// ----------------------------------------------------------------------


	realArr &rr_p = rainRate;

  findTimeBlock(rain,dom.etime);
  TimeSeries rrain = rain;

   Kokkos::parallel_for("rain_interpolation",dom.nCell, KOKKOS_LAMBDA (int iGlob){
	    int ix;
	    int iy;

      // dom.unpackIndices (iGlob, iy, ix);
      // int ii = dom.getHaloExtension(ix,iy);
      unpackIndicesUniformGrid(iGlob, dom.ny, dom.nx, iy, ix);
      int ii = (hc+iy)*(dom.nx+2*hc) + hc+ix;

	    int _x = ix / intervalx;
	    int _y = iy / intervaly;

	    int rain_glob = _x + _y * rainx;

	    real rainValue = interpolatePiecewise(rrain, dom.etime, rain_glob);

	    rr_p(ii) = rainValue;

	    #if SERGHEI_DEBUG_RAINFALL
        std::cerr << GGD "_x : " << _x << " _j: " << _y << " ix: " << ix << ", iy: " << iy << " ~> rainfall " << rr_p (iGlob) << std::endl;
        std::cerr << GGD "rain_glob " << rain_glob << std::endl;
      #endif
    });

	/** basically the same as above but parallel for-ized. this
	    works with MPI. we may think about a switch that uses this
	    portion of code when compiled for CPU.

	Kokkos::parallel_for (dom.nCell, KOKKOS_LAMBDA (int iGlob)
			      {

				int ix; // global x coordinate
				int iy; // global y coordinate
				unpackIndices (iGlob, ny, nx, iy, ix);

				int _x = ix / intervalx;
				int _y = iy / intervaly;

				int rain_glob = _x + _y * rainx;

				real rainValue = rain.interpolate (dom.etime, rain_glob);

				int ii = getIndex (iGlob, dom);
				rr_p (ii) = rainValue;

			      });
	**/

      }
	  #if SERGHEI_DEBUG_RAINFALL
   		std::cerr << GGD "-----------" << std::endl;;
	  #endif
  }

  inline void ComputeEvap (const Domain &dom){
      if(dom.isEvap){
          realArr &rr_e = evapRate;
          findTimeBlock(evap,dom.etime);
          TimeSeries revap = evap;
          Kokkos::parallel_for("evap_interpolation",dom.nCell, KOKKOS_LAMBDA (int iGlob){
              int ix;
              int iy;
              dom.unpackIndices (iGlob, iy, ix);
              int ii = dom.getHaloExtension(ix,iy);
              int evap_glob = 0;
              real evapValue = interpolatePiecewise(revap, dom.etime, evap_glob);
              rr_e(ii) = evapValue;
          });
      }
  }
  
  
    inline void ComputeWind(const Domain &dom)    {
        if (dom.isWind)   {
            realArr &rr_w = windspd;
            realArr &rr_d = winddir;
            findTimeBlock(wind, dom.etime);
            TimeSeries rwind = wind;
            Kokkos::parallel_for("wind_interpolation", dom.nCell, KOKKOS_LAMBDA (int idom){
                int ix, iy;
                dom.unpackIndices (idom, iy, ix);
                int ii = dom.getHaloExtension(ix,iy);
                real spdValue = interpolatePiecewise(rwind, dom.etime, 0);
                rr_w(ii) = spdValue;
                real dirValue = interpolatePiecewise(rwind, dom.etime, 1);
                rr_d(ii) = dirValue;
            });
        }
    }
    

  inline void ComputeSWSourceSink(const State &state, const Domain &dom){
    Kokkos::Timer timer;
    #if SERGHEI_DEBUG_WORKFLOW
    std::cerr << GGD << __PRETTY_FUNCTION__ << std::endl;
    #endif
    ComputeRain(dom);
    ComputeEvap(dom);
    ComputeWind(dom);
    if (inf.model == INF_SINK)  {
        inf.ComputeDrain(state, dom);
    }
    else {
        computeInfiltration(state, dom);
    }
    //no rate correction is necessary here beacuse the rate correction is done in ComputeNewState, according to the new water depth
   // timerRainInf += timer.seconds();
   dom.timers.raininf += timer.seconds();
 }

};



// Subsurface source/sinks
#if SERGHEI_SUBSURFACE_MODEL
class GwSS{

public:
	// subsurface ss directions
	#define XPLUS 1
	#define XMINUS 2
	#define YPLUS 3
	#define YMINUS 4
	#define ZPLUS 5
	#define ZMINUS 6
	// The type of source/sink
    //  0 : Evapotranspiration (from PM equation)
    //  1 : Flux (constant or time series)
    //  2 : Head for drainage (constant or time series)

    int sstype;
    int direction; // direction of source/sink, only needed for drainage ss
    int ndepth;     // number of cells in the vertical direction of the polygon
	int ncellsIT = 0;	// number of internal source/sink cells
	intArr icells; //array of indexes of boundary cells
    realArr ssvals, ssdata;
    real Qinflow, Qoutflow, Cpipe;
    TimeSeries ts;
    TimeSeries evap, tran;
	
	// water stress and root distribution function
    real lai;
	real h1, h2, h3, h4;
	real px, py, pz, xs, ys, zs, xm, ym, zm;
	realArr coef_wat, coef_root;

	MPI_Comm comm;	// communicator for ranks associated to the BC

    void allocateGW (GwDomain const &gdom)  {
        ssdata = realArr ("ssdata", gdom.nCellMem);
        for (int idx = 0; idx < gdom.nCellMem; idx++)   {ssdata(idx) = 0.0;}
		if (sstype == 0)	{
			coef_wat = realArr("wat", ncellsIT);
			coef_root = realArr("root", ncellsIT);
		}
    }

    // find internal cells for applying source/sink conditions
	inline int find_icells(GwState &gw, std::string &id, GwDomain &gdom, Parallel &par, int nPoly, realArr &xPoly, realArr &yPoly, realArr &zPoly){
		int foundInSubdom; // to keep track of which subdomains are associated to this boundary
		std::vector<int> tmpicells; //array of indexes of internal cells
		std::vector<int> tmpgcells; //array of indexes of ghost cells
		std::vector<int> subdomains;	// keeps track of which subdomains are associated to the BC
		// Loop over the entire domain to find internal source/sink cells
        int kmax = 0, kmin = gdom.nz, idx = 0;
        real zm_min = -1e5, zm_max = 0.0;
		for (int kk = 0; kk < gdom.nz; kk++) {
			for (int jj = 0; jj < gdom.ny; jj++) {
				for (int ii = 0; ii < gdom.nx; ii++) {
					int iGlob = (hc+kk)*gdom.nxhc*gdom.nyhc + (hc+jj)*gdom.nxhc + ii + hc;
					foundInSubdom = -1;
		            real xCoord = gdom.xll + ( par.i_beg + ii + 0.5) * gdom.dx;
		            real yCoord = gdom.yll + gdom.ny_glob*gdom.dx - ( par.j_beg + jj + 0.5) * gdom.dx;
					real zCoord = gdom.z(iGlob);
		            if (geometry::isInsidePoly3D(nPoly, xPoly, yPoly, zPoly, xCoord, yCoord, zCoord)){
						tmpicells.push_back(iGlob);
                        if (kk > kmax)  {kmax = kk;}
                        if (kk < kmin)  {kmin = kk;}
                        if (zCoord < zm_max)    {zm_max = zCoord;}
                        if (zCoord > zm_min)    {zm_min = zCoord;}
		            }
				}
			}
		}
        // For now, only consider 1D(z) root distribution, 20240818
        xm = 0.0;   ym = 0.0;   zm = 0.0;
        if (kmax > kmin)    {ndepth = kmax - kmin;  zm = zm_min - zm_max;}
        else {ndepth = 1;}

		ncellsIT=int(tmpicells.size());
		if(ncellsIT > 0) foundInSubdom = par.myrank; // if at least one cell in this subdomain (rank) is in the BC, tag as found

		int ncells_all;
        MPI_Allreduce(&ncellsIT, &ncells_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		int *subdoms;
		subdoms = (int*) malloc(par.nranks * sizeof(int));
		MPI_Allgather(&foundInSubdom,1,MPI_INT,subdoms,1,MPI_INT,MPI_COMM_WORLD);
		for (int i=0; i<par.nranks; i++){
			if (subdoms[i] >= 0)	{subdomains.push_back(subdoms[i]);}
		}
		MPI_Group group, subgroup;
		MPI_Comm_group(MPI_COMM_WORLD,&group);
		MPI_Group_incl(group,subdomains.size(),subdomains.data(),&subgroup);
		MPI_Comm_create(MPI_COMM_WORLD,subgroup,&comm);
		//we need the total internal cells detected by all subdomain to launch an error otherwise
		if(ncells_all>0){
			icells=intArr("icells", ncellsIT);
			#ifdef __NVCC__
				cudaMemcpyAsync( icells.data() , tmpicells.data() , ncellsIT*sizeof(int) , cudaMemcpyHostToDevice );
				cudaDeviceSynchronize();
			#else
				std::memcpy(icells.data(), tmpicells.data(), ncellsIT*sizeof(int));
			#endif
		}
		else{
			if(par.masterproc){
				std::cerr << RERROR << "No internal cells found for subsurface boundary with id '" << id << "'" << std::endl;
			}
			return 0;
		}
		return 1;
	}
	
	// get coefficients for root water uptake declining and root distribution 
	void rootCoef(const GwState &gw, const GwDomain &gdom)	{
		Kokkos::parallel_for("root", ncellsIT, KOKKOS_CLASS_LAMBDA (int idx){
			real c_wat = 1.0, c_root = 1.0, expo, x, y, z;
			int iGlob = icells[idx];
			real h = gw.h(iGlob,1);
			// get coef_wat
			if (h <= h1 && h > h2)	{c_wat = (h - h1) / (h2 - h1);}
			else if (h <= h3 && h > h4)	{c_wat = (h - h4) / (h3 - h4);}
			else if (h > h1 || h <= h4)	{c_wat = 0.0;}
			// get coordinates x, y, z
			x = gdom.x(iGlob);	y = gdom.y(iGlob);	z = gdom.depth(iGlob);
			// get coef_root
			//expo = px/(xm*myfabs(xs-x)) + py/(ym*myfabs(ys-y)) + pz/(zm*myfabs(zs-z));
			//c_root = (1.0-x/xm)*(1.0-y/ym)*(1.0-z/zm)*exp(-expo);
			expo = pz/(zm*myfabs(zs-z));
			c_root = (1.0-z/zm)*exp(-expo);
            if (c_root < 0.0)   {c_root = 0.0;}
			coef_wat(idx) = c_wat;
			coef_root(idx) = c_root;
		});
	}

    inline void applyMatSS(GwState &gw, GwDomain &gdom) {
        if (ncellsIT > 0)   {
        	// ET (Penman-Monteith)
        	if (sstype == 0)	{
        		real qt = -interpolateLinear(tran, gdom.etime);
                real qe = -interpolateLinear(evap, gdom.etime);
				// get root function coefficients 
				rootCoef(gw, gdom);
                Kokkos::parallel_for("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA (int idx){
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
                        gw.coef(idom,7) += gdom.dt * qt * coef_wat(idx) * coef_root(idx) / gdom.dz(iGlob);
                        if (kk == 1)    {
                            gw.coef(idom,7) += gdom.dt * qe / gdom.dz(iGlob);
                        }
                });
        	}
            // Flux Source/Sink
            else if (sstype == 1)    {
                real qbc = interpolateLinear(ts, gdom.etime);
                Kokkos::parallel_for("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA (int idx){
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
                        gw.coef(idom,7) += gdom.dt * qbc;
                });
            }
            // Internal Drainage with Fixed Head
            // Note: For now, this only supports draining in the saturated zone
            else if (sstype == 2)   {
                // real hbc = interpolateLinear(ts, gdom.etime);
                real hbc = ssvals[0];
                if (direction == XPLUS || direction == XMINUS)	{
                	Kokkos::parallel_for("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA (int idx){
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        real flux = 0.0;
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
                  		// get index of cells next to the source/sink
                        int jp = idom+1, jm = idom-1, kp = idom+gdom.nx*gdom.ny, km = idom-gdom.nx*gdom.ny;
                        // treat the source/sink as a pressure boundary
                        gw.coef(jp,7) -= Cpipe * gw.coef(idom,3) * hbc;
                        gw.coef(jm,7) -= Cpipe * gw.coef(idom,4) * hbc;
                        gw.coef(kp,7) -= Cpipe * gw.coef(idom,5) * hbc;
                        gw.coef(km,7) -= Cpipe * gw.coef(idom,6) * hbc;
                        // exclude the source/sink cell from linear system
                        gw.coef(idom,1) = 0.0; gw.coef(idom,3) = 0.0; gw.coef(idom,5) = 0.0;
                        gw.coef(idom,2) = 0.0; gw.coef(idom,4) = 0.0; gw.coef(idom,6) = 0.0;
                        gw.coef(jp,4) = 0.0;	gw.coef(jm,3) = 0.0;
                        gw.coef(kp,6) = 0.0;	gw.coef(km,5) = 0.0;
                	});
                }
                else if (direction == YPLUS || direction == YMINUS)	{
                	Kokkos::parallel_for("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA (int idx){
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
                        int ip = idom+1, im = idom-1, kp = idom+gdom.nx*gdom.ny, km = idom-gdom.nx*gdom.ny;
                        // Seepage as a fixed H condition
                        gw.coef(ip,7) -= Cpipe * gw.coef(idom,1) * hbc;
                        gw.coef(im,7) -= Cpipe * gw.coef(idom,2) * hbc;
                        gw.coef(kp,7) -= Cpipe * gw.coef(idom,5) * hbc;
                        gw.coef(km,7) -= Cpipe * gw.coef(idom,6) * hbc;

                        gw.coef(idom,1) = 0.0; gw.coef(idom,3) = 0.0; gw.coef(idom,5) = 0.0;
                        gw.coef(idom,2) = 0.0; gw.coef(idom,4) = 0.0; gw.coef(idom,6) = 0.0;
                        gw.coef(ip,2) = 0.0;	gw.coef(im,1) = 0.0;
                        gw.coef(kp,6) = 0.0;	gw.coef(km,5) = 0.0;
                	});
                }
                else 	{
                	Kokkos::parallel_for("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA (int idx){
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        real flux = 0.0;
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
                        // get index of cells next to the source/sink
                        int ip = idom+1, im = idom-1, jp = idom+1, jm = idom-1;
                        // treat the source/sink as a pressure boundary
                        gw.coef(jp,7) -= Cpipe * gw.coef(idom,3) * hbc;
                        gw.coef(jm,7) -= Cpipe * gw.coef(idom,4) * hbc;
                        gw.coef(ip,7) -= Cpipe * gw.coef(idom,1) * hbc;
                        gw.coef(im,7) -= Cpipe * gw.coef(idom,2) * hbc;
                        // exclude the source/sink cell from linear system
                        gw.coef(idom,1) = 0.0; gw.coef(idom,3) = 0.0; gw.coef(idom,5) = 0.0;
                        gw.coef(idom,2) = 0.0; gw.coef(idom,4) = 0.0; gw.coef(idom,6) = 0.0;
                        gw.coef(jp,4) = 0.0;	gw.coef(jm,3) = 0.0;
                        gw.coef(ip,2) = 0.0;	gw.coef(im,1) = 0.0;
                	});
                }

            }
        }
    }

    inline void applyWCSS(GwState &gw, GwDomain &gdom) {
    	if (ncellsIT > 0)   {
			Qoutflow = 0.0;
			Qinflow = 0.0;
    		// ET (Penman-Monteith)
        	if (sstype == 0)	{
        		real qt = -interpolateLinear(tran, gdom.etime);
                real qe = -interpolateLinear(evap, gdom.etime);
                Kokkos::parallel_reduce("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA (int idx, real &tmp){
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
                        gw.wc(iGlob,1) += gdom.dt * qt * coef_wat(idx) * coef_root(idx) / gdom.dz(iGlob);
                        tmp += qt * coef_wat(idx) * coef_root(idx) * gdom.dx * gdom.dy;
                        if (kk == 1)    {
                            gw.wc(iGlob,1) += gdom.dt * qe / gdom.dz(iGlob);
                            tmp += qe * gdom.dx * gdom.dy;
                        }
				} , Kokkos::Sum<real>(Qoutflow) );
        	}
            // Flux Source/Sink
            else if (sstype == 1)    {
                real qbc = interpolateLinear(ts, gdom.etime);
                Kokkos::parallel_reduce("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA (int idx, real &tmp){
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
                        gw.wc(iGlob,1) += gdom.dt * qbc;
                        if (qbc > 0)	{tmp += qbc * gdom.dt;}
                        else {tmp -= qbc * gdom.dt;}
                }, Kokkos::Sum<real>(Qinflow) );
            }
            else if (sstype == 2)   {
                // real hbc = interpolateLinear(ts, gdom.etime);
                real hbc = ssvals[0], flux;
                if (direction == XPLUS || direction == XMINUS)	{
                	Kokkos::parallel_reduce("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA (int idx, real &tmp){
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
                        // Calculate the cumulative outflow
                        int jp = iGlob+gdom.nxhc, jm = iGlob-gdom.nxhc, kp = iGlob+gdom.nxhc*gdom.nyhc, km = iGlob-gdom.nxhc*gdom.nyhc;
                        if (hbc < gw.h(jp,1))	{
                        	tmp += Cpipe * gdom.dx * gdom.dz(iGlob) * gw.k(iGlob,1) * (hbc - gw.h(jp,1)) / gdom.dy;
                        }
                        if (hbc < gw.h(jm,1))	{
                        	tmp += Cpipe * gdom.dx * gdom.dz(iGlob) * gw.k(jm,1) * (hbc - gw.h(jm,1)) / gdom.dy;
                        }
                        if (hbc < gw.h(kp,1))	{
                        	tmp += Cpipe * gdom.dx * gdom.dy * gw.k(iGlob,2) * (hbc - gw.h(kp,1)) / gdom.dz(kp);
                        }
                        if (hbc < gw.h(km,1))	{
                        	tmp += Cpipe * gdom.dx * gdom.dy * gw.k(km,2) * (hbc - gw.h(km,1)) / gdom.dz(km);
                        }
                	}, Kokkos::Sum<real>(flux));
                }
                else if (direction == YPLUS || direction == YMINUS)	{
                	Kokkos::parallel_reduce("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA (int idx, real &tmp){
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
                        // Only consider drainage under fully saturated condition
                        int ip = iGlob+1, im = iGlob-1, kp = iGlob+gdom.nxhc*gdom.nyhc, km = iGlob-gdom.nxhc*gdom.nyhc;
                        if (hbc < gw.h(ip,1))	{
                        	tmp += Cpipe * gdom.dy * gdom.dz(iGlob) * gw.k(iGlob,0) * (hbc - gw.h(ip,1)) / gdom.dx;
                        }
                        if (hbc < gw.h(im,1))	{
                        	tmp += Cpipe * gdom.dy * gdom.dz(iGlob) * gw.k(im,0) * (hbc - gw.h(im,1)) / gdom.dx;
                        }
                        if (hbc < gw.h(kp,1))	{
                        	tmp += Cpipe * gdom.dx * gdom.dy * gw.k(iGlob,2) * (hbc - gw.h(kp,1)) / gdom.dz(kp);
                        }
                        if (hbc < gw.h(km,1))	{
                        	tmp += Cpipe * gdom.dx * gdom.dy * gw.k(km,2) * (hbc - gw.h(km,1)) / gdom.dz(km);
                        }
                	}, Kokkos::Sum<real>(flux));
                }
                else	{
                	Kokkos::parallel_reduce("gw_et", ncellsIT, KOKKOS_CLASS_LAMBDA (int idx, real &tmp){
                        int ii, jj, kk, ivg, idom, iGlob = icells[idx];
                        gdom.unpackIndicesHalo(iGlob, kk, jj, ii);
                        idom = (kk-hc)*gdom.nx*gdom.ny + (jj-hc)*gdom.nx + ii - hc;
                        // Only consider drainage under fully saturated condition
                        int ip = iGlob+1, im = iGlob-1, jp = iGlob+gdom.nxhc, jm = iGlob-gdom.nxhc;
                        if (hbc < gw.h(ip,1))	{
                        	tmp += Cpipe * gdom.dy * gdom.dz(iGlob) * gw.k(iGlob,0) * (hbc - gw.h(ip,1)) / gdom.dx;
                        }
                        if (hbc < gw.h(im,1))	{
                        	tmp += Cpipe * gdom.dy * gdom.dz(iGlob) * gw.k(im,0) * (hbc - gw.h(im,1)) / gdom.dx;
                        }
                        if (hbc < gw.h(jp,1))	{
                        	tmp += Cpipe * gdom.dx * gdom.dz(iGlob) * gw.k(iGlob,1) * (hbc - gw.h(jp,1)) / gdom.dy;
                        }
                        if (hbc < gw.h(jm,1))	{
                        	tmp += Cpipe * gdom.dx * gdom.dz(iGlob) * gw.k(jm,1) * (hbc - gw.h(jm,1)) / gdom.dy;
                        }
                	}, Kokkos::Sum<real>(flux));
                }
                Qoutflow = -flux;
            }

        }
    }



};
#endif


class SourceSink{
public:
	std::vector<std::string> id;
    SourceSinkData swss;
    #if SERGHEI_SUBSURFACE_MODEL
    std::vector<GwSS> gwss;
    #endif
};

#endif
