#ifndef _SWSOURCESINK_H_
#define _SWSOURCESINK_H_

#include "define.h"
#include "SArray.h"
#include "State.h"

#define INF_NONE 0
#define INF_CONSTANT 1
#define INF_HORTON 2
#define INF_GREENAMPT 3

//class State;	// forward declaration

/*
  TimeSeries provides a construct/class to store time series.
 */
class TimeSeries{

public:

  int np;     // number of points in time
  int nx = 1; // number of partitions in x direction
  int ny = 1; // number of partitions in y direction

  realArr time;
  realArr value;
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
  // TODO need to program GreenAmpt model
KOKKOS_INLINE_FUNCTION  real greenAmpt (const int ii,const real ) const{
    real infCap = 0.;
    return(infCap);
  }

  realArr infTime;

public:
    int model = -999;
    int nLabels = 0;
    realArr constCap ;
    // Horton
    realArr k;
    realArr fc;
    realArr f0;
    // Green-Ampt
    real ks = -999;
    real psi = -999;
    real dtheta = -999;

    real infDry = 1E-8;  // [L] threshold to consider dry for infiltration purposes

    // This is not a state variable
    // which is why it is here and not in class State.
    // It is necessary for output
    // and because it is a variable in the GreenAmpt model
    realArr infVol; // accumulated infiltration volume
    realArr rate;		// infiltration rate
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
         //			capacity = &InfiltrationModel::greenAmpt;
		 if(par.masterproc){
         	std::cerr << BDASH << "Green-Ampt infiltration capacity" << std::endl;
            std::cerr << RERROR << "Not enabled yet" << std::endl;
		 }
         error++;
         if(ks < 0){
		 	if(par.masterproc){
            	std::cerr << RERROR << "Saturated hydraulic conductivity not found for Green-Ampt infiltration model" << std::endl;
			}
            error++;
         }
         if(psi < 0){
		 	if(par.masterproc){
          		std::cerr << RERROR << "Average suction head not found for Green-Ampt infiltration model" << std::endl;
			}
            error++;
         }
         if(dtheta < 0){
			if(par.masterproc){
               	std::cerr << RERROR << "Water content difference not found for Green-Ampt infiltration model" << std::endl;
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



    inline void ComputeInfiltrationCapacity(const Domain &dom){
      #if SERGHEI_DEBUG_WORKFLOW
        std::cerr << GGD << __PRETTY_FUNCTION__ << std::endl;
      #endif
        if(model){
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
                        break;
                    }
            }
    }
};


class SourceSinkData{

public:

  TimeSeries rain;
  InfiltrationModel inf;
  realArr rainRate;

  //real timerRainInf=0;

  void allocate (Domain const &dom){
    if (dom.isRain)
      rainRate  = realArr ("rainRate", dom.nCellMem);

    if (inf.model)
      inf.allocate(dom);
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

      dom.unpackIndices (iGlob, iy, ix);
      int ii = dom.getHaloExtension(ix,iy);

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

  inline void ComputeSWSourceSink(const State &state, const Domain &dom){
    Kokkos::Timer timer;
    #if SERGHEI_DEBUG_WORKFLOW
    std::cerr << GGD << __PRETTY_FUNCTION__ << std::endl;
    #endif
    ComputeRain(dom);
    inf.ComputeInfiltrationCapacity(dom);
    //no rate correction is necessary here beacuse the rate correction is done in ComputeNewState, according to the new water depth
   // timerRainInf += timer.seconds();
   dom.timers.raininf += timer.seconds();
 }

};

#endif
