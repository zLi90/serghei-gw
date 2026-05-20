//-----------------------------------------------------------------------------
// Reference: EPA SWMM5 (Storm Water Management Model) Version 5.2
// Source: https://github.com/USEPA/Stormwater-Management-Model/blob/develop/src/solver/xsect.c
// Note: This file has been simplified for our use case.
// We have omitted considerations for subcatchments, unit conversions, file outputs, etc.
//-----------------------------------------------------------------------------
#ifndef XSECT_H
#define XSECT_H
#include <cmath>
#include <iostream>
#include <string>
#include "../define.h"
#include "DrainageState.h"


class XXsect{
    public:

     inline int xsect_isOpen(int type)
    {
        if (type == CIRCULAR)
        {
            return 0;
        }
        return 1;
    }

inline double xsect_getAofY(int j, double y, Link& Tlink)
{

        double r = 0.5*Tlink.yFull(j);
        if (y >= Tlink.yFull(j)) return PI*r*r;
        if (y <= 0.0) return 0.0;
        if (y < r)
        {
            double alpha = acos((r - y) / r);
            return alpha * r * r - (r - y) * r * sin(alpha);
        }
        if (y > r){
            double alpha = acos((y - r) / r);
            return PI * r * r - alpha * r * r + (y - r) * r * sin(alpha);          
        }

        
        return PI * r * r / 2.0;

}

    inline double xsect_getWofY(int j, double y, Link& Tlink)
//
//  Input:   xsect = ptr. to a cross section data structure
//           y = depth ft)
//  Output:  returns top width (ft)
//  Purpose: computes xsection's top width at a given depth.
//
{
   
    double r = 0.5*Tlink.yFull(j);
    if (y <= 0.0) return 0.0;
    if (y >= Tlink.yFull(j)) return Tlink.wMax(j);
    if (y < r)
    {
        double alpha = acos((r - y) / r);
        return 2.0 * r * sin(alpha);
    }
    if (y > r)
    {
        double alpha = acos((y - r) / r);
        return 2.0 * r * sin(alpha);
    }
 
    return 2.0 * r;
  
}

    inline double xsect_getRofY(int j, double y, Link& Tlink)
{
    double yNorm = y / Tlink.yFull(j); 
    double area = xsect_getAofY(j, y, Tlink);
    double r = 0.5*Tlink.yFull(j);
    if (y >= Tlink.yFull(j)) return area / (2.0*PI*r);
    if (y <= 0.0) return 0.0;
    if (y < r)
    {
        double alpha = acos((r - y) / r);
        return area / (2*alpha*r);
    };
    if (y > r)
    {
        double alpha = acos((y - r) / r);
        return area / (2*PI*r-2*alpha*r);
    };

    return area / (PI*r);
    }

double circ_getYofA(int j, double a, Link& Tlink)
{
     double r = 0.5*Tlink.yFull(j);
     double y_low, y_mid, y_high, area_mid;
     int max_iter = 1000;
     double tol = 1.0e-3;

     y_low = 0.0;            
     y_high = 2 * r; 


    if (a <= 0)
        return 0.0;
    if (a == PI * r*r)
        return 2.0*r;


    for (int i; i < max_iter; i++)
    {
        y_mid = (y_low + y_high) / 2;
        area_mid = xsect_getAofY(j, y_mid, Tlink);

        if (abs(area_mid - a) < tol)
            return y_mid;
        else if (area_mid < a)
            y_low = y_mid;
        else y_high = y_mid;
    }
    return 0.0;
}

double getYcircular(double alpha)
{
    double theta;
    if ( alpha >= 1.0 ) return 1.0;
    if ( alpha <= 0.0 ) return 0.0;
    if ( alpha <= 1.0e-5 )
    {
        theta = pow(37.6911*alpha, 1./3.);
        return theta * theta / 16.0;
    }
    theta = getThetaOfAlpha(alpha);
    return (1.0 - cos(theta/2.)) / 2.0;
}

double getThetaOfAlpha(double alpha)
{
    int    k;
    double theta, theta1, ap, d;

    if ( alpha > 0.04 ) theta = 1.2 + 5.08 * (alpha - 0.04) / 0.96;
    else theta = 0.031715 - 12.79384 * alpha + 8.28479 * sqrt(alpha);
    theta1 = theta;
    ap  = (2.0*PI) * alpha;
    for (k = 1; k <= 40; k++ )
    {
        d = - (ap - theta + sin(theta)) / (1.0 - cos(theta));
        // --- modification to improve convergence for large theta
        if ( d > 1.0 ) d = sign( 1.0, d );
        theta = theta - d;
        if ( fabs(d) <= 0.0001 ) return theta;
    }
    return theta1;
}
inline double xsect_getAofS(int j, double s, Link& Tlink) 
//
//  Input:   xsect = ptr. to a cross section data structure
//           s = section factor (ft^(8/3))
//  Output:  returns area (ft2)
//  Purpose: computes xsection's area at a given section factor.
//
{
    double psi = s / Tlink.sFull(j);
    if ( s <= 0.0 ) return 0.0;
    if ( s > Tlink.sMax(j) ) s = Tlink.sMax(j);

    return circ_getAofS(j, s, Tlink);
}  

inline double circ_getAofS(int j, double s, Link& Tlink)
{

    double psi = s / Tlink.sFull(j);
    double x_low, x_high, x_mid, area_mid;
    int max_iter = 1000;
    double tol = 1.0e-2;
    if (psi == 0.0) return 0.0;
    if (psi >= 1.0) return Tlink.aFull(j);

    // --- use special function for small s/sFull
    if (psi <= 0.015) return Tlink.aFull(j) * getAcircular(psi);

     x_low = 0.0;            
     x_high = 1; 

    for (int i; i < max_iter; i++)
    {
        x_mid = (x_low + x_high) / 2;
        area_mid = -0.89906279*x_mid*x_mid*x_mid + 1.58653555*x_mid*x_mid + 0.4134494*x_mid - 0.00288671;

        if (abs(area_mid - psi) < tol)
            return Tlink.aFull(j) * x_mid;
        else if (area_mid < psi)
            x_low = x_mid;
        else x_high = x_mid;
    }
    return 0.0;
}

inline double getAcircular(double psi)
{
    double theta;
    if ( psi >= 1.0 ) return 1.0;
    if ( psi <= 0.0 ) return 0.0;
    if ( psi <= 1.0e-6 )
    {
        theta = pow(124.4797*psi, 3./13.);
        return theta*theta*theta / 37.6911;
    }
    theta = getThetaOfPsi(psi);
    return (theta - sin(theta)) / (2.0 * PI);
}

inline double getThetaOfPsi(double psi)
{
    int    k;
    double theta, theta1, ap, tt, tt23, t3, d;

    if      (psi > 0.90)  theta = 4.17 + 1.12 * (psi - 0.90) / 0.176;
    else if (psi > 0.5)   theta = 3.14 + 1.03 * (psi - 0.5) / 0.4;
    else if (psi > 0.015) theta = 1.2 + 1.94 * (psi - 0.015) / 0.485;
    else                  theta = 0.12103 - 55.5075 * psi +
                                  15.62254 * sqrt(psi);
    theta1 = theta;
    ap     = (2.0*PI) * psi;

    for (k = 1; k <= 40; k++)
    {
        theta    = fabs(theta);
        tt       = theta - sin(theta);
        tt23     = pow(tt, 2./3.);
        t3       = pow(theta, 1./3.);
        d        = ap * theta / t3 - tt * tt23;
        d        = d / ( ap*(2./3.)/t3 - (5./3.)*tt23*(1.0-cos(theta)) );
        theta    = theta - d;
        if ( fabs(d) <= 0.0001 ) return theta;
    }
    return theta1;
}



};

#endif
