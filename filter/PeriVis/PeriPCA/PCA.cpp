/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkVectorPCA.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "PCA.h"
#include <math.h>
#include <stdlib.h>
#include <cmath>

#include "linalg.h"



//==============================================================================
void GetInvariants(double* m, double* pqr) 
{
    pqr[0] = -m[0*3+0]*m[1*3+1]*m[2*3+2]
    -m[0*3+1]*m[1*3+2]*m[2*3+0]
    -m[0*3+2]*m[1*3+0]*m[2*3+1]
    +m[2*3+0]*m[1*3+1]*m[0*3+2]
    +m[2*3+1]*m[1*3+2]*m[0*3+0]
    +m[2*3+2]*m[1*3+0]*m[0*3+1];


    pqr[1] = m[1*3+1]*m[2*3+2] - m[2*3+1]*m[1*3+2]
    + m[0*3+0]*m[1*3+1] - m[1*3+0]*m[0*3+1]
    + m[0*3+0]*m[2*3+2] - m[0*3+2]*m[2*3+0];

    pqr[2] = -(m[0*3+0]+m[1*3+1]+m[2*3+2]);
}


int CalculateRoots(double* a, double* r) 
{
    double c1 = a[1] - a[2]*a[2]/3.;
    double c0 = a[0] - a[1]*a[2]/3. + 2./27.*a[2]*a[2]*a[2];
    // Make cubic coefficient 4 and linear coefficient +- 3
    // by substituting y = z*k and multiplying with 4/k^3
    if (c1 == 0) {
        if (c0 == 0)     r[0] =  0;
        else if (c0 > 0) r[0] = -pow(c0, 1./3.);
        else             r[0] =  pow(-c0, 1./3.);
    }
    else {
        bool negc1 = c1 < 0;
        double absc1 = negc1 ? -c1 : c1;
        double k = sqrt(4./3.*absc1);
        double d0 = c0 * 4./(k*k*k);
        // Find the first solution
        if (negc1) {
            if (d0 > 1)       r[0] = -cosh(acosh(d0)/3);
            else if (d0 > -1) r[0] = -cos(acos(d0)/3);
            else              r[0] =  cosh(acosh(-d0)/3);
        }
        else {
            r[0] = -sinh(asinh(d0)/3);
        }
        // Transform back
        r[0] *= k;
    }
    r[0] -= a[2]/3;
    // Other two solutions
    double p = r[0] + a[2];
    double q = r[0]*p + a[1];
    double discrim = p*p - 4*q;
    //    if (forceReal && discrim < 0.0) discrim = 0.0;
    if (discrim >= 0) {
        double root = sqrt(discrim);
        r[1] = (-p - root)/2.;
        r[2] = (-p + root)/2.;
        return 3;
    }
    else {
        double root = sqrt(-discrim);
        r[1] = -p/2;
        r[2] = root/2.;
        return 1;
    }    
}


void CalculateEigenvector(double* m, double lambda, double* evect) {
    double norm[3];
    double cross[9];
    double red[9];

    red[0] =  m[0]-lambda;
    red[1] =  m[1];
    red[2] =  m[2];
    red[3] =  m[3];
    red[4] =  m[4]-lambda;
    red[5] =  m[5];
    red[6] =  m[6];
    red[7] =  m[7];
    red[8] =  m[8]-lambda;

    CalculateCross(&red[3], &red[6], &cross[0]);
    CalculateCross(&red[6], &red[0], &cross[3]);
    CalculateCross(&red[0], &red[3], &cross[6]);

    norm[0] = NormSquared(&cross[0]);
    norm[1] = NormSquared(&cross[3]);
    norm[2] = NormSquared(&cross[6]);
    
    int best = getLargest(norm);
    double len = sqrt(norm[best]);

    if (len > 0) {
        evect[0] = cross[best*3+0] / len;
        evect[1] = cross[best*3+1] / len;
        evect[2] = cross[best*3+2] / len;
    }
    else {
        evect[0] = 0.0;
        evect[1] = 0.0;
        evect[2] = 0.0;
    }
}


void CalculateCross(double* v1, double* v2, double* res) {
    res[0] = v1[1]*v2[2] - v1[2]*v2[1];
    res[1] = v1[2]*v2[0] - v1[0]*v2[2];
    res[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

double NormSquared(double* v) {
    return (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}


int getLargest(double* v) 
{
    if (v[0]>v[1]) {
		if (v[0]>v[2]) {
			return 0;
		} else {
			return 2;
		}
	} else {
		if (v[1]>v[2]) {
			return 1;
		} else {
			return 2;
		}
    }
}
