/*=========================================================================

  Program:   Principal Component Analysis
  Module:    vtkVectorPCA.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __PCA_H__
#define __PCA_H__

inline bool SameSign(float a, float b) {
    return a*b >= 0.0f;
}

inline void swap( int&a, int&b ){
  int swp=a; a=b; b=swp;
};


void GetInvariants(double* m, double* pqr);
int CalculateRoots(double* a, double* r);
void CalculateEigenvector(double* m, double lambda, double* evect);
void CalculateCross(double* v1, double* v2, double* res);
double NormSquared(double* v);
int getLargest(double* v);

#endif