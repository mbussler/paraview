/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkVectorPCA.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __VectorPCA_h
#define __VectorPCA_h

#include "vtkDataSetAlgorithm.h"

class vtkVectorPCA : public vtkDataSetAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkVectorPCA, vtkDataSetAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkVectorPCA *New();

  protected:
    vtkVectorPCA();
    ~vtkVectorPCA();

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  private:
    vtkVectorPCA(const vtkVectorPCA&);  // Not implemented.
    void operator=(const vtkVectorPCA&);  // Not implemented.

    void GetInvariants(double* m, double* pqr); 
    int CalculateRoots(double* a, double* r);
    void CalculateEigenvector(double* m, double lambda, double* evect);
    int getLargest(double* v);
    double NormSquared(double* v);
    void CalculateCross(double* v1, double* v2, double* res);


  };

#endif
