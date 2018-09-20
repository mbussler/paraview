/*=========================================================================

  Program:   Visualization Toolkit
  Module:    PeriPCA.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __PeriPCA_h
#define __PeriPCA_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkIdList;
class vtkPolyData;

class PeriPCA : public vtkPolyDataAlgorithm
{
  
  public:
   
    vtkTypeMacro(PeriPCA, vtkPolyDataAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static PeriPCA *New();

  protected:
    PeriPCA();
    ~PeriPCA();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

	virtual int ProcessRequest(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector);

  private:
    PeriPCA(const PeriPCA&);  // Not implemented.
    void operator=(const PeriPCA&);  // Not implemented.

	vtkSmartPointer<vtkIdList> 
		GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int id);
	void calculatePCA( vtkPoints* points, double eigenValues[3], double* ev=NULL, double* covariance=NULL );

	void GetInvariants(double* m, double* pqr);
	int CalculateRoots(double* a, double* r);
	void CalculateEigenvector(double* m, double lambda, double* evect);
	void CalculateCross(double* v1, double* v2, double* res);
	double NormSquared(double* v);
	int getLargest(double* v);
};

#endif
