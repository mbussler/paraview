/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkLightScattering.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __LightScattering_h
#define __LightScattering_h

#include "vtkPolyDataAlgorithm.h"
#include "linalg.h"
#include "vtkSmartPointer.h"

class vtkIdList;
class vtkPolyData;

class vtkLightScattering : public vtkPolyDataAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkLightScattering, vtkPolyDataAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkLightScattering *New();

	vtkSetMacro( Diffusion, double);
	vtkSetMacro( Stepsize,  double);
	vtkSetMacro( SeedEvery, int);
	vtkSetMacro( SeedParticleStride, int);

	vtkGetMacro( Diffusion, double);
	vtkGetMacro( Stepsize,  double);
	vtkGetMacro( SeedEvery, int);
	vtkGetMacro( SeedParticleStride, int);

	vtkSetMacro( ResetIntensity, bool);
	vtkGetMacro( ResetIntensity, bool);
	vtkSetMacro( InitialSeed, bool);
	vtkGetMacro( InitialSeed, bool);

	vtkSetClampMacro( Iterations, int, 1, 1000);
	vtkGetMacro( Iterations, int);

  protected:
    vtkLightScattering();
    ~vtkLightScattering();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

	vtkSmartPointer<vtkIdList> 
		GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int id);

  private:
    vtkLightScattering(const vtkLightScattering&);  // Not implemented.
    void operator=(const vtkLightScattering&);  // Not implemented.

	virtual int ProcessRequest(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector);

	vtkSmartPointer<vtkIdList> *m_conns;

	double Diffusion;
	double Stepsize;
	int    Iterations;
	int    SeedEvery;
	int    SeedParticleStride;
	bool   ResetIntensity;
	bool   InitialSeed;

};

#endif
