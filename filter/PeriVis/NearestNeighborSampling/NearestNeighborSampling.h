/*=========================================================================

  Program:   Visualization Toolkit
  Module:    NearestNeighborSampling.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __NearestNeighborSampling_h
#define __NearestNeighborSampling_h

#include "vtkImageAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkIdList;
class vtkPoints;
class vtkPolyData;
class vtkDoubleArray;

class NearestNeighborSampling : public vtkImageAlgorithm
{
  
  public:
   
    vtkTypeMacro(NearestNeighborSampling, vtkImageAlgorithm);
    void PrintSelf(ostream &os, vtkIndent indent);
    
    static NearestNeighborSampling *New();

    // Description:
    // Specify i-j-k dimensions on which to sample input points.
    vtkGetVectorMacro(SampleDimensions,int,3);

    // Description:
    // Set the i-j-k dimensions on which to sample the distance function.
    void SetSampleDimensions(int i, int j, int k);

    // Description:
    // Set the i-j-k dimensions on which to sample the distance function.
    void SetSampleDimensions(int dim[3]);

   // Description:
  // These are basically a convenience method that calls SetInputArrayToProcess
  // to set the array used as the input scalars.  The fieldAssociation comes
  // from the vtkDataObject::FieldAssocations enum.  The fieldAttributeType
  // comes from the vtkDataSetAttributes::AttributeTypes enum.
  virtual void SetInputScalars(int fieldAssociation, const char *name);
  virtual void SetInputScalars(int fieldAssociation, int fieldAttributeType);

  protected:
    NearestNeighborSampling();
    ~NearestNeighborSampling();

    virtual int RequestInformation ( vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

	virtual int ProcessRequest(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector);

    // see algorithm for more info
    virtual int FillInputPortInformation(int port, vtkInformation* info);

    double ModelBounds[6];
    int SampleDimensions[3];
    //int FittingFunction;

  private:
    NearestNeighborSampling(const NearestNeighborSampling&);  // Not implemented.
    void operator=(const NearestNeighborSampling&);  // Not implemented.

};

#endif
;