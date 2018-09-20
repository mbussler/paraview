/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRidgt.h

  Copyright (c) Oliver Fernandes

=========================================================================*/

#ifndef __vtkRidge_h
#define __vtkRidge_h

#include "vtkPolyDataAlgorithm.h" //superclass
#include "vtkSmartPointer.h" // compiler errors if this is forward declared

class vtkPolyData;
class vtkTransform;
class vtkInformation;
class vtkInformationVector;
class vtkIterativeClosestPointTransform;

typedef unsigned int uint;

#define DUMPVECTOR(x) (cout << "(" << x[0] << "," << x[1] << "," << x[2] << ")" << endl)

class vtkRidge : public vtkPolyDataAlgorithm
{
 public:

  vtkTypeMacro(vtkRidge, vtkPolyDataAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent);

  static vtkRidge *New();

  vtkGetMacro( isoThreshold, double);
  vtkSetMacro( isoThreshold, double);


protected:
  vtkRidge();
  ~vtkRidge();

  uint dim[3];
  double isoThreshold;

  // Make sure the pipeline knows what type we expect as input
  int FillInputPortInformation( int port, vtkInformation* info );
  // Generate output
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *); //the function that makes this class work with the vtk pipeline
};
#endif
