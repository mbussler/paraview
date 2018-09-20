/*=========================================================================

  Program:   Visualization Toolkit
  Module:    PeriStress.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __PeriStress_h
#define __PeriStress_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkIdList;
class vtkPolyData;

class PeriStress : public vtkPolyDataAlgorithm
{
  
  public:
   
    vtkTypeMacro(PeriStress, vtkPolyDataAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static PeriStress *New();

  protected:
    PeriStress();
    ~PeriStress();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

	virtual int ProcessRequest(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector);

  private:
    PeriStress(const PeriStress&);  // Not implemented.
    void operator=(const PeriStress&);  // Not implemented.

	vtkSmartPointer<vtkIdList> 
		GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int id);

};

#endif
