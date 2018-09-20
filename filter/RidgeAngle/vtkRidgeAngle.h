/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRidgeAngle.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __RidgeAngle_h
#define __RidgeAngle_h

#include "vtkPolyDataAlgorithm.h"

class vtkRidgeAngle : public vtkPolyDataAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkRidgeAngle, vtkPolyDataAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkRidgeAngle *New();
   

  protected:
    vtkRidgeAngle();
    ~vtkRidgeAngle();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

  private:
    vtkRidgeAngle(const vtkRidgeAngle&);  // Not implemented.
    void operator=(const vtkRidgeAngle&);  // Not implemented.
  
};

#endif
