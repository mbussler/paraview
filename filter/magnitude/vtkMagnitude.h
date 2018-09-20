/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMagnitude.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __Magnitude_h
#define __Magnitude_h

#include "vtkDataSetAlgorithm.h"

class vtkMagnitude : public vtkDataSetAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkMagnitude, vtkDataSetAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkMagnitude *New();
   

  protected:
    vtkMagnitude();
    ~vtkMagnitude();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

  private:
    vtkMagnitude(const vtkMagnitude&);  // Not implemented.
    void operator=(const vtkMagnitude&);  // Not implemented.
  
};

#endif
