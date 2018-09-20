/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkVector2Dto3D.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __Vector2Dto3D_h
#define __Vector2Dto3D_h

#include "vtkDataSetAlgorithm.h"

class vtkVector2Dto3D : public vtkDataSetAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkVector2Dto3D, vtkDataSetAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkVector2Dto3D *New();
   

  protected:
    vtkVector2Dto3D();
    ~vtkVector2Dto3D();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

  private:
    vtkVector2Dto3D(const vtkVector2Dto3D&);  // Not implemented.
    void operator=(const vtkVector2Dto3D&);  // Not implemented.
  
};

#endif
