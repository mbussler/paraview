/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTensorNorm.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __EigenAnalysis_h
#define __EigenAnalysis_h

#include "vtkDataSetAlgorithm.h"

class vtkTensorNorm : public vtkDataSetAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkTensorNorm, vtkDataSetAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkTensorNorm *New();
   

  protected:
    vtkTensorNorm();
    ~vtkTensorNorm();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

  private:
    vtkTensorNorm(const vtkTensorNorm&);  // Not implemented.
    void operator=(const vtkTensorNorm&);  // Not implemented.
  
};

#endif
