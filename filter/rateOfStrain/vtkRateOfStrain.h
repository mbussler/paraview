/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRateOfStrain.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __EigenAnalysis_h
#define __EigenAnalysis_h

#include "vtkDataSetAlgorithm.h"

class vtkRateOfStrain : public vtkDataSetAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkRateOfStrain, vtkDataSetAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkRateOfStrain *New();
   

  protected:
    vtkRateOfStrain();
    ~vtkRateOfStrain();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

  private:
    vtkRateOfStrain(const vtkRateOfStrain&);  // Not implemented.
    void operator=(const vtkRateOfStrain&);  // Not implemented.
  
};

#endif
