/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkEigenAnalysis.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __EigenAnalysis_h
#define __EigenAnalysis_h

#include "vtkUnstructuredGridAlgorithm.h"

class vtkEigenAnalysis : public vtkUnstructuredGridAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkEigenAnalysis, vtkUnstructuredGridAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkEigenAnalysis *New();
   

  protected:
    vtkEigenAnalysis();
    ~vtkEigenAnalysis();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

  private:
    vtkEigenAnalysis(const vtkEigenAnalysis&);  // Not implemented.
    void operator=(const vtkEigenAnalysis&);  // Not implemented.
  
};

#endif
