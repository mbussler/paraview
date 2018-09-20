/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCellFlux.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __CellFlux_h
#define __CellFlux_h

#include "vtkDataSetAlgorithm.h"

class vtkCellFlux : public vtkDataSetAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkCellFlux, vtkDataSetAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkCellFlux *New();

  protected:
    vtkCellFlux();
    ~vtkCellFlux();

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  private:
    vtkCellFlux(const vtkCellFlux&);  // Not implemented.
    void operator=(const vtkCellFlux&);  // Not implemented.
    
    double calculateArea();
  
};

#endif
