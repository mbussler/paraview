/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCellSize.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __CellSize_h
#define __CellSize_h

#include "vtkDataSetAlgorithm.h"

class vtkCellSize : public vtkDataSetAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkCellSize, vtkDataSetAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkCellSize *New();
   

  protected:
    vtkCellSize();
    ~vtkCellSize();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

  private:
    vtkCellSize(const vtkCellSize&);  // Not implemented.
    void operator=(const vtkCellSize&);  // Not implemented.
  
};

#endif
