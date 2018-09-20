/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkConvertToPolydata.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __ConvertToPolydata_h
#define __ConvertToPolydata_h

#include "vtkPolyDataAlgorithm.h"

class vtkConvertToPolydata : public vtkPolyDataAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkConvertToPolydata, vtkPolyDataAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkConvertToPolydata *New();
   

  protected:
    vtkConvertToPolydata();
    ~vtkConvertToPolydata();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int FillInputPortInformation(int port, vtkInformation *info);

  private:
    vtkConvertToPolydata(const vtkConvertToPolydata&);  // Not implemented.
    void operator=(const vtkConvertToPolydata&);  // Not implemented.
  
};

#endif
