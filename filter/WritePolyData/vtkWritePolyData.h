/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkWritePolyData.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __WritePolyData_h
#define __WritePolyData_h

#include "vtkPolyDataAlgorithm.h"

class vtkWritePolyData : public vtkPolyDataAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkWritePolyData, vtkPolyDataAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkWritePolyData *New();
   
    vtkGetStringMacro(OutputFileName);
    vtkSetStringMacro(OutputFileName);
    
  protected:
    vtkWritePolyData();
    ~vtkWritePolyData();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

    // The filename to write to
    char *OutputFileName;

  private:
    vtkWritePolyData(const vtkWritePolyData&);  // Not implemented.
    void operator=(const vtkWritePolyData&);  // Not implemented.

    virtual int RequestUpdateExtent(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector);


};

#endif
