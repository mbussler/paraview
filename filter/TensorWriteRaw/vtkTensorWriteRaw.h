/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTensorWriteRaw.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __EigenAnalysis_h
#define __EigenAnalysis_h

#include "vtkImageAlgorithm.h"

class vtkTensorWriteRaw : public vtkImageAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkTensorWriteRaw, vtkImageAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkTensorWriteRaw *New();
   
    vtkGetStringMacro(DataFileName);
    vtkSetStringMacro(DataFileName);
    
  protected:
    vtkTensorWriteRaw();
    ~vtkTensorWriteRaw();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

    // The filename to write to
    char *DataFileName;

  private:
    vtkTensorWriteRaw(const vtkTensorWriteRaw&);  // Not implemented.
    void operator=(const vtkTensorWriteRaw&);  // Not implemented.


};

#endif
