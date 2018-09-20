/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkAddTime.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __AddTime_h
#define __AddTime_h

#include "vtkDataSetAlgorithm.h"

class vtkAddTime : public vtkDataSetAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkAddTime, vtkDataSetAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkAddTime *New();
   

  protected:
    vtkAddTime();
    ~vtkAddTime();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

  private:
    vtkAddTime(const vtkAddTime&);  // Not implemented.
    void operator=(const vtkAddTime&);  // Not implemented.
  
};

#endif
