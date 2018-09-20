/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkDivergence.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __Divergence_h
#define __Divergence_h

#include "vtkDataSetAlgorithm.h"

class vtkDivergence : public vtkDataSetAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkDivergence, vtkDataSetAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkDivergence *New();
   
    // Description:
    // Turn on/off scalar calculation of absolute values. 
    vtkSetMacro(AbsoluteValue,int);
    vtkGetMacro(AbsoluteValue,int);
    vtkBooleanMacro(AbsoluteValue,int);

  protected:
    vtkDivergence();
    ~vtkDivergence();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

    int AbsoluteValue; // Boolean controls whether to calculate absolute values

  private:
    vtkDivergence(const vtkDivergence&);  // Not implemented.
    void operator=(const vtkDivergence&);  // Not implemented.
  
};

#endif
