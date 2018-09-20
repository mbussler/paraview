/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkVectorCalc.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __VectorCalc_h
#define __VectorCalc_h

#include "vtkUnstructuredGridAlgorithm.h"

class vtkVectorCalc : public vtkUnstructuredGridAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkVectorCalc, vtkUnstructuredGridAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkVectorCalc *New();
   
    // Description:
    // Turn on/off scalar calculation of absolute values. 
    vtkSetMacro(AbsoluteValue,int);
    vtkGetMacro(AbsoluteValue,int);
    vtkBooleanMacro(AbsoluteValue,int);

    // Description:
    // Select function to calculate on given input data
    vtkSetMacro(Function,int);
    vtkGetMacro(Function,int);

  protected:
    vtkVectorCalc();
    ~vtkVectorCalc();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

    int Function; // Function to calculate on given input data
    int AbsoluteValue; // Boolean controls whether to calculate absolute values

  private:
    vtkVectorCalc(const vtkVectorCalc&);  // Not implemented.
    void operator=(const vtkVectorCalc&);  // Not implemented.
  
};

#endif
