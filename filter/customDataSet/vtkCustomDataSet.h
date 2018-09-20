/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCustomDataSet.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __CustomDataSet_h
#define __CustomDataSet_h

#include "vtkDataSetAlgorithm.h"

class vtkCustomDataSet : public vtkDataSetAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkCustomDataSet, vtkDataSetAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkCustomDataSet *New();
  
    // Description:
    // Get the input to the filter.
    vtkDataObject *GetInput();
    
    // Description:
    // Set the flag indicating that Eigenvectors are possibly flipped during evaluation.
    vtkSetMacro(flipEigenVectorsOnEval,int);
    vtkGetMacro(flipEigenVectorsOnEval,int);
    vtkBooleanMacro(flipEigenVectorsOnEval,int);

  protected:
    vtkCustomDataSet();
    ~vtkCustomDataSet();

    virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *); //generate output data
    virtual int RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    virtual int FillInputPortInformation(int port, vtkInformation *info);
    virtual int RequestDataObject(vtkInformation *, vtkInformationVector **,
                                  vtkInformationVector *);

    // flag to indicate that Eigenvectors are possibly flipped during evaluation.
    int flipEigenVectorsOnEval;
    
  private:
    vtkCustomDataSet(const vtkCustomDataSet&);  // Not implemented.
    void operator=(const vtkCustomDataSet&);  // Not implemented.
  
};

#endif
