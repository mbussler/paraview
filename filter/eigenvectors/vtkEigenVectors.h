/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkEigenVectors.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __EigenVectors_h
#define __EigenVectors_h

#include "vtkDataSetAlgorithm.h"

class vtkEigenVectors : public vtkDataSetAlgorithm
{
  
  public:
   
    vtkTypeMacro(vtkEigenVectors, vtkDataSetAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkEigenVectors *New();

    vtkSetMacro(EigenvalueMethod, int);
    vtkGetMacro(EigenvalueMethod, int);

    // Description
    // Turn on/off vorticity computation at streamline points
    // (necessary for generating proper stream-ribbons using the
    // vtkRibbonFilter.
    vtkSetMacro(Output3EV, bool);
    vtkGetMacro(Output3EV, bool);

    vtkSetMacro(FixEVs, bool);
    vtkGetMacro(FixEVs, bool);
      
    vtkSetMacro(AnisotrophyMeasure, bool);
    vtkGetMacro(AnisotrophyMeasure, bool);
    
  protected:
    vtkEigenVectors();
    ~vtkEigenVectors();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

  private:
    vtkEigenVectors(const vtkEigenVectors&);  // Not implemented.
    void operator=(const vtkEigenVectors&);  // Not implemented.

    int EigenvalueMethod;
    bool Output3EV;
    bool FixEVs;
    bool AnisotrophyMeasure;
};

#endif
