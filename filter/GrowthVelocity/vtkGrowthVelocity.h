/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGrowthVelocity.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __GrowthVelocity_h
#define __GrowthVelocity_h

#include "vtkDataSetAlgorithm.h"

class vtkGrowthVelocity : public vtkDataSetAlgorithm
{
  
  public:
   
   vtkSetMacro(Time, double);
   vtkGetMacro(Time, double);
   vtkSetMacro(Range, double);
   vtkGetMacro(Range, double);
   vtkSetMacro(Buckets, int);
   vtkGetMacro(Buckets, int);
   
    vtkTypeMacro(vtkGrowthVelocity, vtkDataSetAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static vtkGrowthVelocity *New();
   

  protected:
    vtkGrowthVelocity();
    ~vtkGrowthVelocity();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

  private:
    vtkGrowthVelocity(const vtkGrowthVelocity&);  // Not implemented.
    void operator=(const vtkGrowthVelocity&);  // Not implemented.

    virtual int FillOutputPortInformation(int port, vtkInformation* info);


  double Time;
  double Range;
  int Buckets;
  
};

#endif
