/*=========================================================================

  Program:   Visualization Toolkit
  Module:    VonMiesesStress.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __VonMiesesStress_h
#define __VonMiesesStress_h

#include "vtkDataSetAlgorithm.h"

class VonMiesesStress : public vtkDataSetAlgorithm
{
  
  public:
   
    vtkTypeMacro(VonMiesesStress, vtkDataSetAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static VonMiesesStress *New();
   

  protected:
    VonMiesesStress();
    ~VonMiesesStress();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

  private:
    VonMiesesStress(const VonMiesesStress&);  // Not implemented.
    void operator=(const VonMiesesStress&);  // Not implemented.
  
};

#endif
