/*=========================================================================

  Program:   Visualization Toolkit
  Module:    PowerCrust.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __ConvertToPolydata_h
#define __ConvertToPolydata_h

#include "vtkPolyDataAlgorithm.h"

#include <vtkSmartPointer.h>
#include <vector>
#include <string>
class vtkPolyData;

class PowerCrust : public vtkPolyDataAlgorithm
{
  
  public:
   
    vtkTypeMacro(PowerCrust, vtkPolyDataAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // 
    static PowerCrust *New();

protected:
    PowerCrust();
    ~PowerCrust();

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int FillOutputPortInformation( int port, vtkInformation* info );

private:
    PowerCrust(const PowerCrust&);  // Not implemented.
    void operator=(const PowerCrust&);  // Not implemented.

};

#endif
