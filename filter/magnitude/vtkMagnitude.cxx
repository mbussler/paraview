/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMagnitude.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkMagnitude.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkDataSet.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkTransform.h"
#include "vtkSmartPointer.h"

vtkStandardNewMacro(vtkMagnitude);

//==============================================================================
vtkMagnitude::vtkMagnitude()
{
  this->SetNumberOfInputPorts(1);
  
  // by default, process active point tensors
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                                                  vtkDataSetAttributes::TENSORS);


}
//==============================================================================
vtkMagnitude::~vtkMagnitude()
{
  // clean up
}

//==============================================================================
int vtkMagnitude::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkDataSet *input = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataSet *output = vtkDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDataArray *inData;
  double data[9];
  vtkIdType numPts, inPtId, ptIncr, i;
  int j;

  
  vtkDebugMacro(<<"Calculating magnitude");

  vtkPointData *outPD = output->GetPointData();
  inData = this->GetInputArrayToProcess(0, inputVector);
  numPts = input->GetNumberOfPoints();

  if ( !inData || numPts < 1 )
    {
    vtkErrorMacro(<<"No data!");
    return 1;
    }

  vtkSmartPointer<vtkDoubleArray> magn = vtkSmartPointer<vtkDoubleArray>::New();
  magn->SetName("magnitude");
  magn->SetNumberOfComponents(1);
  magn->SetNumberOfTuples(numPts);

  //
  // Traverse all Input points, calculate and store magnitude
  //

  int nNumComponents = inData->GetNumberOfComponents();
  
  for (inPtId=0; inPtId < numPts; inPtId++)
    {
    // Translation is postponed

    // compute magnitude of data as root from sum of squared
    double m = 0.0;
    for (i=0; i<nNumComponents; i++)
    {
      double v = inData->GetComponent(inPtId, i);
      m += v*v;
    }
    m = sqrt(m);
   
    // store to data array
    magn->SetValue(inPtId, m);
    
    }
  
  vtkDebugMacro(<<"Calculated " << numPts <<" Magnitude values");

  //
  // Update output and release memory
  //

  output->ShallowCopy(input);

  vtkPointData *pd = output->GetPointData();
  pd->AddArray(magn);
  
  output->Squeeze();

  return 1;
}

//==============================================================================
void vtkMagnitude::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
