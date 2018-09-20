/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkDivergence.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkDivergence.h"
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

#include <stdlib.h>

vtkStandardNewMacro(vtkDivergence);

//==============================================================================
vtkDivergence::vtkDivergence()
{
  this->SetNumberOfInputPorts(1);
  
  // by default, process active point tensors
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::TENSORS);

  this->AbsoluteValue = 0;
}
//==============================================================================
vtkDivergence::~vtkDivergence()
{
  // clean up
}

//==============================================================================
int vtkDivergence::RequestData(
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

  vtkDataArray *inTensors;
  double jacobian[9];
  vtkIdType numPts, inPtId, ptIncr, i;
  int j;

  
  vtkDebugMacro(<<"Calculating Divergence");

  vtkPointData *outPD = output->GetPointData();
  inTensors = this->GetInputArrayToProcess(0, inputVector);
  numPts = input->GetNumberOfPoints();

  if ( !inTensors || numPts < 1 )
    {
    vtkErrorMacro(<<"No data!");
    return 1;
    }

  vtkSmartPointer<vtkDoubleArray> div = vtkSmartPointer<vtkDoubleArray>::New();
  div->SetName("divergence");
  div->SetNumberOfComponents(1);
  div->SetNumberOfTuples(numPts);

  //
  // Traverse all Input points, calculate and store divergence
  //

  for (inPtId=0; inPtId < numPts; inPtId++)
    {
    // Translation is postponed

    inTensors->GetTuple(inPtId, jacobian);

    // compute divergence as trace of Jacobian
    double trace = 0.0;
    trace += jacobian[0];
    trace += jacobian[4];
    trace += jacobian[8];

    if( this->AbsoluteValue )
    {
      trace = abs(trace);
    }

    // store to data array
    div->SetValue(inPtId, trace);
    
    }
  
  vtkDebugMacro(<<"Calculated " << numPts <<" Divergence values");

  //
  // Update output and release memory
  //

  output->ShallowCopy(input);

  vtkPointData *pd = output->GetPointData();
  pd->AddArray(div);
  
  output->Squeeze();

  return 1;
}

//==============================================================================
void vtkDivergence::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
