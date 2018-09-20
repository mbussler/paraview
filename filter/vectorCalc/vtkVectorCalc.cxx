/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkVectorCalc.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkVectorCalc.h"
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

vtkStandardNewMacro(vtkVectorCalc);

//==============================================================================
vtkVectorCalc::vtkVectorCalc()
{
  this->SetNumberOfInputPorts(1);
  
  // by default, process active point tensors
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::TENSORS);

  // by default, process active point scalars
  this->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::SCALARS);

  this->AbsoluteValue = 0;
}
//==============================================================================
vtkVectorCalc::~vtkVectorCalc()
{
  // clean up
}

//==============================================================================
int vtkVectorCalc::RequestData(
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
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDataArray *inTensors;
  double tensor[9];
  vtkIdType numPts, inPtId, ptIncr, i;
  int j;

  
  vtkPointData *outPD = output->GetPointData();
  inTensors = this->GetInputArrayToProcess(0, inputVector);
  numPts = input->GetNumberOfPoints();

  if ( !inTensors || numPts < 1 )
    {
    vtkErrorMacro(<<"No data!");
    return 1;
    }

  vtkSmartPointer<vtkDoubleArray> result = vtkSmartPointer<vtkDoubleArray>::New();

  result->SetName("result");
  result->SetNumberOfTuples(numPts);


  switch( this->Function ) {
    case 0:  // Magnitude
      vtkDebugMacro(<<"Calculating Magnitude");
      result->SetNumberOfComponents(1);
      for (inPtId=0; inPtId < numPts; inPtId++)
      {
        inTensors->GetTuple(inPtId, tensor);

        // compute magnitude of tensor
        double v, m = 0.0;
        for (i=0; i<9; i++)
        {
          v = tensor[i];
          m += v*v;
        }

        m = sqrt(m);
   
        // store to data array
        result->SetValue(inPtId, m);
      }
      break;

    case 1: // Divergence
      vtkDebugMacro(<<"Calculating Divergence");
      result->SetNumberOfComponents(1);
      for (inPtId=0; inPtId < numPts; inPtId++)
      {
        inTensors->GetTuple(inPtId, tensor);

        // compute divergence as trace of Jacobian
        double trace = 0.0;
        trace += tensor[0];
        trace += tensor[4];
        trace += tensor[8];

        if( this->AbsoluteValue )
        {
          trace = abs(trace);
        }

        // store to data array
        result->SetValue(inPtId, trace);
      }
      break;

    case 2: // Mean Direction 
      vtkDebugMacro(<<"Calculating Mean Direction");
      result->SetNumberOfComponents(3);
      break;
  }


  vtkDebugMacro(<<"Calculated " << numPts <<" VectorCalc values");

  //
  // Update output and release memory
  //

  output->ShallowCopy(input);

  vtkPointData *pd = output->GetPointData();
  pd->AddArray(result);
  
  output->Squeeze();

  return 1;
}

//==============================================================================
void vtkVectorCalc::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
