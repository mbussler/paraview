/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCustomDataSet.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkCustomDataSet.h"

#include "vtkDataObjectToDataSetFilter.h"

#include "vtkUnstructuredGrid.h"
#include "vtkCellArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkObjectFactory.h"
#include "vtkQuadraticTetra.h"
#include "vtkCustomQuadraticTetra.h"

vtkStandardNewMacro(vtkCustomDataSet);

//==============================================================================
vtkCustomDataSet::vtkCustomDataSet()
{
  this->SetNumberOfInputPorts(1);
  vtkDataSet *output = vtkUnstructuredGrid::New();
  this->GetExecutive()->SetOutputData(0, output);
  
  output->ReleaseData();
  output->Delete();
}
//==============================================================================
vtkCustomDataSet::~vtkCustomDataSet()
{
  // clean up
}

//----------------------------------------------------------------------------
vtkDataObject *vtkCustomDataSet::GetInput()
{
  if (this->GetNumberOfInputConnections(0) < 1)
    {
    return NULL;
    }

  return this->GetExecutive()->GetInputData(0, 0);
}

//----------------------------------------------------------------------------
int vtkCustomDataSet::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkDataObject *input = inInfo->Get(vtkDataObject::DATA_OBJECT());

  return 1;
}
//==============================================================================
int vtkCustomDataSet::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkUnstructuredGrid *input = vtkUnstructuredGrid::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if( input ) {
    vtkDebugMacro(<<"Converting input to CustomDataSet");

    output->ShallowCopy(input);

    vtkCellArray* inCells = input->GetCells();
    for( vtkIdType cellId = 0; cellId < inCells->GetNumberOfCells(); cellId++)
    {
    }
    
    vtkDebugMacro(<<"Finished!");
    
    return 1;
  }
  else 
  {
    vtkErrorMacro(<<"Unsupported dataset type!");

    return 0;  
  }
}
//----------------------------------------------------------------------------
int vtkCustomDataSet::FillInputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
  return 1;
}

//----------------------------------------------------------------------------
int vtkDataObjectToDataSetFilter::RequestDataObject(
  vtkInformation *,
  vtkInformationVector **,
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkDataSet *output = vtkDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  if (!output || (output->GetDataObjectType() != VTK_UNSTRUCTURED_GRID))
  {
    output = vtkUnstructuredGrid::New();
    if (output)
      {
      outInfo->Set(vtkDataObject::DATA_OBJECT(), output);
      output->Delete();
      }
    }
  return 1;
}

//==============================================================================
void vtkCustomDataSet::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
