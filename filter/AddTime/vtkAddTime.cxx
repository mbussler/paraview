/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkAddTime.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkAddTime.h"
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

vtkStandardNewMacro(vtkAddTime);

//==============================================================================
vtkAddTime::vtkAddTime()
{
  this->SetNumberOfInputPorts(1);
  
  // by default, process active point tensors
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                                                  vtkDataSetAttributes::TENSORS);


}
//==============================================================================
vtkAddTime::~vtkAddTime()
{
  // clean up
}

//==============================================================================
int vtkAddTime::RequestData(
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

  
  vtkDebugMacro(<<"Calculating AddTime");

  vtkIdType numPts = input->GetNumberOfPoints();

  vtkSmartPointer<vtkDoubleArray> timeArray = vtkSmartPointer<vtkDoubleArray>::New();
  timeArray->SetName("AddTime");
  timeArray->SetNumberOfComponents(1);
  timeArray->SetNumberOfTuples(numPts);

  double upTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

  //
  // Traverse all Input points, calculate and store AddTime
  //

  for (vtkIdType inPtId=0; inPtId < numPts; inPtId++)
    {

    // store to data array
    timeArray->SetValue(inPtId, upTime);
    
    }
  
  vtkDebugMacro(<<"Calculated " << numPts <<" AddTime values");

  //
  // Update output and release memory
  //

  output->ShallowCopy(input);

  vtkPointData *pd = output->GetPointData();
  pd->AddArray(timeArray);

  pd->Modified();
  
  output->Squeeze();

  return 1;
}

//==============================================================================
void vtkAddTime::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
