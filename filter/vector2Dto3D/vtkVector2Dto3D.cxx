/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkVector2Dto3D.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkVector2Dto3D.h"
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

vtkStandardNewMacro(vtkVector2Dto3D);

//==============================================================================
vtkVector2Dto3D::vtkVector2Dto3D()
{
  this->SetNumberOfInputPorts(1);
  
  // by default, process active point tensors
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::VECTORS);

}
//==============================================================================
vtkVector2Dto3D::~vtkVector2Dto3D()
{
  // clean up
}

//==============================================================================
int vtkVector2Dto3D::RequestData(
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

   vtkIdType numPts, inPtId;
 
   
   vtkDataArray *inVectors = this->GetInputArrayToProcess(0, inputVector);
   numPts = input->GetNumberOfPoints();
 
   if ( !inVectors|| numPts < 1 )
     {
     vtkErrorMacro(<<"No data!");
     return 1;
     }
     
   if ( inVectors->GetNumberOfComponents() != 2 )
     {
     vtkErrorMacro("Input array must be a vector with 2 components.");
     return 0;
     }
     
   vtkSmartPointer<vtkDoubleArray> result = vtkSmartPointer<vtkDoubleArray>::New();
 
   result->SetName("vector3D");
   result->SetNumberOfComponents(3);
   result->SetNumberOfTuples(numPts);
   
   double vecIn[2], vecOut[3];
   for (inPtId=0; inPtId < numPts; inPtId++)
     {
       inVectors->GetTuple( inPtId, vecIn);
       vecOut[0] = vecIn[0];
       vecOut[1] = vecIn[1];
       vecOut[2] = 0.0;

       // store to data array
       result->SetTupleValue( inPtId, vecOut);
     }

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
void vtkVector2Dto3D::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
