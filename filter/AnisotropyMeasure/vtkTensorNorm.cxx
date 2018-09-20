/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTensorNorm.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkTensorNorm.h"
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

#include "linalg.h"

vtkStandardNewMacro(vtkTensorNorm);

//==============================================================================
vtkTensorNorm::vtkTensorNorm()
{
  this->SetNumberOfInputPorts(1);
  
  // by default, process active point tensors
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::TENSORS);

}
//==============================================================================
vtkTensorNorm::~vtkTensorNorm()
{
  // clean up
}

//==============================================================================
int vtkTensorNorm::RequestData(
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

  vtkDataArray *inTensor;
  double tensor[9];
  vtkIdType numPts, inPtId, ptIncr, i;
  int j;

  
  vtkDebugMacro(<<"Calculating TensorNorm");

  vtkPointData *outPD = output->GetPointData();
  inTensor = this->GetInputArrayToProcess(0, inputVector);
  numPts = input->GetNumberOfPoints();

  if ( !inTensor || numPts < 1 )
    {
    vtkErrorMacro(<<"No data!");
    return 1;
    }
  if ( inTensor->GetNumberOfComponents() != 9 )
    {
    vtkErrorMacro("Input array must be a gradient tensor with 9 components.");
    return 0;
    }

  // allocate mem for output data arrays
  vtkSmartPointer<vtkDoubleArray> normArray = vtkSmartPointer<vtkDoubleArray>::New();
  normArray->SetName("Spektralnorm");
  normArray->SetNumberOfComponents(1);
  normArray->SetNumberOfTuples(numPts);

  //
  // Traverse all Input points, calculate and store symmetric and antisymmetric tensors
  //

  for (inPtId=0; inPtId < numPts; inPtId++)
    {

    inTensor->GetTuple(inPtId, tensor);
    mat3 m;

    for (j=0; j<3; j++)
    {
      for (i=0; i<3; i++)
      {
        // compute A' * A
        m[i][j] = tensor[j+3*i] * tensor[i+3*j];
      }
    }

    vec3 lambda;
    mat3eigenvalues(m, lambda);
    double lambdaMax = std::max<double>( lambda[0], std::max<double>( lambda[1], lambda[2]));

    //copy tensors
    normArray->SetTupleValue(inPtId, &lambdaMax );
    }
  
  vtkDebugMacro(<<"Calculated " << numPts <<" tensor norm values");

  //
  // Update output and release memory
  //

  output->ShallowCopy(input);

  vtkPointData *pd = output->GetPointData();
  pd->AddArray(normArray);
 
  output->Squeeze();

  return 1;
}

//==============================================================================
void vtkTensorNorm::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
