/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRateOfStrain.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkRateOfStrain.h"
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

vtkStandardNewMacro(vtkRateOfStrain);

//==============================================================================
vtkRateOfStrain::vtkRateOfStrain()
{
  this->SetNumberOfInputPorts(1);
  
  // by default, process active point tensors
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::TENSORS);

  // by default, process active point scalars
  this->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::SCALARS);

}
//==============================================================================
vtkRateOfStrain::~vtkRateOfStrain()
{
  // clean up
}

//==============================================================================
int vtkRateOfStrain::RequestData(
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

  vtkDataArray *inJacobi;
  double jacobi[9];
  vtkIdType numPts, inPtId, ptIncr, i;
  int j;

  
  // allocate mem for temporal matrices
  double em[9], rm[9];

  vtkDebugMacro(<<"Calculating RateOfStrain");

  vtkPointData *outPD = output->GetPointData();
  inJacobi = this->GetInputArrayToProcess(0, inputVector);
  numPts = input->GetNumberOfPoints();

  if ( !inJacobi || numPts < 1 )
    {
    vtkErrorMacro(<<"No data!");
    return 1;
    }
  if ( inJacobi->GetNumberOfComponents() != 9 )
    {
    vtkErrorMacro("Input array must be a gradient tensor with 9 components.");
    return 0;
    }

  // allocate mem for output data arrays
  vtkSmartPointer<vtkDoubleArray> E_tens = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> R_tens = vtkSmartPointer<vtkDoubleArray>::New();
  E_tens->SetName("symmetric strain rate");
  R_tens->SetName("antisymmetric strain rate");
  E_tens->SetNumberOfComponents(9);
  R_tens->SetNumberOfComponents(9);
  E_tens->SetNumberOfTuples(numPts);
  R_tens->SetNumberOfTuples(numPts);

  //
  // Traverse all Input points, calculate and store symmetric and antisymmetric tensors
  //

  for (inPtId=0; inPtId < numPts; inPtId++)
    {

    inJacobi->GetTuple(inPtId, jacobi);

    for (j=0; j<3; j++)
    {
      for (i=0; i<3; i++)
      {
        // compute symmetric part of jacobi
        em[i+3*j] = ( jacobi[i+3*j] + jacobi[j+3*i] ) / 2.0;
        // compute antisymmetric part of jacobi
        rm[i+3*j] = ( jacobi[i+3*j] - jacobi[j+3*i] ) / 2.0;
      }
    }

    //copy tensors
    E_tens->SetTupleValue(inPtId, em);
    R_tens->SetTupleValue(inPtId, rm);
    }
  
  vtkDebugMacro(<<"Generated " << numPts <<" Eigenvectors");

  //
  // Update output and release memory
  //

  output->ShallowCopy(input);

  vtkPointData *pd = output->GetPointData();
  pd->AddArray(E_tens);
  pd->AddArray(R_tens);
 
  output->Squeeze();

  return 1;
}

//==============================================================================
void vtkRateOfStrain::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
