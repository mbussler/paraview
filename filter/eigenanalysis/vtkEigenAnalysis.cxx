/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkEigenAnalysis.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkEigenAnalysis.h"
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

vtkStandardNewMacro(vtkEigenAnalysis);

//==============================================================================
vtkEigenAnalysis::vtkEigenAnalysis()
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
vtkEigenAnalysis::~vtkEigenAnalysis()
{
  // clean up
}

//==============================================================================
int vtkEigenAnalysis::RequestData(
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

  
  double *m[3], w[3], *v[3];
  double m0[3], m1[3], m2[3];
  double v0[3], v1[3], v2[3];
  double xv[3], yv[3], zv[3];

  // set up working matrices
  m[0] = m0; m[1] = m1; m[2] = m2;
  v[0] = v0; v[1] = v1; v[2] = v2;

  vtkDebugMacro(<<"Calculating Eigen vectors");

  vtkPointData *outPD = output->GetPointData();
  inTensors = this->GetInputArrayToProcess(0, inputVector);
  numPts = input->GetNumberOfPoints();

  if ( !inTensors || numPts < 1 )
    {
    vtkErrorMacro(<<"No data!");
    return 1;
    }

  // only copy scalar data through
  vtkDoubleArray* ev0 = vtkDoubleArray::New();
  vtkDoubleArray* ev1 = vtkDoubleArray::New();
  vtkDoubleArray* ev2 = vtkDoubleArray::New();
  ev0->SetName("ev0");
  ev1->SetName("ev1");
  ev2->SetName("ev2");
  ev0->SetNumberOfComponents(3);
  ev1->SetNumberOfComponents(3);
  ev2->SetNumberOfComponents(3);
  ev0->SetNumberOfTuples(numPts);
  ev1->SetNumberOfTuples(numPts);
  ev2->SetNumberOfTuples(numPts);

  vtkDoubleArray* ew0 = vtkDoubleArray::New();
  vtkDoubleArray* ew1 = vtkDoubleArray::New();
  vtkDoubleArray* ew2 = vtkDoubleArray::New();
  ew0->SetName("ew0");
  ew1->SetName("ew1");
  ew2->SetName("ew2");
  ew0->SetNumberOfComponents(1);
  ew1->SetNumberOfComponents(1);
  ew2->SetNumberOfComponents(1);
  ew0->SetNumberOfTuples(numPts);
  ew1->SetNumberOfTuples(numPts);
  ew2->SetNumberOfTuples(numPts);
  
  //
  // Traverse all Input points, calculate and store eigen vectors
  //

  for (inPtId=0; inPtId < numPts; inPtId++)
    {
    // Translation is postponed

    inTensors->GetTuple(inPtId, tensor);

    // compute orientation vectors and scale factors from tensor
    for (j=0; j<3; j++)
      {
      for (i=0; i<3; i++)
        {
        m[i][j] = tensor[i+3*j];
        }
      }
    vtkMath::Jacobi(m, w, v);

    //copy eigenvectors
    ev0->SetTupleValue(inPtId, v0);
    ev1->SetTupleValue(inPtId, v1);
    ev2->SetTupleValue(inPtId, v2);
    ew0->SetTupleValue(inPtId, &w[0]);
    ew1->SetTupleValue(inPtId, &w[1]);
    ew2->SetTupleValue(inPtId, &w[2]);
    }
  
  vtkDebugMacro(<<"Generated " << numPts <<" Eigenvectors");

  //
  // Update output and release memory
  //

  output->ShallowCopy(input);

  vtkPointData *pd = output->GetPointData();
  pd->AddArray(ev0);
  pd->AddArray(ev1);
  pd->AddArray(ev2);
  pd->AddArray(ew0);
  pd->AddArray(ew1);
  pd->AddArray(ew2);
  
  ev0->Delete();
  ev1->Delete();
  ev2->Delete();

  ew0->Delete();
  ew1->Delete();
  ew2->Delete();
  
  output->Squeeze();

  return 1;
}

//==============================================================================
void vtkEigenAnalysis::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
