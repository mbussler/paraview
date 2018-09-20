/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkRidgeAngle.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkRidgeAngle.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
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
#include "vtkCellData.h"
#include "vtkTriangle.h"
#include "vtkPolyDataNormals.h"

vtkStandardNewMacro(vtkRidgeAngle);

//==============================================================================
vtkRidgeAngle::vtkRidgeAngle()
{
}
//==============================================================================
vtkRidgeAngle::~vtkRidgeAngle()
{
  // clean up
}

//==============================================================================
int vtkRidgeAngle::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  
  vtkDebugMacro(<<"Calculating RidgeAngle");

  vtkIdType numCells = input->GetNumberOfCells();

  vtkSmartPointer<vtkDoubleArray> RidgeAngleArray = vtkSmartPointer<vtkDoubleArray>::New();
  RidgeAngleArray->SetName("RidgeAngle");
  RidgeAngleArray->SetNumberOfComponents(1);
  RidgeAngleArray->SetNumberOfTuples(numCells);


  //
  // Traverse all cells ...
  //

  vtkSmartPointer<vtkPolyDataNormals> calcNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
  calcNormals->SetInputData(input);
  calcNormals->ComputeCellNormalsOn();
  calcNormals->ComputePointNormalsOff();
  calcNormals->SetSplitting(1);
  calcNormals->SetConsistency(1);
  calcNormals->SetNonManifoldTraversal(0);
  calcNormals->Update();

  vtkDataArray* cellNormals = calcNormals->GetOutput()->GetCellData()->GetArray("Normals");

  vtkSmartPointer<vtkIdList> cellPts = vtkSmartPointer<vtkIdList>::New();

  for (vtkIdType cellId=0; cellId < numCells; cellId++)
    {
        double angle = 0.0;
        double normal[3];

        cellNormals->GetTuple(cellId, normal);


        //vtkCell* cell = input->GetCell(cellId);
        //if( cell && cell->GetCellType() == VTK_TRIANGLE)
        //{
        //    double normal[3];
        //    double p1[3], p2[3], p3[3];

        //    input->GetPoint( cell->GetPointId(0), p1);
        //    input->GetPoint( cell->GetPointId(1), p2);
        //    input->GetPoint( cell->GetPointId(2), p3);

        //    vtkTriangle::ComputeNormal(p1, p2, p3, normal);

        double v[3] = {1,0,0};
        double dot = vtkMath::Dot(v,normal);

        angle = 180.0*acos(dot)/3.141592654;
        if( angle >= 0)
            angle -= 90.0;
        else if( angle < 0)
            angle += 90.0;

        //    //if( angle > 180.0)
        //    //    angle -= 180.0;
        //    //else if( angle < -180.0)
        //    //    angle += 180.0;
        //}

        angle = abs(angle);

        RidgeAngleArray->SetValue(cellId, angle);
    }

    
  vtkDebugMacro(<<"Calculated " << numCells <<" RidgeAngle values");

  //
  // Update output and release memory
  //

  output->ShallowCopy(input);
  output->GetCellData()->AddArray(RidgeAngleArray);
  output->GetCellData()->Modified();

  output->Squeeze();

  return 1;
}

//==============================================================================
void vtkRidgeAngle::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
