/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCellFlux.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkCellFlux.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkCell.h"
#include "vtkCellData.h"
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
#include "vtkTriangle.h"

#include <stdlib.h>

vtkStandardNewMacro(vtkCellFlux);

//==============================================================================
vtkCellFlux::vtkCellFlux()
{
  this->SetNumberOfInputPorts(1);
  
  // by default, process active point tensors
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::VECTORS);

}
//==============================================================================
vtkCellFlux::~vtkCellFlux()
{
  // clean up
}

//==============================================================================
int vtkCellFlux::RequestData(
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

  vtkIdType cellId, ptId, faceId, facePtId;
  vtkIdType numCells, numPts, numFaces;
  vtkPointData *inPD=input->GetPointData();
  vtkCellData *outCD=output->GetCellData();
  int maxCellSize=input->GetMaxCellSize();
  vtkIdList *cellPts, *facePts;
  vtkCell *cell, *face;
  vtkDataArray *inVelocities;
  double vel[3], faceVel[3], fcoords[3], x[3];
  double *weights;
  int subId;

  double pointCoordsA[3], pointCoordsB[3], AC[3];
  double base, height, area, cellFlux;
  
  if ( (numCells=input->GetNumberOfCells()) < 1 )
    {
    vtkDebugMacro(<<"No input cell data!");
    return 1;
    }

  inVelocities = this->GetInputArrayToProcess(0, inputVector);
  numPts = input->GetNumberOfPoints();

  if ( !inVelocities || numPts < 1 )
  {
    vtkErrorMacro(<<"No data!");
    return 1;
  }

  vtkSmartPointer<vtkDoubleArray> flux = vtkSmartPointer<vtkDoubleArray>::New();
  flux->SetName("flux");
  flux->SetNumberOfTuples(numCells);
  flux->SetNumberOfComponents(1);
  
  weights = new double [input->GetMaxCellSize()];

  vtkDebugMacro(<<"Calculating Cell Flux");

  // 1) iterate over all cells
  // 2) foreach cell: 
  // 2.1) iterate over cell faces and sum up interpolated velocity at face center
  // 3) calc flux as sum of face fluxes
  
  
  int abort=0;
  vtkIdType progressInterval = numCells/10 + 1;
  int hasEmptyCells = 0;
  for (cellId=0; cellId < numCells && !abort; cellId++)
  {
    if ( ! (cellId % progressInterval) )
    {
      vtkDebugMacro(<<"Processing #" << cellId);
      this->UpdateProgress (0.5*cellId/numCells);
      abort = this->GetAbortExecute();
    }
    
    cell = input->GetCell(cellId);

    numFaces = cell->GetNumberOfFaces();
    if (cell->GetCellType() != VTK_EMPTY_CELL && numFaces > 0)
    {  

      cell = input->GetCell(cellId);
      double cellCenter[3];
      cell->GetParametricCenter(cellCenter);
      cellFlux = 0.0;
      
      for( faceId=0; faceId < numFaces; faceId++) 
      {
        // interpolate face velocity and store in faceVel
        face = cell->GetFace( faceId);
        subId = face->GetParametricCenter(fcoords);
        face->EvaluateLocation(subId, fcoords, x, weights); 
        
        //calculate face normal
        double faceNormal[3];
        vtkMath::Subtract( fcoords, cellCenter, faceNormal); // faceNormal = fcoords - cellCenter
        vtkMath::Normalize(faceNormal);
        
        faceVel[0] = faceVel[1] = faceVel[2] = 0.0;
        
        // calculate face area
        double vfaceArea[3] = {0.0, 0.0, 0.0};
        int numPoints = face->GetNumberOfPoints();
        
        for( facePtId = 0; facePtId < numPoints; facePtId++)
        {
          // interpolate velocity at face center
          ptId = face->GetPointId(facePtId);
          inVelocities->GetTuple(ptId, vel);
          faceVel[0] += weights[facePtId] * vel[0];
          faceVel[1] += weights[facePtId] * vel[1];
          faceVel[2] += weights[facePtId] * vel[2];
          
          // split face into triangles and sum up area
          cell->GetPoints()->GetPoint(ptId, pointCoordsA);
          cell->GetPoints()->GetPoint((ptId+1)%numPoints, pointCoordsB );

          // calculate triangle normal and area
          vtkTriangle::ComputeNormal(pointCoordsA, pointCoordsB, fcoords, faceNormal);
          area = vtkTriangle::TriangleArea(pointCoordsA, pointCoordsB, fcoords);

          // dA_triangle
          vtkMath::MultiplyScalar( faceNormal, area);

          // dA_total += dA_triangle
          vtkMath::Add(faceNormal, vfaceArea, vfaceArea);
        }
        
        // calculate face flux as face normal * face area
        cellFlux += vtkMath::Dot(faceVel, faceNormal);        
      }
      
      flux->SetTuple1(cellId, cellFlux);
    }
    else
    {
      hasEmptyCells = 1;
    }
  }  

  vtkDebugMacro(<<"Calculated " << numPts <<" CellFlux values");

  //
  // Update output and release memory
  //

  output->ShallowCopy(input);

  outCD->AddArray(flux);
  
  output->Squeeze();

  return 1;
}

//==============================================================================
void vtkCellFlux::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
