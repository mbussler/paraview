/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCellSize.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkCellSize.h"
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
#include "vtkCellData.h"
#include "vtkTriangle.h"

vtkStandardNewMacro(vtkCellSize);

//==============================================================================
vtkCellSize::vtkCellSize()
{
}
//==============================================================================
vtkCellSize::~vtkCellSize()
{
  // clean up
}

//==============================================================================
int vtkCellSize::RequestData(
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

  
  vtkDebugMacro(<<"Calculating CellSize");

  vtkIdType numCells = input->GetNumberOfCells();

  vtkSmartPointer<vtkDoubleArray> cellSizeArray = vtkSmartPointer<vtkDoubleArray>::New();
  cellSizeArray->SetName("CellSize");
  cellSizeArray->SetNumberOfComponents(1);
  cellSizeArray->SetNumberOfTuples(numCells);

  vtkSmartPointer<vtkDoubleArray> areaArray = vtkSmartPointer<vtkDoubleArray>::New();
  areaArray->SetName("Area");
  areaArray->SetNumberOfComponents(1);
  areaArray->SetNumberOfTuples(numCells);

  //
  // Traverse all cells ...
  //

  vtkSmartPointer<vtkIdList> cellPts = vtkSmartPointer<vtkIdList>::New();

  for (vtkIdType cellId=0; cellId < numCells; cellId++)
    {
        double maxlength = 0.0;
        double area = 0;

        vtkCell* cell = input->GetCell(cellId);
        if( cell )
        {
            for( int i=0; i<cell->GetNumberOfEdges(); i++)
            {
                vtkCell* edge = cell->GetEdge(i);
                if( edge && edge->GetNumberOfPoints() == 2 ) {
                    double p1[3], p2[3];
                    input->GetPoint(edge->GetPointId(0), p1);
                    input->GetPoint(edge->GetPointId(1), p2);
                    double length = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
                    maxlength = std::max<double>(length, maxlength);
                }
            }
            if( cell->GetCellType() == VTK_TRIANGLE )
            {
                double p1[3], p2[3], p3[3];
                input->GetPoint(cell->GetPointId(0), p1);
                input->GetPoint(cell->GetPointId(1), p2);
                input->GetPoint(cell->GetPointId(2), p3);
                area = vtkTriangle::TriangleArea(p1, p2, p3);
            }
        }

        cellSizeArray->SetValue(cellId, maxlength);
        areaArray->SetValue(cellId,area);
    }

    
  vtkDebugMacro(<<"Calculated " << numCells <<" CellSize values");

  //
  // Update output and release memory
  //

  output->ShallowCopy(input);
  output->GetCellData()->AddArray(cellSizeArray);
  output->GetCellData()->AddArray(areaArray);
  output->GetCellData()->Modified();

  output->Squeeze();

  return 1;
}

//==============================================================================
void vtkCellSize::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
