/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGrowthVelocity.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkGrowthVelocity.h"
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
#include "vtkTable.h"

vtkStandardNewMacro(vtkGrowthVelocity);

//==============================================================================
vtkGrowthVelocity::vtkGrowthVelocity()
{
    this->SetNumberOfOutputPorts(2);
}

int vtkGrowthVelocity::FillOutputPortInformation(int port, vtkInformation* info)
{
    if ( port == 0 ) {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet" );
        return 1;
    }
    if ( port == 1 ) {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable" );
        return 1;
    }
    return 0;
}

//==============================================================================
vtkGrowthVelocity::~vtkGrowthVelocity()
{
  // clean up
}

//==============================================================================
int vtkGrowthVelocity::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkInformation *outInfo1 = outputVector->GetInformationObject(1);

  // get the input and output
  vtkDataSet *input = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataSet *output = vtkDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkTable *outTable = vtkTable::
      SafeDownCast(outInfo1->Get(vtkDataObject::DATA_OBJECT()));

  
  vtkDebugMacro(<<"Calculating GrowthVelocity");

  vtkIdType numCells = input->GetNumberOfCells();

  vtkDoubleArray* cellGrowth = vtkDoubleArray::SafeDownCast(
      input->GetCellData()->GetArray("CellGrowth"));

  if( !cellGrowth ) {
      return 0;
  }

  vtkSmartPointer<vtkDoubleArray> velocityArray = vtkSmartPointer<vtkDoubleArray>::New();
  velocityArray->SetName("Velocity");
  velocityArray->SetNumberOfComponents(1);
  velocityArray->SetNumberOfTuples(numCells);

  int numRanges = this->Buckets;
  std::vector<double> ranges;
  for( int i=0; i<numRanges; i++){
      ranges.push_back(this->Range*i);
  }

  std::vector<double> areas;
  areas.resize(numRanges, 0.0);

  //
  // Traverse all cells ...
  //

  double v_max = 0.0;

  for (vtkIdType cellId=0; cellId < numCells; cellId++)
    {
        double area = 0;

        vtkCell* cell = input->GetCell(cellId);
        if( cell && cell->GetCellType() == VTK_TRIANGLE )
        {
            double p1[3], p2[3], p3[3];
            input->GetPoint(cell->GetPointId(0), p1);
            input->GetPoint(cell->GetPointId(1), p2);
            input->GetPoint(cell->GetPointId(2), p3);
            area = vtkTriangle::TriangleArea(p1, p2, p3);
        }

        double growthRate = cellGrowth->GetValue(cellId);
        double velocity = growthRate / this->Time;
        v_max = std::max<double>(v_max, velocity);

        int bucket = 0;
        while( velocity > ranges[bucket] && bucket < numRanges-1) {
            bucket++;
        }

        areas[bucket] += area;

        velocityArray->SetValue(cellId, velocity);
    }

  vtkSmartPointer<vtkDoubleArray> rangesArray = 
      vtkSmartPointer<vtkDoubleArray>::New();
  rangesArray->SetName("Velocity Ranges");
  vtkSmartPointer<vtkDoubleArray> areasArray = 
      vtkSmartPointer<vtkDoubleArray>::New();
  areasArray->SetName("Area");
  
  outTable->AddColumn(rangesArray);
  outTable->AddColumn(areasArray);
  outTable->SetNumberOfRows(numRanges);

  // output..
  std::cout << "\n************************\n";

  for( int i=0; i<numRanges; i++)
  {
      std::cout << ranges[i] << "\t" << areas[i] << "\n";
      outTable->SetValue(i, 0, ranges[i]);
      outTable->SetValue(i, 1, areas[i]);
  }
  std::cout << "\nV_max: " << v_max << "\n";
  std::cout << "\n************************\n";

  
  //
  // Update output and release memory
  //

  output->ShallowCopy(input);
  output->GetCellData()->AddArray(velocityArray);
  output->Squeeze();

  return 1;
}

//==============================================================================
void vtkGrowthVelocity::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
