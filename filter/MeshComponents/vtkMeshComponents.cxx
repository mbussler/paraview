
#include "vtkMeshComponents.h"

#include "vtkCellData.h"
#include "vtkCompositeDataIterator.h"
#include "vtkDataSet.h"
#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkCompositeDataPipeline.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkDistancePolyDataFilter.h"
#include "vtkMath.h"
#include <set>
#include <map>
#include <algorithm>


vtkStandardNewMacro(vtkMeshComponents);

vtkMeshComponents::vtkMeshComponents()
{
}

vtkMeshComponents::~vtkMeshComponents()
{
}

int vtkMeshComponents::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector)
{
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    
    vtkPolyData *input = vtkPolyData::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output = vtkPolyData::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    output->ShallowCopy(input);

    vtkSmartPointer<vtkIntArray> labels;

    labels = vtkSmartPointer<vtkIntArray>::New();
    labels->SetName("labels");
    labels->SetNumberOfComponents(1);
    output->GetCellData()->AddArray(labels);


    output->BuildCells();
    output->BuildLinks();

    // connected components of new cells
    vtkCellArray* cells = output->GetPolys();
    int numCells = cells->GetNumberOfCells();

    std::vector<bool> visited;
    visited.resize(numCells, false);

    labels->SetNumberOfValues(numCells);

    for( vtkIdType cellId=0; cellId<numCells; cellId++ )
    {
        labels->SetValue(cellId, 0);
    }

    int startCell = 0;
    int currentLabel = 1;
    std::vector<vtkIdType> stack;

    while( startCell < numCells)
    {
        stack.push_back(startCell);
        visited[startCell] = true;
        labels->SetValue( startCell, currentLabel);

        while( !stack.empty())
        {
            int currentCell = stack.back();
            stack.pop_back();

            std::set<vtkIdType> cellNeighbors;
            getCellNeighbors( output, currentCell, cellNeighbors );

            std::set<vtkIdType>::iterator cellNeighborsIter;
            for( cellNeighborsIter = cellNeighbors.begin(); cellNeighborsIter != cellNeighbors.end(); cellNeighborsIter++)
            {
                int pos = *cellNeighborsIter;
                if( !visited[pos]) {
                    stack.push_back(pos);
                    labels->SetValue(pos, currentLabel);
                    visited[pos] = true;
                }
            }
        }
        while( startCell < numCells && visited[startCell]) {
            startCell++;
        }
        currentLabel++;
    }

    labels->Modified();

    output->Squeeze();

    return 1;
}

void vtkMeshComponents::getCellNeighbors( vtkPolyData* mesh, vtkIdType cellId, std::set<vtkIdType>& cellNeighbors)
{    
    vtkSmartPointer<vtkIdList> cellPointIds =
        vtkSmartPointer<vtkIdList>::New();
    mesh->GetCellPoints(cellId, cellPointIds);

    for(vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++)
    {
        vtkSmartPointer<vtkIdList> idList =
            vtkSmartPointer<vtkIdList>::New();
        idList->InsertNextId(cellPointIds->GetId(i));

        //add the other edge point
        if(i+1 == cellPointIds->GetNumberOfIds())
        {
            idList->InsertNextId(cellPointIds->GetId(0));
        }
        else
        {
            idList->InsertNextId(cellPointIds->GetId(i+1));
        }

        //get the neighbors of the cell
        vtkSmartPointer<vtkIdList> neighborCellIds =
            vtkSmartPointer<vtkIdList>::New();

        mesh->GetCellNeighbors(cellId, idList, neighborCellIds);

        for(vtkIdType j = 0; j < neighborCellIds->GetNumberOfIds(); j++)
        {
            cellNeighbors.insert( neighborCellIds->GetId(j) );
        }
    }
}
