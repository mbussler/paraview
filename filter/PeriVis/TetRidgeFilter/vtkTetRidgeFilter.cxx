#include "vtkTetRidgeFilter.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkGenericCell.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkLine.h"
#include "vtkMergePoints.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkTriangle.h"
#include "vtkIncrementalPointLocator.h"
#include "vtkTransform.h"

#include <math.h>
#include <algorithm>
#include "vtkUnstructuredGrid.h"

#include <map>
#include <set>

/* performance measure */
#include "timer.h"
#include <QElapsedTimer>

vtkStandardNewMacro(vtkTetRidgeFilter);

template <typename T>
inline void swap( T& a, T& b) {
    T swp=a; a=b; b=swp;
};

//----------------------------------------------------------------------------
vtkTetRidgeFilter::vtkTetRidgeFilter()
{
    // by default process active point scalars
    this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
        vtkDataSetAttributes::SCALARS);

}

//----------------------------------------------------------------------------
vtkTetRidgeFilter::~vtkTetRidgeFilter()
{
}

//----------------------------------------------------------------------------
//
// Clip through data generating surface.
//
int vtkTetRidgeFilter::RequestData(
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

    vtkIdType cellId, i, updateTime;
    vtkSmartPointer<vtkCellArray> newTris;
    vtkSmartPointer<vtkPoints> newPoints;

    vtkIdType estimatedSize, numCells=input->GetNumberOfCells();
    vtkIdType numPts=input->GetNumberOfPoints();
    vtkPoints *inPts=input->GetPoints();
    int numberOfPoints;
    vtkPointData *inPD=input->GetPointData(), *outPD = output->GetPointData();
    vtkCellData *inCD=input->GetCellData(), *outCD = output->GetCellData();
    vtkCellData *outClippedCD = NULL;

    vtkDataArray *scalarArray = this->GetInputArrayToProcess(0, inputVector);
    //if( !scalarArray )
    //{
    //    vtkErrorMacro("No input array selected");
    //    return 0;
    //}

    vtkDebugMacro(<< "Calculating iso surface");

    // Initialize self; create output objects
    //
    if ( !input || numPts < 1 || inPts == NULL )
    {
        vtkDebugMacro(<<"No data.");
        return 1;
    }

    // new point data
    newPoints = vtkSmartPointer<vtkPoints>::New();
    newPoints->DeepCopy(input->GetPoints());
    output->SetPoints(newPoints);
    output->GetPointData()->DeepCopy(input->GetPointData());

    // cells
    newTris = vtkSmartPointer<vtkCellArray>::New();

    //output->ShallowCopy(input);

    QElapsedTimer timer;
    timer.start();

    vtkCellArray* polys = input->GetPolys();

    if( polys && polys->GetNumberOfCells() > 0)
    {
        //// copy cell data arrays
        //int numCellArrays = inCD->GetNumberOfArrays();
        //for( int i=0; i<numCellArrays; i++ ){
        //    int type = inCD->GetArray(i)->GetDataType();
        //    vtkDataArray *arr = vtkDataArray::CreateDataArray(type);
        //    arr->SetNumberOfComponents(inCD->GetArray(i)->GetNumberOfComponents());
        //    arr->SetName(inCD->GetArray(i)->GetName());
        //    outCD->AddArray(arr);
        //    arr->Delete();
        //}

        vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();
        polys->InitTraversal();
        //vtkIdType cellId = 0;
        while(polys->GetNextCell(pts))
        {
            bool valid = true;
            if(scalarArray)
            {
                for( int i=0; i<pts->GetNumberOfIds(); i++)
                {
                    double val = scalarArray->GetTuple1(pts->GetId(i));
                    valid &= (val >= this->MinDataValue && val <= this->MaxDataValue);
                }
            }
            if( valid) 
            {
                newTris->InsertNextCell(pts);
                // copy cell data
                double data[9];
                //for( int i=0; i<numCellArrays; i++ ){
                //    inCD->GetArray(i)->GetTuple(cellId, data);
                //    outCD->GetArray(i)->InsertNextTuple(data);
                //}
            }
            //cellId++;
        }
    }
    output->SetPolys(newTris);

    //vtkSmartPointer<vtkIntArray> labels = vtkSmartPointer<vtkIntArray>::New();
    //labels->SetName("labels");
    //labels->SetNumberOfComponents(1);

    std::vector<std::vector<vtkIdType> > components;

    if( this->RegionThreshold > 0) 
    {
        output->BuildCells();
        output->BuildLinks();

        // connected components of new cells
        int numCells = output->GetNumberOfPolys();
        std::vector<bool> visited;
        visited.resize(numCells, false);
        //labels->SetNumberOfTuples(numCells);

        int startCell = 0;
        int currentLabel = 1;
        std::vector<vtkIdType> stack;

        while( startCell < numCells)
        {
            stack.push_back(startCell);
            visited[startCell] = true;
            //labels->SetValue(startCell, currentLabel);

            std::vector<vtkIdType> component;
            component.push_back( startCell );

            while( !stack.empty())
            {
                int currentCell = stack.back();
                stack.pop_back();

                std::vector<vtkIdType> cellNeighbors;
                getCellNeighbors( output, currentCell, cellNeighbors );

                std::vector<vtkIdType>::iterator cellNeighborsIter;
                for( cellNeighborsIter = cellNeighbors.begin(); cellNeighborsIter != cellNeighbors.end(); cellNeighborsIter++)
                {
                    int neighbor = *cellNeighborsIter;
                    if( !visited[neighbor]) {
                        stack.push_back(neighbor);
                        component.push_back(neighbor);
                        visited[neighbor] = true;
                        //labels->SetValue(neighbor, currentLabel);
                    }
                }
            }
            while( startCell < numCells && visited[startCell]) {
                startCell++;
            }
            currentLabel++;
            components.push_back(component);
        }

        newTris = vtkSmartPointer<vtkCellArray>::New();

        for( std::vector<std::vector<vtkIdType> >::iterator componentsIter = components.begin();
            componentsIter != components.end(); componentsIter++)
        {
            int size = componentsIter->size();
            if( size >= this->RegionThreshold )
            {
                for( std::vector<vtkIdType>::iterator componentIter = componentsIter->begin();
                    componentIter != componentsIter->end(); componentIter++)
                {
                    vtkSmartPointer<vtkIdList> cellPointIds =
                        vtkSmartPointer<vtkIdList>::New();
                    output->GetCellPoints(*componentIter, cellPointIds);

                    newTris->InsertNextCell(cellPointIds);
                }
            }
        }
        output->SetPolys(newTris);
    }


    write_timer("ridgeFilter", "filter", timer.elapsed());

    output->BuildCells();
    //output->GetCellData()->AddArray(labels);
    
    output->Squeeze();

    return 1;
}

void vtkTetRidgeFilter::getCellNeighbors( vtkPolyData* mesh, vtkIdType cellId, std::vector<vtkIdType>& cellNeighbors)
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
            cellNeighbors.push_back( neighborCellIds->GetId(j) );
        }
    }

}
