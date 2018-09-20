
#include "vtkTimeGeometry.h"

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


vtkStandardNewMacro(vtkTimeGeometry);

vtkTimeGeometry::vtkTimeGeometry()
{
    growMesh = NULL;
    currentTime = 0;
}

vtkTimeGeometry::~vtkTimeGeometry()
{
}

int vtkTimeGeometry::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector)
{
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    
    vtkPolyData *input = vtkPolyData::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output = vtkPolyData::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // get the requested update times
    double upTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    
    int     numTimes = 0;
    double *inTimes  = NULL;
    bool noTimesteps = false;
    if( inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS())) {
        inTimes = inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        numTimes = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        noTimesteps = numTimes <= 1;
    } else {
        noTimesteps = true;
    }

    if( noTimesteps ) {
        vtkErrorMacro("No timesteps to grow.");
        return 0;
    }

    // init/reset current mesh and ref point vector
    int resetTimestep = std::max<int>(std::min<int>(this->ResetOnTimestep, numTimes-1), 0);
    if( growMesh.GetPointer() == NULL || upTime <= inTimes[resetTimestep])
    {
        growMesh = vtkSmartPointer<vtkPolyData>::New();
        growMesh->SetPoints(vtkPoints::New());
        growMesh->Allocate(); // init cell data arrays
        currentTime = upTime;
    }

    if( upTime <= currentTime) {
        output->DeepCopy(growMesh);
        return 1;
    }
    else
    {
        currentTime = upTime;
    }

    vtkPointData *inPD = input->GetPointData();
    vtkCellData *inCD = input->GetCellData();
    int numArrays = inPD->GetNumberOfArrays();
    int numCellArrays = inCD->GetNumberOfArrays();

    vtkPointData *outPD = growMesh->GetPointData();
    vtkCellData  *outCD = growMesh->GetCellData();

    // copy point data arrays
    for( int i=0; i<numArrays; i++ )
    {
        if( !outPD->HasArray(inPD->GetArray(i)->GetName())) 
        {
            int type = inPD->GetArray(i)->GetDataType();
            vtkDataArray *arr = vtkDataArray::CreateDataArray(type);
            arr->SetNumberOfComponents(inPD->GetArray(i)->GetNumberOfComponents());
            arr->SetName(inPD->GetArray(i)->GetName());
            outPD->AddArray(arr);
            arr->Delete();
        }
    }

    // copy cell data arrays
    for( int i=0; i<numCellArrays; i++ )
    {
        if( !outCD->HasArray(inCD->GetArray(i)->GetName())) 
        {
            int type = inCD->GetArray(i)->GetDataType();
            vtkDataArray *arr = vtkDataArray::CreateDataArray(type);
            arr->SetNumberOfComponents(inCD->GetArray(i)->GetNumberOfComponents());
            arr->SetName(inCD->GetArray(i)->GetName());
            outCD->AddArray(arr);
            arr->Delete();
        }
    }

    std::vector<std::pair<vtkDataArray*, vtkDataArray*> > arraysPD;
    for( int i=0; i<numArrays; i++ )
    {
        std::pair<vtkDataArray*, vtkDataArray*> arr;
        arr.first = inPD->GetArray(i);
        arr.second = outPD->GetArray(arr.first->GetName());
        arraysPD.push_back(arr);
    }

    std::vector<std::pair<vtkDataArray*, vtkDataArray*> > arraysCD;
    for( int i=0; i<numCellArrays; i++ )
    {
        std::pair<vtkDataArray*, vtkDataArray*> arr;
        arr.first = inCD->GetArray(i);
        arr.second = outCD->GetArray(arr.first->GetName());
        arraysCD.push_back(arr);
    }

    vtkIdType numPts = input->GetNumberOfPoints();
    std::vector<vtkIdType> newPtIds;
    newPtIds.resize(numPts, -1);

    vtkSmartPointer<vtkPoints> points = growMesh->GetPoints();

    for( vtkIdType ptId=0; ptId<numPts; ptId++)
    {
        double pos[3];
        input->GetPoint(ptId, pos);

        // insert point into output
        vtkIdType newId = points->InsertNextPoint(pos);
        newPtIds[ptId] = newId;

        // copy point data
        for( int i=0; i<numArrays; i++ ){
            double data[9];
            arraysPD[i].first->GetTuple(ptId, data);
            arraysPD[i].second->InsertNextTuple(data);
        }
    }
    
    vtkSmartPointer<vtkIdList> newCellPtIds;

    // iterate over all cells, insert cells that do not belong to surface
    for( vtkIdType cellId=0; cellId<input->GetNumberOfCells(); cellId++) 
    {
        vtkCell* cell = input->GetCell(cellId);

        if( cell ){
            vtkIdList* cellPtIds = cell->GetPointIds();
            int numPtIds = cellPtIds->GetNumberOfIds();
            newCellPtIds = vtkSmartPointer<vtkIdList>::New();

            bool active = true;
            for( int i=0; i<numPtIds; i++) {
                int newPtId = newPtIds[cellPtIds->GetId(i)];
                newCellPtIds->InsertNextId( newPtId );
                active &= (newPtId > 0);
            }

            if( active ) {
                growMesh->InsertNextCell(cell->GetCellType(), newCellPtIds);

                // copy cell data
                for( int i=0; i<numCellArrays; i++ ){
                    double data[9];
                    arraysCD[i].first->GetTuple(cellId, data);
                    arraysCD[i].second->InsertNextTuple(data);
                }
            }

        }
    }

    output->ShallowCopy(growMesh);
    output->Squeeze();

    return 1;
}
