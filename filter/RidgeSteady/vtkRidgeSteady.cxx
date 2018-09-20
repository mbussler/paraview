
#include "vtkRidgeSteady.h"

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
#include "vtkHausdorffDistancePointSetFilter.h"
#include <set>
#include <map>


vtkStandardNewMacro(vtkRidgeSteady);

vtkRidgeSteady::vtkRidgeSteady()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);

    oldMesh = NULL;
}

vtkRidgeSteady::~vtkRidgeSteady()
{
}

int vtkRidgeSteady::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0) {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    }
    return 1;
}

int vtkRidgeSteady::FillOutputPortInformation(int port, vtkInformation* info)
{
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
}

int vtkRidgeSteady::RequestData(vtkInformation *request,
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

    int numPts = input->GetNumberOfPoints();
    int numCells = input->GetNumberOfPolys();

    vtkPointData *inPD = input->GetPointData();

    // init/reset current mesh and ref point vector
    int resetTimestep = std::max<int>(std::min<int>(this->ResetOnTimestep, numTimes-1), 0);
    if( oldMesh == NULL ||  upTime <= inTimes[resetTimestep])
    {
        output->ShallowCopy(input);
        output->Squeeze();

        vtkSmartPointer<vtkDoubleArray> timeValues = vtkSmartPointer<vtkDoubleArray>::New();
        timeValues->SetNumberOfComponents(1);
        timeValues->SetNumberOfTuples(output->GetNumberOfPoints());
        timeValues->SetName("Time");
        for( int i=0; i<output->GetNumberOfPoints(); i++)
            timeValues->SetValue(i, upTime);
        output->GetPointData()->AddArray(timeValues);

        oldMesh = vtkSmartPointer<vtkPolyData>::New();
        oldMesh->ShallowCopy(output);
        
        return 1;
    }

    // Initialize output
    output->ShallowCopy(input);
    vtkPointData *outPD = output->GetPointData();

    // store current time in array
    vtkSmartPointer<vtkDoubleArray> timeValues = vtkSmartPointer<vtkDoubleArray>::New();
    timeValues->SetNumberOfComponents(1);
    timeValues->SetNumberOfTuples(output->GetNumberOfPoints());
    timeValues->SetName("Time");
    for( int i=0; i<output->GetNumberOfPoints(); i++)
        timeValues->SetValue(i, upTime);
    outPD->AddArray(timeValues);

    // copy extra points from old to current
    vtkSmartPointer<vtkPoints> points = output->GetPoints();

    // calculate per point distance of old mesh to input mesh
    vtkSmartPointer<vtkHausdorffDistancePointSetFilter> distanceFilter =
        vtkSmartPointer<vtkHausdorffDistancePointSetFilter>::New();
    distanceFilter->SetTargetDistanceMethod(
        vtkHausdorffDistancePointSetFilter::POINT_TO_POINT /*POINT_TO_CELL*/);

    distanceFilter->SetInputData(0, oldMesh );
    distanceFilter->SetInputData(1, input );
    distanceFilter->Update();

    vtkPolyData* res = distanceFilter->GetOutput(0); // dist old -> input
    vtkDataArray* distArray = res->GetPointData()->GetArray("Distance");
    vtkIdTypeArray* nearestIdArray = vtkIdTypeArray::SafeDownCast(
        res->GetPointData()->GetArray("NearestPointId"));
    

    int numArrays = oldMesh->GetPointData()->GetNumberOfArrays();
    std::vector<int> newPtIds;
    newPtIds.resize(oldMesh->GetNumberOfPoints(), -1);

    if( distArray )
    {
        for( vtkIdType ptId=0; ptId<oldMesh->GetNumberOfPoints(); ptId++)
        {
            double dist = distArray->GetTuple1(ptId);
            
            // check if point is apart from input mesh
            if( dist >= this->MinDistance) 
            {
                // insert point into output
                double pos[3];
                res->GetPoint(ptId, pos);
                vtkIdType newId = points->InsertNextPoint(pos);
                newPtIds[ptId] = newId;

                // copy point data
                for( int i=0; i<numArrays; i++ ){
                    double data[9];
                    oldMesh->GetPointData()->GetArray(i)->GetTuple(ptId, data);
                    outPD->GetArray(i)->InsertNextTuple(data);
                }
                //timeValues->InsertNextValue(upTime);
            }
        }
    }

    // copy cells
    vtkSmartPointer<vtkCellArray> cells = output->GetPolys();
    vtkCellArray* tris_ref = res->GetPolys();
    tris_ref->InitTraversal();
    vtkSmartPointer<vtkIdList> cellPts = vtkSmartPointer<vtkIdList>::New();
    vtkIdType cellId = 0;

    while( tris_ref->GetNextCell(cellPts))
    {
        // check if all cell points are active
        vtkSmartPointer<vtkIdList> newCellPtIds = vtkSmartPointer<vtkIdList>::New();

        bool active = true;
        if( this->MergeRidges )
            active = false;

        for( int i=0; i<cellPts->GetNumberOfIds(); i++)
        {
            int newPtId = newPtIds[cellPts->GetId(i)];
            newCellPtIds->InsertNextId( newPtId );
            if( this->MergeRidges )
                active |= (newPtId>=0); // at least one point is new
            else
                active &= (newPtId>=0); // all points needs to be new
        }

        if( active ){
            if( this->MergeRidges ) {
                // grab missing points from other mesh
                for( int i=0; i<newCellPtIds->GetNumberOfIds(); i++){
                    if( newCellPtIds->GetId(i) < 0)
                        newCellPtIds->SetId(i, nearestIdArray->GetValue(cellPts->GetId(i)));
                }
            }
            cells->InsertNextCell(newCellPtIds);
        }
        cellId++;
    }

    oldMesh->ShallowCopy(output);
    output->Squeeze();

    return 1;
}

