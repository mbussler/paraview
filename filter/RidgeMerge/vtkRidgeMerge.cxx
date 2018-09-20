
#include "vtkRidgeMerge.h"

#include "vtkCellData.h"
#include "vtkCompositeDataIterator.h"
#include "vtkDataSet.h"
#include "vtkDoubleArray.h"
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


vtkStandardNewMacro(vtkRidgeMerge);

vtkRidgeMerge::vtkRidgeMerge()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);

    this->CacheData = true;
    this->NumberOfCacheEntries = this->MergeRange;
}

vtkRidgeMerge::~vtkRidgeMerge()
{

}

int vtkRidgeMerge::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0) {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    }
    return 1;
}

int vtkRidgeMerge::FillOutputPortInformation(int port, vtkInformation* info)
{
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
}

int vtkRidgeMerge::RequestDataObject(vtkInformation *, 
                                       vtkInformationVector ** inputVector, 
                                       vtkInformationVector * outputVector)
{
    if (this->GetNumberOfInputPorts() == 0 || this->GetNumberOfOutputPorts() == 0)
    {
        return 1;
    }

    vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    if (!inInfo)
    {
        return 0;
    }
    vtkDataObject *input = inInfo->Get(vtkDataObject::DATA_OBJECT());

    if (input)
    {
        // for each output
        for(int i=0; i < this->GetNumberOfOutputPorts(); ++i)
        {
            vtkInformation* info = outputVector->GetInformationObject(i);
            vtkDataObject *output = info->Get(vtkDataObject::DATA_OBJECT());

            if (!output || !output->IsA(input->GetClassName()))
            {
                vtkDataObject* newOutput = input->NewInstance();
                info->Set(vtkDataObject::DATA_OBJECT(), newOutput);
                newOutput->Delete();
            }
        }
        return 1;
    }
    return 0;
}

int vtkRidgeMerge::RequestInformation( vtkInformation *vtkNotUsed(request),
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector)
{
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

    //
    // find time on input
    //
    int     numTimes = 0;
    double *inTimes  = NULL;
    double  outRange[2];

    if (inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS()))
    {
        inTimes =
            inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        numTimes =
            inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

        //outRange[0] = inTimes[0];
        //outRange[1] = inTimes[numTimes-1];
        //outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),
        //    outRange,2);
    }

    // Can we continue
    if (numTimes<this->MergeRange)
    {
        vtkErrorMacro(<<"Not enough input time steps for merge operation");
        return 0;
    }

    //outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
    //    inTimes, numTimes);

    outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

    return 1;
}

int vtkRidgeMerge::RequestUpdateExtent( vtkInformation *vtkNotUsed(request),
                                          vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector)
{
    // get the info objects
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

    // find the required input time steps and request them
    if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
    {
        // get the update times
        //double upTime =  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

        // get the available input times
        double *inTimes =
            inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        int numInTimes =
            inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

        int start = this->StartTimestep>0 ? this->StartTimestep : numInTimes;

        // only if the input is not continuous should we do anything
        if (inTimes && numInTimes >= this->MergeRange && 
            start >= this->MergeRange && start <= numInTimes )
        {
            //compute request times f
            std::vector<double> inUpTimes;
            int numInUpTimes = this->MergeRange;

            for( int i=this->MergeRange; i>0; i--) {
                inUpTimes.push_back( inTimes[start-i] );
            }

            inInfo->Set(vtkMultiTimeStepAlgorithm::UPDATE_TIME_STEPS(), &inUpTimes[0],numInUpTimes);
        }
    }

    this->NumberOfCacheEntries = this->MergeRange;

    return 1;
}

int vtkRidgeMerge::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector)
{
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkPolyData *output = vtkPolyData::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // get the requested update times
    double upTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

    vtkMultiBlockDataSet *inData = vtkMultiBlockDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    int numTimeSteps  = inData->GetNumberOfBlocks();

    // below the range
    if (numTimeSteps > 0)
    {
        std::vector<vtkPolyData*> data;

        // read in data, oldest timestep first
        for( int i=numTimeSteps-1; i>=0; i-- ){
            
            vtkPolyData* mesh = vtkPolyData::SafeDownCast(inData->GetBlock(i));
            data.push_back( mesh );

            if( mesh==NULL )
            {
                vtkErrorMacro("Null data set");
                return 0;
            }
        }

        output->ShallowCopy(data[0]);
        int numArrays = output->GetPointData()->GetNumberOfArrays();
        
        // time id array
        vtkIntArray* arr = vtkIntArray::New();
        arr->SetNumberOfComponents(1);
        arr->SetNumberOfTuples(output->GetNumberOfPoints());
        arr->SetName("Time id");


        for( vtkIdType ptId=0; ptId<output->GetNumberOfPoints(); ptId++){
            arr->SetValue(ptId, 0);
        }

        vtkSmartPointer<vtkHausdorffDistancePointSetFilter> distanceFilter =
            vtkSmartPointer<vtkHausdorffDistancePointSetFilter>::New();
        distanceFilter->SetTargetDistanceMethod(vtkHausdorffDistancePointSetFilter::POINT_TO_POINT /*POINT_TO_CELL*/);

        for( int i=1; i<numTimeSteps; i++ ){

            // calculate distance to current mesh
            vtkPolyData* mesh = data[i];
            distanceFilter->SetInputData(0, mesh );
            distanceFilter->SetInputData(1, output );
            distanceFilter->Update();

            vtkPolyData* res = distanceFilter->GetOutput(0);
            vtkDataArray* distArray = res->GetPointData()->GetArray("Distance");
            vtkIdTypeArray* nearestIdArray = vtkIdTypeArray::SafeDownCast(
                res->GetPointData()->GetArray("NearestPointId"));
            res->BuildLinks();

            std::set<vtkIdType> cellIds;
            std::vector<vtkIdType> newPtIds;
            newPtIds.resize(res->GetNumberOfPoints(), -1);

            if( distArray ){
                for( vtkIdType ptId=0; ptId<res->GetNumberOfPoints(); ptId++){
                    double dist = distArray->GetTuple1(ptId);
                    if( dist > this->MinDistance ){
                        // add point to output
                        double pos[3];
                        res->GetPoint(ptId, pos);
                        vtkIdType newId = output->GetPoints()->InsertNextPoint(pos);
                        newPtIds[ptId] = newId;
                        arr->InsertNextTuple1(i);

                        // copy point data
                        for( int i=0; i<numArrays; i++ ){
                            double data[9];
                            res->GetPointData()->GetArray(i)->GetTuple(ptId, data);
                            output->GetPointData()->GetArray(i)->InsertNextTuple(data);
                        }

                        // store cells of current point
                        vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
                        res->GetPointCells(ptId, cells);
                        for( vtkIdType i = 0; i<cells->GetNumberOfIds(); i++){
                            cellIds.insert( cells->GetId(i));
                        }
                    }
                }
            }
            
            //vtkSmartPointer<vtkCellArray> tris = vtkSmartPointer<vtkCellArray>::New();
            //tris->DeepCopy(output->GetPolys());
            vtkSmartPointer<vtkIdList> newCellPtIds;
            for( std::set<vtkIdType>::iterator iter = cellIds.begin(); iter != cellIds.end(); iter++){
                vtkIdType cellId = *iter;
                vtkCell* cell = res->GetCell(cellId);
                if( cell && cell->GetCellType() == VTK_TRIANGLE )
                {
                    vtkIdList* cellPtIds = cell->GetPointIds();
                    newCellPtIds = vtkSmartPointer<vtkIdList>::New();
                    newCellPtIds->InsertNextId( newPtIds[cellPtIds->GetId(0)]);
                    newCellPtIds->InsertNextId( newPtIds[cellPtIds->GetId(1)]);
                    newCellPtIds->InsertNextId( newPtIds[cellPtIds->GetId(2)]);

                    // if point was discarded, use nearest point in output mesh
                    //if(nearestIdArray) {
                    //    if( newCellPtIds->GetId(0) < 0 )
                    //        newCellPtIds->SetId(0, nearestIdArray->GetValue(cellPtIds->GetId(0)));
                    //    if( newCellPtIds->GetId(1) < 0 )
                    //        newCellPtIds->SetId(1, nearestIdArray->GetValue(cellPtIds->GetId(1)));
                    //    if( newCellPtIds->GetId(2) < 0 )
                    //        newCellPtIds->SetId(2, nearestIdArray->GetValue(cellPtIds->GetId(2)));
                    //}

                    //check if cell contains deleted point ids
                    if( newCellPtIds->GetId(0) > -1 &&
                        newCellPtIds->GetId(1) > -1 &&
                        newCellPtIds->GetId(2) > -1 )
                    {
                        //tris->InsertNextCell(newCellPtIds);
                        output->GetPolys()->InsertNextCell(newCellPtIds);
                    }
                }
            }
            //output->SetPolys(tris);
        }
        
        output->GetPointData()->AddArray(arr);
        arr->Delete();
    }

    // @TODO remove this when we move to new time framework
    // stamp the new temporal dataset with a time key (old style of management)
    //outData->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(),  upTime);

    //vtkSmartPointer<vtkDoubleArray> originalTimes = vtkSmartPointer<vtkDoubleArray>::New();
    //originalTimes->SetName("OriginalTimeSteps");
    //originalTimes->SetNumberOfTuples(numTimeSteps);
    //for(int i=0; i<numTimeSteps; i++)
    //{
    //    originalTimes->SetValue(i, inData->GetBlock(i)->GetInformation()->Get(vtkDataObject::DATA_TIME_STEP()));
    //}
    //outData->GetFieldData()->AddArray(originalTimes);

    return 1;
}

