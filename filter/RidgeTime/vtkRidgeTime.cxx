
#include "vtkRidgeTime.h"

#include "vtkCellData.h"
#include "vtkCompositeDataIterator.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkCompositeDataPipeline.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkTransform.h"
#include "vtkMath.h"


vtkStandardNewMacro(vtkRidgeTime);

void vtkRidgeTime::SetMergeRange(int in, int out)
{
    int range[2];
    range[0] = in;
    range[1] = out;
    SetMergeRange(range);
}

void vtkRidgeTime::SetMergeRange(int range[2])
{
    if( range[0] != this->MergeRange[0] ||
        range[1] != this->MergeRange[1] )
    {
        if( range[0] < 0)
            range[0] = 0;
        if( range[1] < 0)
            range[1] = 0;

        if( range[0] > range[1] )
            range[0] = range[1];

        this->MergeRange[0] = range[0];
        this->MergeRange[1] = range[1];
    }
}

void vtkRidgeTime::SetTranslation(double x, double y, double z)
{
    double translation[3];
    translation[0] = x;
    translation[1] = y;
    translation[2] = z;
    this->SetTranslation(translation);
}

void vtkRidgeTime::SetTranslation(double xyz[3])
{
    if( xyz[0] != this->Translation[0] ||
        xyz[1] != this->Translation[1] ||
        xyz[2] != this->Translation[2] )
    {
        for( int i=0; i<3; i++ ){
            this->Translation[i] = xyz[i];
        }

        this->Modified();
    }
}

vtkRidgeTime::vtkRidgeTime()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);

    this->CacheData = true;
    this->NumberOfCacheEntries = 0;

    this->Translation[0] = 0.0;
    this->Translation[1] = 0.0;
    this->Translation[2] = 0.0;

    this->MergeRange[0] = 0;
    this->MergeRange[1] = 2;
}

vtkRidgeTime::~vtkRidgeTime()
{

}

int vtkRidgeTime::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0) {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    }
    return 1;
}

int vtkRidgeTime::FillOutputPortInformation(int port, vtkInformation* info)
{
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
}

int vtkRidgeTime::RequestDataObject(vtkInformation *, 
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

int vtkRidgeTime::RequestInformation( vtkInformation *vtkNotUsed(request),
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

    if (inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS()))
    {
        inTimes =
            inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        numTimes =
            inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    }

    // Can we continue
    if (numTimes<(this->MergeRange[1]-this->MergeRange[0]+1))
    {
        vtkErrorMacro(<<"Not enough input time steps for merge operation");
        return 0;
    }

    // we want time static output
    outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

    return 1;
}

int vtkRidgeTime::RequestUpdateExtent( vtkInformation *vtkNotUsed(request),
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

        // only if the input is not continuous should we do anything
        if (inTimes && 
            this->MergeRange[0] >= 0 &&
            this->MergeRange[1] < numInTimes )
        {
            //compute request times f
            std::vector<double> inUpTimes;
            int numInUpTimes = this->MergeRange[1]-this->MergeRange[0]+1;

            for( int i=this->MergeRange[0]; i<=this->MergeRange[1]; i++) {
                inUpTimes.push_back( inTimes[i] );
            }

            inInfo->Set(vtkMultiTimeStepAlgorithm::UPDATE_TIME_STEPS(), &inUpTimes[0],numInUpTimes);

            this->NumberOfCacheEntries = numInUpTimes;
        }
    }

    return 1;
}

int vtkRidgeTime::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector)
{
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // get the requested update times
    double upTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

    vtkMultiBlockDataSet *inData = vtkMultiBlockDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    int numTimeSteps  = inData->GetNumberOfBlocks();

    // below the range
    if (numTimeSteps > 0)
    {
        std::vector<vtkUnstructuredGrid*> data;

        // read in data, oldest timestep last
        for( int i=0; i<numTimeSteps; i++ ){
            
            vtkUnstructuredGrid* mesh = vtkUnstructuredGrid::SafeDownCast(inData->GetBlock(i));
            data.push_back( mesh );

            if( mesh==NULL )
            {
                vtkErrorMacro("Null data set");
                return 0;
            }
        }

        //output->ShallowCopy(data[0]);
        //vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
        //vtkSmartPointer<vtkCellArray> newCells = vtkSmartPointer<vtkCellArray>::New();

        output->SetPoints(vtkPoints::New());
        output->Allocate();

        // time id array
        vtkIntArray* arr = vtkIntArray::New();
        arr->SetNumberOfComponents(1);
        arr->SetNumberOfTuples(output->GetNumberOfPoints());
        arr->SetName("Time id");

        // copy point data arrays
        vtkPointData *inPD=data[0]->GetPointData(), 
                     *outPD = output->GetPointData();

        int numArrays = inPD->GetNumberOfArrays();
        for( int i=0; i<numArrays; i++ ){
            int type = inPD->GetArray(i)->GetDataType();
            vtkDataArray *arr = vtkDataArray::CreateDataArray(type);
            arr->SetNumberOfComponents(inPD->GetArray(i)->GetNumberOfComponents());
            arr->SetName(inPD->GetArray(i)->GetName());
            outPD->AddArray(arr);
            arr->Delete();
        }

        // copy cell data arrays
        vtkCellData *inCD=data[0]->GetCellData(), 
                    *outCD = output->GetCellData();
        int numCellArrays = inCD->GetNumberOfArrays();
        for( int i=0; i<numCellArrays; i++ ){
            int type = inCD->GetArray(i)->GetDataType();
            vtkDataArray *arr = vtkDataArray::CreateDataArray(type);
            arr->SetNumberOfComponents(inCD->GetArray(i)->GetNumberOfComponents());
            arr->SetName(inCD->GetArray(i)->GetName());
            outCD->AddArray(arr);
            arr->Delete();
        }

        //for( vtkIdType ptId=0; ptId<output->GetNumberOfPoints(); ptId++){
        //    arr->SetValue(ptId, 0);
        //}

        double transVec[3] = { 0,0,0 };

        for( int i=0; i<numTimeSteps; i++ ){

            // calculate distance to current mesh
            vtkUnstructuredGrid* dataSet = data[i];

            //int numArrays = dataSet->GetPointData()->GetNumberOfArrays();
            //int numCellArrays = dataSet->GetCellData()->GetNumberOfArrays();
            
            std::vector<vtkIdType> newPtIds;
            newPtIds.resize( dataSet->GetNumberOfPoints(), -1);

            for( vtkIdType ptId=0; ptId<dataSet->GetNumberOfPoints(); ptId++)
            {
                // add point to output
                double pos[3];
                dataSet->GetPoint(ptId, pos);
                pos[0] += transVec[0];
                pos[1] += transVec[1];
                pos[2] += transVec[2];
                vtkIdType newId = output->GetPoints()->InsertNextPoint(pos);
                newPtIds[ptId] = newId;
                arr->InsertNextTuple1(i);

                // copy point data
                for( int i=0; i<numArrays; i++ ){
                    double data[9];
                    dataSet->GetPointData()->GetArray(i)->GetTuple(ptId, data);
                    output->GetPointData()->GetArray(i)->InsertNextTuple(data);
                }
            }
            
            vtkMath::Add(transVec, this->Translation, transVec);

            //vtkSmartPointer<vtkCellArray> cells = output->GetCells();
            //vtkCellArray* cells_ref = dataSet->GetCells();
            //cells_ref->InitTraversal();

            //vtkSmartPointer<vtkIdList> cellPts = vtkSmartPointer<vtkIdList>::New();
            //vtkIdType cellid = 0;
            //while( cells_ref->GetNextCell(cellPts))
            for( vtkIdType cellId=0; cellId<dataSet->GetNumberOfCells(); cellId++)
            {
                vtkSmartPointer<vtkIdList> newCellPtIds = vtkSmartPointer<vtkIdList>::New();
                vtkCell* cell = dataSet->GetCell(cellId);

                if( cell ) {
                    for( int i=0; i<cell->GetNumberOfPoints(); i++)
                    {
                        int newPtId = newPtIds[cell->GetPointId(i)];
                        newCellPtIds->InsertNextId( newPtId );
                    }

                    //newCells->InsertNextCell(newCellPtIds);
                    //cellTypes->InsertNextValue(cell->GetCellType());
                    //cells->InsertNextCell(newCellPtIds);
                    output->InsertNextCell(cell->GetCellType(), newCellPtIds);

                    for( int i=0; i<numCellArrays; i++ ){
                        double data[9];
                        dataSet->GetCellData()->GetArray(i)->GetTuple(cellId, data);
                        output->GetCellData()->GetArray(i)->InsertNextTuple(data);
                    }
                }

                //for( int i=0; i<cellPts->GetNumberOfIds(); i++)
                //{
                //    int newPtId = newPtIds[cellPts->GetId(i)];
                //    newCellPtIds->InsertNextId( newPtId );
                //}

                //cells->InsertNextCell(newCellPtIds);

                // copy cell data
                //for( int i=0; i<numCellArrays; i++ ){
                //    double data[9];
                //    dataSet->GetCellData()->GetArray(i)->GetTuple(cellid, data);
                //    output->GetCellData()->GetArray(i)->InsertNextTuple(data);
                //}

                //cellid++;
            }
            
            this->SetProgress(i/(double)numTimeSteps);
        }
        
        output->GetPointData()->AddArray(arr);
        output->Squeeze();
        arr->Delete();
    }

    return 1;
}

