
#include "vtkRidgeCompare.h"

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
#include "vtkDistancePolyDataFilter.h"
#include "vtkHausdorffDistancePointSetFilter.h"


vtkStandardNewMacro(vtkRidgeCompare);

vtkRidgeCompare::vtkRidgeCompare()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(2);

    this->CacheData = true;
    this->NumberOfCacheEntries = 2;
}

vtkRidgeCompare::~vtkRidgeCompare()
{

}

int vtkRidgeCompare::FillInputPortInformation(int port, vtkInformation* info)
{
    if (port == 0) {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    }
    else if (port == 1) {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    }
    return 1;
}

int vtkRidgeCompare::FillOutputPortInformation(int port, vtkInformation* info)
{
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
}

int vtkRidgeCompare::RequestDataObject(vtkInformation *, 
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

int vtkRidgeCompare::RequestInformation( vtkInformation *vtkNotUsed(request),
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

        outRange[0] = inTimes[0];
        outRange[1] = inTimes[numTimes-1];
        outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),
            outRange,2);
    }

    // Can we continue
    if (numTimes<2)
    {
        vtkErrorMacro(<<"Not enough input time steps for interpolation");
        return 0;
    }

    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
        inTimes, numTimes);

    return 1;
}

int vtkRidgeCompare::RequestUpdateExtent( vtkInformation *vtkNotUsed(request),
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
        double upTime =  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

        // get the available input times
        double *inTimes =
            inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        int numInTimes =
            inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

        // only if the input is not continuous should we do anything
        if (inTimes)
        {
            //compute request times f
            double inUpTimes[2];
            int numInUpTimes(0);

            // for each requested time mark the required input times
            numInUpTimes= 0;
            // below the range
            if (upTime <= inTimes[0])
            {
                inUpTimes[numInUpTimes++] = inTimes[0];
            }
            // above the range?
            else if (upTime > inTimes[numInTimes-1])
            {
                inUpTimes[numInUpTimes++] = inTimes[numInTimes-1];
            }
            // in the middle
            else
            {
                int i = 0;
                while (upTime > inTimes[i])
                {
                    ++i;
                }
                inUpTimes[numInUpTimes++] = inTimes[i-1];
                inUpTimes[numInUpTimes++] = inTimes[i];
            }

            inInfo->Set(vtkMultiTimeStepAlgorithm::UPDATE_TIME_STEPS(), inUpTimes,numInUpTimes);
        }
    }
    return 1;
}

int vtkRidgeCompare::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector)
{
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkInformation *outInfo1 = outputVector->GetInformationObject(1);
    vtkDataObject *outData = vtkDataObject::GetData(outInfo);
    vtkPolyData   *outRef  = vtkPolyData::SafeDownCast(outInfo1->Get(vtkDataObject::DATA_OBJECT()));

    // get the requested update times
    double upTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

    vtkMultiBlockDataSet *inData = vtkMultiBlockDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    int numTimeSteps  = inData->GetNumberOfBlocks();

    // below the range
    if (numTimeSteps==1)
    {
        // pass the lowest data
        outData->ShallowCopy(inData->GetBlock(0));
    }
    else
    {
        vtkDataObject* data0 = inData->GetBlock(0);
        vtkDataObject* data1 = inData->GetBlock(1);
        if (data0==NULL && data1==NULL)
        {
            vtkErrorMacro("Null data set");
            return 0;
        }
        // interpolate i-1 and i
        double t0 = data0->GetInformation()->Get(vtkDataObject::DATA_TIME_STEP());
        double t1 = data1->GetInformation()->Get(vtkDataObject::DATA_TIME_STEP());
        double Ratio  = (upTime-t0)/(t1 - t0);
        
        vtkPolyData *inds0 = vtkPolyData::SafeDownCast(data0);
        vtkPolyData *inds1 = vtkPolyData::SafeDownCast(data1);

        //vtkSmartPointer<vtkDistancePolyDataFilter> distanceFilter =
        //    vtkSmartPointer<vtkDistancePolyDataFilter>::New();

        vtkSmartPointer<vtkHausdorffDistancePointSetFilter> distanceFilter =
            vtkSmartPointer<vtkHausdorffDistancePointSetFilter>::New();

        distanceFilter->SetInputData( 0, inds1);
        distanceFilter->SetInputData( 1, inds0 );
        //distanceFilter->SignedDistanceOff();
        //distanceFilter->ComputeSecondDistanceOn();
        distanceFilter->SetTargetDistanceMethod(this->TargetDistanceMethod);

        distanceFilter->Update();

        outData->ShallowCopy(distanceFilter->GetOutput(0));
        outRef->ShallowCopy(distanceFilter->GetOutput(1));
        
        vtkDataArray* distRefArray = outRef->GetPointData()->GetArray("Distance");
        if( distRefArray )
            distRefArray->SetName("Distance_ref");
    }

    // @TODO remove this when we move to new time framework
    // stamp the new temporal dataset with a time key (old style of management)
    outData->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(),  upTime);

    vtkSmartPointer<vtkDoubleArray> originalTimes = vtkSmartPointer<vtkDoubleArray>::New();
    originalTimes->SetName("OriginalTimeSteps");
    originalTimes->SetNumberOfTuples(numTimeSteps);
    for(int i=0; i<numTimeSteps; i++)
    {
        originalTimes->SetValue(i, inData->GetBlock(i)->GetInformation()->Get(vtkDataObject::DATA_TIME_STEP()));
    }
    outData->GetFieldData()->AddArray(originalTimes);

    return 1;
}

