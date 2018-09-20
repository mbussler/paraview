/*=========================================================================

Program:   Visualization Toolkit
Module:    NearestNeighborSampling.h

Copyright (c) Michael Bußler

=========================================================================*/

#include "NearestNeighborSampling.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkDataSet.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkTransform.h"
#include "vtkSmartPointer.h"
#include "vtkMath.h"
#include "vtkPolyData.h"
#include "vtkImageData.h"
#include "vtkKdTreePointLocator.h"

#include "linalg.h"

#include <stdio.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>

vtkStandardNewMacro(NearestNeighborSampling);


//==============================================================================
bool NearestNeighborSamplingFilterHasArray(vtkFieldData *fieldData,
    vtkDataArray *array)
{
    int numarrays = fieldData->GetNumberOfArrays();
    for (int i = 0; i < numarrays; i++)
    {
        if (array == fieldData->GetArray(i))
        {
            return true;
        }
    }
    return false;
}

//==============================================================================
NearestNeighborSampling::NearestNeighborSampling()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

    this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, vtkDataSetAttributes::SCALARS);

    this->ModelBounds[0] = -1.0;
    this->ModelBounds[1] =  1.0;
    this->ModelBounds[2] = -1.0;
    this->ModelBounds[3] =  1.0;
    this->ModelBounds[4] = -1.0;
    this->ModelBounds[5] =  1.0;

    this->SampleDimensions[0] = 50;
    this->SampleDimensions[1] = 50;
    this->SampleDimensions[2] = 50;
}
//==============================================================================
NearestNeighborSampling::~NearestNeighborSampling()
{
}

//-----------------------------------------------------------------------------
void NearestNeighborSampling::SetInputScalars(int fieldAssociation, const char *name)
{
    if (   (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_POINTS))
    {
        vtkErrorMacro("Input Array must be associated with points.");
        return;
    }

    this->SetInputArrayToProcess(0, 0, 0, fieldAssociation, name);
}

//-----------------------------------------------------------------------------
void NearestNeighborSampling::SetInputScalars(int fieldAssociation,
    int fieldAttributeType)
{
    if (   (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_POINTS) )
    {
        vtkErrorMacro("Input Array must be associated with points.");
        return;
    }

    this->SetInputArrayToProcess(0, 0, 0, fieldAssociation, fieldAttributeType);
}

//==============================================================================
int NearestNeighborSampling::RequestInformation ( 
    vtkInformation * vtkNotUsed(request),
    vtkInformationVector ** inputVector,
    vtkInformationVector *outputVector)
{
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation* outInfo = outputVector->GetInformationObject(0);

    vtkDataSet *input = vtkDataSet::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));

    // use input bounds for sampling
    input->GetBounds( this->ModelBounds);

    int i;
    double ar[3], origin[3];

    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(),
        0, this->SampleDimensions[0]-1,
        0, this->SampleDimensions[1]-1,
        0, this->SampleDimensions[2]-1);

    for (i=0; i < 3; i++)
    {
        origin[i] = this->ModelBounds[2*i];
        if ( this->SampleDimensions[i] <= 1 )
        {
            ar[i] = 1;
        }
        else
        {
            ar[i] = (this->ModelBounds[2*i+1] - this->ModelBounds[2*i])
                / (this->SampleDimensions[i] - 1);
        }
    }
    outInfo->Set(vtkDataObject::ORIGIN(),origin,3);
    outInfo->Set(vtkDataObject::SPACING(),ar,3);

    return 1;
}

//==============================================================================
int NearestNeighborSampling::RequestData(
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
	vtkImageData *output = vtkImageData::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkIdType numPts, inPtId, inCellId;

	vtkDebugMacro(<<"Resampling using Nearest Neighbor");

	numPts = input->GetNumberOfPoints();

	if ( numPts < 1 )
	{
		vtkErrorMacro(<<"No data!");
		return 1;
	}

	vtkDebugMacro(<<"Initialise..");

    // Evaluate function on rectangular grid
    //

	// init output
	output->SetExtent( 0, this->SampleDimensions[0]-1,
		0, this->SampleDimensions[1]-1,
		0, this->SampleDimensions[2]-1 );

	vtkIdType numNewPts = output->GetNumberOfPoints();

    // copy point data arrays
    vtkPointData *inPD=input->GetPointData();
    vtkPointData *outPD = output->GetPointData();
    int numArrays = inPD->GetNumberOfArrays();
    for( int i=0; i<numArrays; i++ ) {

        int type = inPD->GetArray(i)->GetDataType();
        vtkDataArray *arr = vtkDataArray::CreateDataArray(type);
        arr->SetNumberOfComponents(inPD->GetArray(i)->GetNumberOfComponents());
        arr->SetName(inPD->GetArray(i)->GetName());
        arr->SetNumberOfTuples(numNewPts);

        outPD->AddArray(arr);
        arr->Delete();
    }

	// point location
	vtkSmartPointer<vtkKdTreePointLocator> tree = vtkSmartPointer<vtkKdTreePointLocator>::New();
	tree->SetDataSet(input);
	tree->BuildLocator();

	for( vtkIdType ptId=0; ptId<numNewPts; ptId++) {

		double pos[3];
		output->GetPoint(ptId, pos);
		vtkIdType nn = tree->FindClosestPoint(pos);

        // copy point data
        for( int i=0; i<numArrays; i++ ){
            double data[9];
            inPD->GetArray(i)->GetTuple(nn, data);
            outPD->GetArray(i)->SetTuple(ptId, data);
        }

		this->UpdateProgress((ptId/(float)numNewPts));
	}

	return 1;
}

//==============================================================================
void NearestNeighborSampling::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);

    os << indent << "Sample Dimensions: (" << this->SampleDimensions[0] << ", "
        << this->SampleDimensions[1] << ", "
        << this->SampleDimensions[2] << ")\n";

    os << indent << "ModelBounds: \n";
    os << indent << "  Xmin,Xmax: (" << this->ModelBounds[0] << ", " << this->ModelBounds[1] << ")\n";
    os << indent << "  Ymin,Ymax: (" << this->ModelBounds[2] << ", " << this->ModelBounds[3] << ")\n";
    os << indent << "  Zmin,Zmax: (" << this->ModelBounds[4] << ", " << this->ModelBounds[5] << ")\n";
}

int NearestNeighborSampling::ProcessRequest(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	return vtkImageAlgorithm::ProcessRequest(request, inputVector, outputVector);
}

int NearestNeighborSampling::FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
}

// Set the i-j-k dimensions on which to sample the distance function.
void NearestNeighborSampling::SetSampleDimensions(int i, int j, int k)
{
    int dim[3];

    dim[0] = i;
    dim[1] = j;
    dim[2] = k;

    this->SetSampleDimensions(dim);
}

// Set the i-j-k dimensions on which to sample the distance function.
void NearestNeighborSampling::SetSampleDimensions(int dim[3])
{
    int dataDim, i;

    vtkDebugMacro(<< " setting SampleDimensions to (" << dim[0] << ","
        << dim[1] << "," << dim[2] << ")");

    if ( dim[0] != this->SampleDimensions[0] ||
        dim[1] != this->SampleDimensions[1] ||
        dim[2] != this->SampleDimensions[2] )
    {
        if ( dim[0]<1 || dim[1]<1 || dim[2]<1 )
        {
            vtkErrorMacro (<< "Bad Sample Dimensions, retaining previous values");
            return;
        }

        for (dataDim=0, i=0; i<3 ; i++)
        {
            if (dim[i] > 1)
            {
                dataDim++;
            }
        }

        if ( dataDim  < 3 )
        {
            vtkErrorMacro(<<"Sample dimensions must define a volume!");
            return;
        }

        for ( i=0; i<3; i++)
        {
            this->SampleDimensions[i] = dim[i];
        }

        this->Modified();
    }
}


