/*=========================================================================

Program:   Visualization Toolkit
Module:    PeriMultiFit.h

Copyright (c) Michael Bußler

=========================================================================*/

#include "PeriMultiFit.h"
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

vtkStandardNewMacro(PeriMultiFit);

inline bool SameSign(float a, float b) {
	return a*b >= 0.0f;
}

inline void swap( int&a, int&b ){
	int swp=a; a=b; b=swp;
};

enum FitFuncOrder {
    Linear = 0,
    Quadratic = 1,
	NoHyperbolic = 2,
    maxCoeff = 10 // number of coefficients needed for highest order
};


//==============================================================================
bool PeriMultiFitFilterHasArray(vtkFieldData *fieldData,
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
PeriMultiFit::PeriMultiFit()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

    this->SetInputScalars(vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,
        vtkDataSetAttributes::SCALARS);

    this->ModelBounds[0] = -1.0;
    this->ModelBounds[1] =  1.0;
    this->ModelBounds[2] = -1.0;
    this->ModelBounds[3] =  1.0;
    this->ModelBounds[4] = -1.0;
    this->ModelBounds[5] =  1.0;

    this->SampleDimensions[0] = 50;
    this->SampleDimensions[1] = 50;
    this->SampleDimensions[2] = 50;

	this->FitConstantTerm		= true;
	this->FitLinearTerms		= true;
	this->FitHyperbolicTerms 	= true;
	this->FitQuadraticTerms		= true;
}
//==============================================================================
PeriMultiFit::~PeriMultiFit()
{
}

//-----------------------------------------------------------------------------
void PeriMultiFit::SetInputScalars(int fieldAssociation, const char *name)
{
    if (   (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_POINTS)
        && (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_CELLS)
        && (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS) )
    {
        vtkErrorMacro("Input Array must be associated with points or cells.");
        return;
    }

    this->SetInputArrayToProcess(0, 0, 0, fieldAssociation, name);
}

//-----------------------------------------------------------------------------
void PeriMultiFit::SetInputScalars(int fieldAssociation,
    int fieldAttributeType)
{
    if (   (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_POINTS)
        && (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_CELLS)
        && (fieldAssociation != vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS) )
    {
        vtkErrorMacro("Input Array must be associated with points or cells.");
        return;
    }

    this->SetInputArrayToProcess(0, 0, 0, fieldAssociation, fieldAttributeType);
}

//==============================================================================
int PeriMultiFit::RequestInformation ( 
    vtkInformation * vtkNotUsed(request),
    vtkInformationVector ** inputVector,
    vtkInformationVector *outputVector)
{
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation* outInfo = outputVector->GetInformationObject(0);

    vtkPolyData *input = vtkPolyData::SafeDownCast(
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
int PeriMultiFit::RequestData(
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
	vtkImageData *output = vtkImageData::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkDataArray *array = this->GetInputArrayToProcess(0, inputVector);

    if (array == NULL)
    {
        vtkErrorMacro("No input array.");
        return 0;
    }
    if (array->GetNumberOfComponents() == 0)
    {
        vtkErrorMacro("Input array must have at least one component.");
        return 0;
    }

    int fieldAssociation;
    if (PeriMultiFitFilterHasArray(input->GetPointData(), array))
    {
        fieldAssociation = vtkDataObject::FIELD_ASSOCIATION_POINTS;
    }
    else if (PeriMultiFitFilterHasArray(input->GetCellData(), array))
    {
        fieldAssociation = vtkDataObject::FIELD_ASSOCIATION_CELLS;
    }
    else
    {
        vtkErrorMacro("Input arrays do not seem to be either point or cell arrays.");
        return 0;
    }

	vtkIdType numPts, inPtId, inCellId;

	vtkDebugMacro(<<"Calculating Multi Quadratic Least Squares");

	numPts = input->GetNumberOfPoints();

	if ( numPts < 1 )
	{
		vtkErrorMacro(<<"No data!");
		return 1;
	}

	vtkDebugMacro(<<"Initialise..");

	bool functionChecked = FitConstantTerm || FitLinearTerms || FitHyperbolicTerms || FitQuadraticTerms;
	if( !functionChecked )
	{
		vtkErrorMacro(<<"No terms for function checked!");
		return 1;
	}

    input->BuildLinks(); // build up connectivity lookup table

	// allocate mem for output data arrays
	// array stores intensity per point
	vtkSmartPointer<vtkDoubleArray> coeffs;

	coeffs = vtkDoubleArray::New();
	coeffs->SetName("Coefficients");

	int numCoeffs =   FitConstantTerm    * 1
					+ FitLinearTerms     * 3
					+ FitHyperbolicTerms * 3
					+ FitQuadraticTerms  * 3;

    coeffs->SetNumberOfComponents(numCoeffs);

	vtkDebugMacro(<<"Calculate PCA..");

    // allocate mem for result
    coeffs->SetNumberOfTuples(numPts);

	// calculate multi quadratic fitting per particle
	for (inPtId=0; inPtId < numPts; inPtId++) {

		double pos[3], center[3];
		input->GetPoint(inPtId, center);

		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkDoubleArray> values = vtkSmartPointer<vtkDoubleArray>::New();
		//points->InsertNextPoint(pos);

        const vtkSmartPointer<vtkIdList>& bondedPoints = 
                GetConnectedVertices( inPtId, input, fieldAssociation, array, values);

        for( vtkIdType i=0; i<bondedPoints->GetNumberOfIds(); i++ )
        {
            input->GetPoint( bondedPoints->GetId(i), pos);            
            vtkMath::Subtract(pos,center,pos); // shift to origin
            points->InsertNextPoint(pos);
        }

		double coeff[maxCoeff];
		LeastSquaresFitting(points, values, coeff);
		coeffs->SetTuple(inPtId, coeff);

        this->UpdateProgress(0.8*(inPtId/(float)numPts));
	}

    // Evaluate fitted function on rectangular grid
    //
    vtkIdType numNewPts = this->SampleDimensions[0] * this->SampleDimensions[1] * this->SampleDimensions[2];

    vtkDoubleArray* newScalars = vtkDoubleArray::New();
    newScalars->SetNumberOfComponents(1);
    newScalars->SetNumberOfTuples(numNewPts);
    newScalars->SetName( "Fitted Function" );

    vtkDoubleArray* gradient = vtkDoubleArray::New();
    gradient->SetNumberOfComponents(3);
    gradient->SetNumberOfTuples(numNewPts);
    gradient->SetName( "Gradient" );

    // init output
    output->SetExtent( 0, this->SampleDimensions[0]-1,
					   0, this->SampleDimensions[1]-1,
					   0, this->SampleDimensions[2]-1 );

    // point location
    vtkSmartPointer<vtkKdTreePointLocator> tree = vtkSmartPointer<vtkKdTreePointLocator>::New();
    tree->SetDataSet(input);
    tree->BuildLocator();

    for( vtkIdType ptId=0; ptId<numNewPts; ptId++) {

		double pos[3], center[3];
		double coeff[maxCoeff];
        
		output->GetPoint(ptId, pos);
        vtkIdType nn = tree->FindClosestPoint(pos);

        coeffs->GetTuple(nn, coeff);
		input->GetPoint(nn, center);
		vtkMath::Subtract(pos,center,pos);

        double val;
		double grad[3];

		val = evaluateFunction(pos, coeff);
		evaluateGradient(pos, coeff, grad);

        newScalars->SetValue(ptId, val);
		gradient->SetTuple(ptId, grad);

        this->UpdateProgress(0.8+0.2*(ptId/(float)numNewPts));

    }

	output->GetPointData()->AddArray(newScalars);
    output->GetPointData()->AddArray(gradient);

	return 1;
}

//==============================================================================
vtkSmartPointer<vtkIdList> PeriMultiFit::GetConnectedVertices( 
    int id, vtkSmartPointer<vtkPolyData> mesh, 
    int fieldAssocation, vtkDataArray *array, 
    vtkSmartPointer<vtkDoubleArray> values )
{
	vtkSmartPointer<vtkIdList> connectedVertices =
		vtkSmartPointer<vtkIdList>::New();

	//get all cells that vertex 'id' is a part of
	vtkSmartPointer<vtkIdList> cellIdList =
		vtkSmartPointer<vtkIdList>::New();

    mesh->GetPointCells(id, cellIdList);

	for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
	{
		vtkSmartPointer<vtkIdList> pointIdList =
			vtkSmartPointer<vtkIdList>::New();

        vtkIdType cellId = cellIdList->GetId(i);

		mesh->GetCellPoints(cellId, pointIdList);

        vtkIdType bondedId = pointIdList->GetId(0);
        if( bondedId == id)
            bondedId = pointIdList->GetId(1);

        connectedVertices->InsertNextId(bondedId);

        double value = 0.0;
        if( fieldAssocation == vtkDataObject::FIELD_ASSOCIATION_POINTS ) {
            value = array->GetComponent(bondedId, 0);
        } else /*if (fieldAssocation == vtkDataObject::FIELD_ASSOCIATION_CELLS)*/ {
            value = array->GetComponent(cellId, 0);
        }
        values->InsertNextValue(value);
	}

	return connectedVertices;
}



void PeriMultiFit::LeastSquaresFitting(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkDoubleArray> values, double coeff[])
{
	double chisq;
	gsl_matrix *X, *cov;
	gsl_vector *v, *w, *c;

	vtkIdType nPoints = points->GetNumberOfPoints();
	//float wt[nPoints];

	size_t conOff =   0;
	size_t linOff =   conOff + FitConstantTerm    * 1;
	size_t hypOff =   linOff + FitLinearTerms     * 3;
	size_t quaOff =   hypOff + FitHyperbolicTerms * 3;
	size_t nCoeff =   quaOff + FitQuadraticTerms  * 3;

	X = gsl_matrix_alloc (nPoints, nCoeff);
	v = gsl_vector_alloc (nPoints);
	w = gsl_vector_alloc (nPoints);

	c = gsl_vector_alloc (nCoeff);
	cov = gsl_matrix_alloc (nCoeff, nCoeff);

	for (size_t i = 0; i<nPoints ; i++)  {

		double pos[3];
		points->GetPoint(i, pos);

		// constant term
		if( this->FitConstantTerm)
			gsl_matrix_set (X, i, 0, 1.0);

		// linear terms
		if( this->FitLinearTerms )
		{
			gsl_matrix_set (X, i, linOff,   pos[0]);
			gsl_matrix_set (X, i, linOff+1, pos[1]);
			gsl_matrix_set (X, i, linOff+2, pos[2]);
		}

		// mixed terms
		if( this->FitHyperbolicTerms) {
			gsl_matrix_set (X, i, hypOff,   pos[0]*pos[1]);
			gsl_matrix_set (X, i, hypOff+1, pos[0]*pos[2]);
			gsl_matrix_set (X, i, hypOff+2, pos[1]*pos[2]);
		}

		// quadratic terms
		if( this->FitQuadraticTerms ){
			gsl_matrix_set (X, i, quaOff,   pos[0]*pos[0]);
			gsl_matrix_set (X, i, quaOff+1, pos[1]*pos[1]);
			gsl_matrix_set (X, i, quaOff+2, pos[2]*pos[2]);
		}

		double d = values->GetValue(i);
		gsl_vector_set (v, i, d);   // function value at x,y,z
		gsl_vector_set (w, i, 1.0); // weights
	}

	gsl_multifit_linear_workspace * work 
		= gsl_multifit_linear_alloc (nPoints, nCoeff);
	gsl_multifit_wlinear (X, w, v, c, cov,
		&chisq, work);
	gsl_multifit_linear_free (work);


	for (size_t i=0; i<nCoeff; i++) {
		coeff[i] = gsl_vector_get(c, i);
	}

	gsl_matrix_free (X);
	gsl_vector_free (v);
	gsl_vector_free (w);
	gsl_vector_free (c);
	gsl_matrix_free (cov);
}

double PeriMultiFit::evaluateFunction(double pos[], double c[])
{
	double x = pos[0],
		   y = pos[1],
		   z = pos[2];

	size_t conOff =   0;
	size_t linOff =   conOff + FitConstantTerm    * 1;
	size_t hypOff =   linOff + FitLinearTerms     * 3;
	size_t quaOff =   hypOff + FitHyperbolicTerms * 3;
	size_t nCoeff =   quaOff + FitQuadraticTerms  * 3;

	return    FitConstantTerm    * c[conOff]
			+ FitLinearTerms     * c[linOff]*x   + c[linOff+1]*y   + c[linOff+2]*z  
			+ FitHyperbolicTerms * c[hypOff]*x*y + c[hypOff+1]*x*z + c[hypOff+2]*y*z
			+ FitQuadraticTerms  * c[quaOff]*x*x + c[quaOff+1]*y*y + c[quaOff+2]*z*z;
}

void PeriMultiFit::evaluateGradient(double pos[], double c[], double grad[])
{
	double x = pos[0],
		   y = pos[1],
		   z = pos[2];

	size_t conOff =   0;
	size_t linOff =   conOff + FitConstantTerm    * 1;
	size_t hypOff =   linOff + FitLinearTerms     * 3;
	size_t quaOff =   hypOff + FitHyperbolicTerms * 3;
	size_t nCoeff =   quaOff + FitQuadraticTerms  * 3;

	grad[0] =   (FitLinearTerms)     * c[linOff+0]                     
		      + (FitHyperbolicTerms) * c[hypOff+0]*y  + c[hypOff+1]*z 
			  + (FitQuadraticTerms)  * c[quaOff+0]*x                 ;
									 
	grad[1] =   (FitLinearTerms)     * c[linOff+1]                     
			  + (FitHyperbolicTerms) * c[hypOff+0]*x  + c[hypOff+2]*z 
			  + (FitQuadraticTerms)  * c[quaOff+1]*y                 ;
									 
	grad[2] =   (FitLinearTerms)     * c[linOff+2]                     
		      + (FitHyperbolicTerms) * c[hypOff+1]*x  + c[hypOff+2]*y 
		      + (FitQuadraticTerms)  * c[quaOff+2]*z                 ;

}

//==============================================================================
void PeriMultiFit::PrintSelf(ostream &os, vtkIndent indent)
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

int PeriMultiFit::ProcessRequest(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	return vtkImageAlgorithm::ProcessRequest(request, inputVector, outputVector);
}

int PeriMultiFit::FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
}

// Set the i-j-k dimensions on which to sample the distance function.
void PeriMultiFit::SetSampleDimensions(int i, int j, int k)
{
    int dim[3];

    dim[0] = i;
    dim[1] = j;
    dim[2] = k;

    this->SetSampleDimensions(dim);
}

// Set the i-j-k dimensions on which to sample the distance function.
void PeriMultiFit::SetSampleDimensions(int dim[3])
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


