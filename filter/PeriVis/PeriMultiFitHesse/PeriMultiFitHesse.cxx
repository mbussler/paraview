/*=========================================================================

Program:   Visualization Toolkit
Module:    PeriMultiFitHesse.h

Copyright (c) Michael Bußler

=========================================================================*/

#include "PeriMultiFitHesse.h"
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
#include "vtkUnstructuredGrid.h"
#include "vtkKdTreePointLocator.h"
#include <vtkDelaunay3D.h>

#include "linalg.h"

#include <stdio.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>

vtkStandardNewMacro(PeriMultiFitHesse);

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

#define H(x,y) H_f_a[3*x+y]

//==============================================================================
bool PeriMultiFitHesseFilterHasArray(vtkFieldData *fieldData,
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
PeriMultiFitHesse::PeriMultiFitHesse()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);

    this->SetInputScalars(vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,
        vtkDataSetAttributes::SCALARS);

    this->FitConstantTerm		= true;
    this->FitLinearTerms		= true;
    this->FitHyperbolicTerms 	= true;
    this->FitQuadraticTerms		= true;

    this->Interpolate = true;
}
//==============================================================================
PeriMultiFitHesse::~PeriMultiFitHesse()
{
}

//-----------------------------------------------------------------------------
void PeriMultiFitHesse::SetInputScalars(int fieldAssociation, const char *name)
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
void PeriMultiFitHesse::SetInputScalars(int fieldAssociation,
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
int PeriMultiFitHesse::RequestData(
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
    vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
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
    if (PeriMultiFitHesseFilterHasArray(input->GetPointData(), array))
    {
        fieldAssociation = vtkDataObject::FIELD_ASSOCIATION_POINTS;
    }
    else if (PeriMultiFitHesseFilterHasArray(input->GetCellData(), array))
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

        this->UpdateProgress(0.5*(inPtId/(float)numPts));
    }

    vtkDebugMacro(<<"Initialize output..");

    vtkSmartPointer<vtkDelaunay3D> delaunay3D =
        vtkSmartPointer<vtkDelaunay3D>::New();
    delaunay3D->SetInputData( input );
    delaunay3D->Update();
    //output =  delaunay3D->GetOutput();

    output->ShallowCopy(delaunay3D->GetOutput());


    // Evaluate fitted function
    vtkIdType numNewPts = output->GetNumberOfPoints();

    vtkDoubleArray* taylor = vtkDoubleArray::New();
    taylor->SetNumberOfComponents(1);
    taylor->SetNumberOfTuples(numNewPts);
    taylor->SetName( "Taylor Function" );

    vtkDoubleArray* hesse = vtkDoubleArray::New();
    hesse->SetNumberOfComponents(9);
    hesse->SetNumberOfTuples(numNewPts);
    hesse->SetName( "Hesse" );

    for( vtkIdType ptId=0; ptId<numNewPts; ptId++) {

        double x[3];
        output->GetPoint(ptId, x); // position of function evaluation

        // get point in original dataset, use as origin
        double center[3],
               coeff[maxCoeff];

        coeffs->GetTuple(ptId, coeff);
        vtkMath::Subtract(x,center,x); // use center as origin

        double a[] = {0,0,0}; // this is the new center

        double f_a;
        double grad_f_a[3];
        double H_f_a[9];

        f_a = evaluateFunction( a, coeff);
        evaluateGradient( a, coeff, grad_f_a );
        evaluateHessian(  a, coeff, H_f_a );
        hesse->SetTuple( ptId, H_f_a);


        // estimate (x-a)
        vtkMath::Subtract(x,a,x); // use center as origin

        double ps[] = { x[0]*x[0], x[1]*x[1], x[2]*x[2],  // precalculate some polygons
                        x[0]*x[1], x[0]*x[2], x[1]*x[2] };

        double val =  f_a 
                    + vtkMath::Dot( x, grad_f_a) // (x-a)^T * grad(f(a))
                    + 1.0/2.0 * // 1/2!
                    (
                           H(0,0)*ps[0]  +    H(1,1)*ps[1]  +    H(2,2)*ps[2]
                      + 2*(H(0,1)*ps[3]) + 2*(H(0,2)*ps[4]) + 2*(H(1,2)*ps[5])
                    ); // (x-a)^T * H_f(a) * (x-a)

        taylor->SetValue( ptId, val);	
        this->UpdateProgress(0.5+0.5*(ptId/(float)numNewPts));

    }

    output->GetPointData()->AddArray(taylor);
    output->GetPointData()->AddArray(hesse);

    return 1;
}

//==============================================================================
vtkSmartPointer<vtkIdList> PeriMultiFitHesse::GetConnectedVertices( 
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



void PeriMultiFitHesse::LeastSquaresFitting(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkDoubleArray> values, double coeff[])
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
            gsl_matrix_set (X, i, hypOff,   pos[0]*pos[1]); //xy
            gsl_matrix_set (X, i, hypOff+1, pos[0]*pos[2]); //xz
            gsl_matrix_set (X, i, hypOff+2, pos[1]*pos[2]); //yz
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

double PeriMultiFitHesse::evaluateFunction(double pos[], double c[])
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

void PeriMultiFitHesse::evaluateGradient(double pos[], double c[], double grad[])
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

void PeriMultiFitHesse::evaluateHessian(double a[], double c[], double H_f_a[])
{
    double x = a[0],
           y = a[1],
           z = a[2];

    size_t conOff =   0;								   // 0 (max)
    size_t linOff =   conOff + FitConstantTerm    * 1;	   // 1
    size_t hypOff =   linOff + FitLinearTerms     * 3;	   // 4
    size_t quaOff =   hypOff + FitHyperbolicTerms * 3;	   // 7
    size_t nCoeff =   quaOff + FitQuadraticTerms  * 3;	   // 10

    H(0,0) = (FitQuadraticTerms) * c[quaOff+0]; 
    H(1,1) = (FitQuadraticTerms) * c[quaOff+1]; 
    H(2,2) = (FitQuadraticTerms) * c[quaOff+2]; 
    H(0,1) = (FitHyperbolicTerms) * c[hypOff+0];
    H(1,0) = H(0,1);							
    H(0,2) = (FitHyperbolicTerms) * c[hypOff+1];
    H(2,0) = H(0,2);						    
    H(1,2) = (FitHyperbolicTerms) * c[hypOff+2];
    H(2,1) = H(0,2);						    
}

//==============================================================================
void PeriMultiFitHesse::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}

int PeriMultiFitHesse::ProcessRequest(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
    return vtkUnstructuredGridAlgorithm::ProcessRequest(request, inputVector, outputVector);
}

int PeriMultiFitHesse::FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
}


