/*=========================================================================

Program:   Visualization Toolkit
Module:    PeriMultiFitTaylor.h

Copyright (c) Michael Bußler

=========================================================================*/

#include "PeriMultiFitTaylor.h"
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
#include "vtkDelaunay3D.h"
#include "vtkCellTreeLocator.h"
#include "vtkGenericCell.h"

#include "linalg.h"

#include <stdio.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>

vtkStandardNewMacro(PeriMultiFitTaylor);

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
bool PeriMultiFitTaylorFilterHasArray(vtkFieldData *fieldData,
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
PeriMultiFitTaylor::PeriMultiFitTaylor()
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

    this->Interpolate = false;

    this->OutputGradient = true;
    this->OutputHessian = true;

    this->FitAlternativeFunction = false;
}
//==============================================================================
PeriMultiFitTaylor::~PeriMultiFitTaylor()
{
}

//-----------------------------------------------------------------------------
void PeriMultiFitTaylor::SetInputScalars(int fieldAssociation, const char *name)
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
void PeriMultiFitTaylor::SetInputScalars(int fieldAssociation,
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
int PeriMultiFitTaylor::RequestInformation ( 
    vtkInformation * vtkNotUsed(request),
    vtkInformationVector ** inputVector,
    vtkInformationVector *outputVector)
{
    //vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation* outInfo = outputVector->GetInformationObject(0);

    //vtkPolyData *input = vtkPolyData::SafeDownCast(
    //    inInfo->Get(vtkDataObject::DATA_OBJECT()));

    // use input bounds for sampling
    //input->GetBounds( this->ModelBounds);

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
int PeriMultiFitTaylor::RequestData(
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
    if (PeriMultiFitTaylorFilterHasArray(input->GetPointData(), array))
    {
        fieldAssociation = vtkDataObject::FIELD_ASSOCIATION_POINTS;
    }
    else if (PeriMultiFitTaylorFilterHasArray(input->GetCellData(), array))
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

        if( this->FitAlternativeFunction )
        {
            size_t nPoints = points->GetNumberOfPoints();
            for( vtkIdType i=0; i<nPoints; i++)
            {
                double val = values->GetValue(i);
                if( val == 0.0 ) // insert midpoint of broken bond
                {
                    values->SetValue(i, 1.0);

                    double newPos[3];
                    points->GetPoint(i, newPos);
                    vtkMath::MultiplyScalar2D(newPos, 0.5);
                    points->InsertNextPoint(newPos);
                    values->InsertNextValue(0.0);
                }
            }
        }

        double coeff[maxCoeff];
        LeastSquaresFitting(points, values, coeff);
        coeffs->SetTuple(inPtId, coeff);

        this->UpdateProgress(0.5*(inPtId/(float)numPts));
    }

    // Evaluate fitted function on rectangular grid
    //
    vtkIdType numNewPts = this->SampleDimensions[0] * this->SampleDimensions[1] * this->SampleDimensions[2];

    vtkDoubleArray* Taylor = vtkDoubleArray::New();
    Taylor->SetNumberOfComponents(1);
    Taylor->SetNumberOfTuples(numNewPts);
    Taylor->SetName( "Taylor Function" );

    vtkDoubleArray* Const = vtkDoubleArray::New();
    Const->SetNumberOfComponents(1);
    Const->SetNumberOfTuples(numNewPts);
    Const->SetName( "Function" );

    vtkDoubleArray* Grad = NULL;
    if( this->OutputGradient) {
        Grad = vtkDoubleArray::New();
        Grad->SetNumberOfComponents(3);
        Grad->SetNumberOfTuples(numNewPts);
        Grad->SetName( "Gradient" );
    }
    vtkDoubleArray* Hesse = NULL;
    if( this->OutputHessian ){
        Hesse = vtkDoubleArray::New();
        Hesse->SetNumberOfComponents(9);
        Hesse->SetNumberOfTuples(numNewPts);
        Hesse->SetName( "Hesse" );
    }

    // init output
    double origin[3], spacing[3];
    ComputeModelBounds(origin, spacing);

    output->SetExtent( 0, this->SampleDimensions[0]-1,
                       0, this->SampleDimensions[1]-1,
                       0, this->SampleDimensions[2]-1 );
    output->SetOrigin(origin);
    output->SetSpacing(spacing);

    if( !this->Interpolate) 
    {

        // evaluate taylor polynom using nearest neighbor point location
        vtkSmartPointer<vtkKdTreePointLocator> tree = vtkSmartPointer<vtkKdTreePointLocator>::New();
        tree->SetDataSet(input);
        tree->BuildLocator();

        for( vtkIdType ptId=0; ptId<numNewPts; ptId++) {

            double x[3];
            output->GetPoint(ptId, x); // position of function evaluation

            // get point in original dataset, use as origin
            double center[3],
                coeff[maxCoeff];
            vtkIdType nn = tree->FindClosestPoint(x);
            input->GetPoint(nn, center);
            coeffs->GetTuple(nn, coeff);
            vtkMath::Subtract(x,center,x); // use center as origin

            double a[] = {0,0,0}; // this is the new center

            double f_a;
            double grad_f_a[3];
            double H_f_a[9];

            // constant term f(a)
            f_a = evaluateFunction( a, coeff);

            // gradient(f(a))
            evaluateGradient( a, coeff, grad_f_a );

            // Hessian(f(a))
            evaluateHessian(  a, coeff, H_f_a );

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

            Taylor->SetValue( ptId, val);

            // evaluate fitted function at position x
            f_a = evaluateFunction( x, coeff);
            Const->SetValue( ptId, f_a );
            if( this->OutputGradient ) {
                evaluateGradient( x, coeff, grad_f_a );
                Grad->SetTuple( ptId, grad_f_a );
            }
            if( this->OutputHessian) {
                evaluateHessian(  x, coeff, H_f_a );
                Hesse->SetTuple( ptId, H_f_a);
            }

            this->UpdateProgress(0.5+0.5*(ptId/(float)numNewPts));
        }
    } else {

        // create delaunay triangulation
        vtkSmartPointer<vtkDelaunay3D> delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
        delaunay->SetInputData(input);
        delaunay->Update();
        vtkSmartPointer<vtkUnstructuredGrid> grid = delaunay->GetOutput();


        vtkSmartPointer<vtkCellTreeLocator> loc = vtkSmartPointer<vtkCellTreeLocator>::New();
        loc->SetDataSet(grid);
        loc->AutomaticOn();
        loc->BuildLocator();

        for( vtkIdType ptId=0; ptId<numNewPts; ptId++) {

            double x[3];
            output->GetPoint(ptId, x); // position of function evaluation

            double pcoords[3];
            double weights[4];
            vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
            vtkIdType cellId = loc->FindCell(x, 0.0, cell, pcoords, weights);

            double total_taylor = 0;
            double total_f_a = 0;
            double total_grad_f_a[3] = {0,0,0};
            double total_H_f_a[9]    = {0,0,0, 0,0,0, 0,0,0};

            if( cellId != -1 )
            {

                for( int i = 0; i<cell->GetNumberOfPoints(); i++)
                {

                    // evaluate taylor polynom of cell points at x
                    vtkIdType cellPtId = cell->GetPointId(i);

                    // get point in original dataset, use as origin
                    double center[3],
                            coeff[maxCoeff];

                    input->GetPoint(cellPtId, center);
                    coeffs->GetTuple(cellPtId, coeff);
                    vtkMath::Subtract(x,center,x); // use center as origin

                    double a[] = {0,0,0}; // this is the new center
                    double f_a;
                    double grad_f_a[3];
                    double H_f_a[9];

                    // constant term f(a)
                    f_a = evaluateFunction( a, coeff);

                    // gradient(f(a))
                    evaluateGradient( a, coeff, grad_f_a );

                    // Hessian(f(a))
                    evaluateHessian(  a, coeff, H_f_a );

                    // estimate (x-a)
                    vtkMath::Subtract(x,a,x); // use center as origin

                    double ps[] = { x[0]*x[0], x[1]*x[1], x[2]*x[2],  // precalculate some polygons
                        x[0]*x[1], x[0]*x[2], x[1]*x[2] };

                    double val = f_a 
                        + vtkMath::Dot( x, grad_f_a) // (x-a)^T * grad(f(a))
                        + 1.0/2.0 * // 1/2!
                        (
                                H(0,0)*ps[0]  +    H(1,1)*ps[1]  +    H(2,2)*ps[2]
                        + 2*(H(0,1)*ps[3]) + 2*(H(0,2)*ps[4]) + 2*(H(1,2)*ps[5])
                        ); // (x-a)^T * H_f(a) * (x-a)

                    // store interpolated results
                    total_taylor += weights[i]*val;
                    total_f_a    += weights[i]*f_a;

                    for( int j=0; j<3; j++)
                        total_grad_f_a[j] += weights[i]*grad_f_a[j];

                    for( int j=0; j<9; j++)
                        total_H_f_a[j] += weights[i]*H_f_a[j];

                }// end for cell pts
            } // end if cell

            Const->SetValue( ptId, total_f_a );
            Taylor->SetValue( ptId, total_taylor);	

            if( this->OutputGradient )
                Grad->SetTuple( ptId, total_grad_f_a );

            if( this->OutputHessian)
                Hesse->SetTuple( ptId, total_H_f_a);

            this->UpdateProgress(0.5+0.5*(ptId/(float)numNewPts));
        }
    }

    output->GetPointData()->AddArray(Const);
    if( this->OutputGradient )
        output->GetPointData()->AddArray(Grad);
    if( this->OutputHessian )
        output->GetPointData()->AddArray(Hesse);
    output->GetPointData()->AddArray(Taylor);

    return 1;
}

//==============================================================================
vtkSmartPointer<vtkIdList> PeriMultiFitTaylor::GetConnectedVertices( 
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


void PeriMultiFitTaylor::LeastSquaresFitting(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkDoubleArray> values, double coeff[])
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

    int status = 0;
    gsl_set_error_handler_off(); // hey gsl, please, don't crash all the time!!

    gsl_multifit_linear_workspace * work 
        = gsl_multifit_linear_alloc (nPoints, nCoeff);
    status = 
        gsl_multifit_wlinear (X, w, v, c, cov, &chisq, work);

    if( status ) // check if sth went wrong
    {
        for (size_t i=0; i<nCoeff; i++)
            coeff[i] = 0.0;
    } else {
        gsl_multifit_linear_free (work);

        for (size_t i=0; i<nCoeff; i++)
            coeff[i] = gsl_vector_get(c, i);
    }

    gsl_matrix_free (X);
    gsl_vector_free (v);
    gsl_vector_free (w);
    gsl_vector_free (c);
    gsl_matrix_free (cov);
}

double PeriMultiFitTaylor::evaluateFunction(double pos[], double c[])
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

void PeriMultiFitTaylor::evaluateGradient(double pos[], double c[], double grad[])
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

void PeriMultiFitTaylor::evaluateHessian(double a[], double c[], double H_f_a[])
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
void PeriMultiFitTaylor::PrintSelf(ostream &os, vtkIndent indent)
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

int PeriMultiFitTaylor::ProcessRequest(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
    return vtkImageAlgorithm::ProcessRequest(request, inputVector, outputVector);
}

int PeriMultiFitTaylor::FillInputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
}

// Set the i-j-k dimensions on which to sample the distance function.
void PeriMultiFitTaylor::SetSampleDimensions(int i, int j, int k)
{
    int dim[3];

    dim[0] = i;
    dim[1] = j;
    dim[2] = k;

    this->SetSampleDimensions(dim);
}

// Set the i-j-k dimensions on which to sample the distance function.
void PeriMultiFitTaylor::SetSampleDimensions(int dim[3])
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

// Compute ModelBounds from input geometry.
void PeriMultiFitTaylor::ComputeModelBounds(double origin[3], double spacing[3])
{

    // compute model bounds if not set previously
    vtkDataSet *ds = vtkDataSet::SafeDownCast(this->GetInput());
    ds->GetBounds( this->ModelBounds);

    // Set volume origin and data spacing
    for (int i=0; i<3; i++)
    {
        origin[i] = this->ModelBounds[2*i];
        spacing[i] = (this->ModelBounds[2*i+1] - this->ModelBounds[2*i])
            / (this->SampleDimensions[i] - 1);
    }
}
