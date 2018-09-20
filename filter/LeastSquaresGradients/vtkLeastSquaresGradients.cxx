#include "vtkLeastSquaresGradients.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
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
#include "vtkKdTree.h"
#include "vtkDoubleArray.h"
#include "vtkUnstructuredGrid.h"

#include <math.h>
#include <algorithm>

#include "linalg.h"

/* performance measure */
#include "timer.h"
#include <QElapsedTimer>

vtkStandardNewMacro(vtkLeastSquaresGradients);

//----------------------------------------------------------------------------
vtkLeastSquaresGradients::vtkLeastSquaresGradients()
{
    this->SetNumberOfOutputPorts(1);

    // by default process active point scalars
    this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
       vtkDataSetAttributes::SCALARS);

    this->Radius = 1.0;
}

//----------------------------------------------------------------------------
vtkLeastSquaresGradients::~vtkLeastSquaresGradients()
{
}

//----------------------------------------------------------------------------
//
// Clip through data generating surface.
//
int vtkLeastSquaresGradients::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and output
    vtkPointSet *input = vtkPointSet::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPointSet *output = vtkPointSet::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkIdType cellId, i, updateTime;
    vtkPoints *cellPts;
    vtkDataArray *clipScalars;
    vtkFloatArray *cellScalars;
    vtkCellArray *newTris, *newLines, *newPolys, *connList=NULL;
    vtkPoints *newPoints;
    double s;
    vtkIdType estimatedSize, numCells=input->GetNumberOfCells();
    vtkIdType numPts=input->GetNumberOfPoints();
    vtkPoints *inPts=input->GetPoints();
    int numberOfPoints;
    vtkPointData *inPD=input->GetPointData(), *outPD = output->GetPointData();
    vtkCellData *inCD=input->GetCellData(), *outCD = output->GetCellData();
    vtkCellData *outClippedCD = NULL;

    vtkDataArray *scalarArray = this->GetInputArrayToProcess(0, inputVector);
    if( !scalarArray )
    {
        vtkErrorMacro("No input array selected");
        return 0;
    }

    vtkDebugMacro(<< "Calculating Least Squares Gradients");

    // Initialize self; create output objects
    //
    if ( !input || numPts < 1 || inPts == NULL )
    {
        vtkDebugMacro(<<"No data.");
        return 1;
    }

    output->ShallowCopy(input);

    // allocate mem for output data arrays
    vtkSmartPointer<vtkDoubleArray> gradient = vtkSmartPointer<vtkDoubleArray>::New();
    gradient->SetName( "Gradient" );
    gradient->SetNumberOfComponents(3);
    gradient->SetNumberOfTuples(numPts);
    output->GetPointData()->AddArray(gradient);

    vtkSmartPointer<vtkDoubleArray> hesse = vtkSmartPointer<vtkDoubleArray>::New();
    hesse->SetName( "Hesse" );
    hesse->SetNumberOfComponents(9);
    hesse->SetNumberOfTuples(numPts);
    output->GetPointData()->AddArray(hesse);

    vtkSmartPointer<vtkKdTree> kd = vtkSmartPointer<vtkKdTree>::New();
    kd->BuildLocatorFromPoints(inPts);

    std::vector<vtkSmartPointer<vtkIdList> > neighbors;

    double center[3], pos[3];

    QElapsedTimer timer;
    timer.start();

    // calculate gradient in a least squares fashion
    for (vtkIdType point = 0; point < numPts; point++)
    {
        input->GetPoint(point, center);
        vtkSmartPointer<vtkIdList> neighPoints = vtkSmartPointer<vtkIdList>::New();
        kd->FindPointsWithinRadius(Radius, center, neighPoints);
        neighbors.push_back(neighPoints);

        mat3 m = {{0}, {0}, {0}};
        vec3 rhs = {0};
        double fi = scalarArray->GetTuple1(point);
        int numNodes = neighPoints->GetNumberOfIds();

        // estimate weights of neighbors
        for( vtkIdType neighPtId = 0; neighPtId < numNodes; neighPtId++ )
        {
            vtkIdType ptId = neighPoints->GetId(neighPtId);
            input->GetPoint(ptId, pos);
            vtkMath::Subtract(pos, center, pos);

            m[0][0] += pos[0]*pos[0]; m[0][1] += pos[0]*pos[1]; m[0][2] += pos[0]*pos[2];
            m[1][0] += pos[1]*pos[0]; m[1][1] += pos[1]*pos[1]; m[1][2] += pos[1]*pos[2];
            m[2][0] += pos[2]*pos[0]; m[2][1] += pos[2]*pos[1]; m[2][2] += pos[2]*pos[2];
            
            double f = scalarArray->GetTuple1(ptId) - fi;
            rhs[0] += f*pos[0]; 
            rhs[1] += f*pos[1]; 
            rhs[2] += f*pos[2];
        }
        double det = mat3det(m);
        if( det == 0) det = 1;
        vec3 grad;

        grad[0] = vec3det( rhs, m[1], m[2]) / det;	// dV0/dx
        grad[1] = vec3det( m[0], rhs, m[2]) / det;	// dV0/dy
        grad[2] = vec3det( m[0], m[1], rhs) / det;	// dV0/dz

        for (int v=0; v<3; v++) {
            if (grad[v] > FLT_MAX) {
                grad[v] = FLT_MAX;
            }
            else if (grad[v] < -FLT_MAX) {
                grad[v] = -FLT_MAX;
            }
        }

        gradient->SetTuple(point, grad);
    }

    write_timer("LeastSquaresGradient", "gradient", timer.elapsed());

    timer.restart();

    // calculate hesse in a least squares fashion
    for (vtkIdType point = 0; point < numPts; point++)
    {
        input->GetPoint(point, center);
        vtkSmartPointer<vtkIdList> neighPoints = neighbors[point];

        mat3 m = {{0}, {0}, {0}};
        vec3 rhsX = { 0 };
        vec3 rhsY = { 0 };
        vec3 rhsZ = { 0 };
        vec3 gradX;
        gradient->GetTuple(point, gradX);
        int numNodes = neighPoints->GetNumberOfIds();

        // estimate weights of neighbors
        for( vtkIdType neighPtId = 0; neighPtId < numNodes; neighPtId++ )
        {
            vtkIdType ptId = neighPoints->GetId(neighPtId);
            input->GetPoint(ptId, pos);
            vtkMath::Subtract(pos, center, pos);

            m[0][0] += pos[0]*pos[0]; m[0][1] += pos[0]*pos[1]; m[0][2] += pos[0]*pos[2];
            m[1][0] += pos[1]*pos[0]; m[1][1] += pos[1]*pos[1]; m[1][2] += pos[1]*pos[2];
            m[2][0] += pos[2]*pos[0]; m[2][1] += pos[2]*pos[1]; m[2][2] += pos[2]*pos[2];
            
            // Relative function values of neighbors
            vec3 grad;
            gradient->GetTuple(ptId, grad);
            double gradientX = grad[0] - gradX[0];
            double gradientY = grad[1] - gradX[1];
            double gradientZ = grad[2] - gradX[2];
            
            rhsX[0] += gradientX*pos[0]; rhsX[1] += gradientX*pos[1]; rhsX[2] += gradientX*pos[2];
            rhsY[0] += gradientY*pos[0]; rhsY[1] += gradientY*pos[1]; rhsY[2] += gradientY*pos[2];
            rhsZ[0] += gradientZ*pos[0]; rhsZ[1] += gradientZ*pos[1]; rhsZ[2] += gradientZ*pos[2];
        }
        double det = mat3det(m);
        if( det == 0) det = 1;
        double h[9];

        h[0] = vec3det(rhsX, m[1], m[2]) / det;	// dV0/dx
        h[1] = vec3det(m[0], rhsX, m[2]) / det;	// dV0/dy
        h[2] = vec3det(m[0], m[1], rhsX) / det;	// dV0/dz

        h[3] = vec3det(rhsY, m[1], m[2]) / det;	// dV0/dx
        h[4] = vec3det(m[0], rhsY, m[2]) / det;	// dV0/dy
        h[5] = vec3det(m[0], m[1], rhsY) / det;	// dV0/dz

        h[6] = vec3det(rhsZ, m[1], m[2]) / det;	// dV0/dx
        h[7] = vec3det(m[0], rhsZ, m[2]) / det;	// dV0/dy
        h[8] = vec3det(m[0], m[1], rhsZ) / det;	// dV0/dz

        hesse->SetTuple(point, h);
    }

    write_timer("LeastSquaresGradient", "hessian", timer.elapsed());

    output->Squeeze();
//
    return 1;
}
