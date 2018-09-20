#include "vtkGaussianSmooth.h"

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

vtkStandardNewMacro(vtkGaussianSmooth);

//----------------------------------------------------------------------------
vtkGaussianSmooth::vtkGaussianSmooth()
{
    this->SetNumberOfOutputPorts(1);

    // by default process active point scalars
    this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
       vtkDataSetAttributes::SCALARS);
}

//----------------------------------------------------------------------------
vtkGaussianSmooth::~vtkGaussianSmooth()
{
}

//----------------------------------------------------------------------------
//
// Clip through data generating surface.
//
int vtkGaussianSmooth::RequestData(
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

    vtkDebugMacro(<< "Calculating Gaussian Smoothing");

    // Initialize self; create output objects
    //
    if ( !input || numPts < 1 || inPts == NULL )
    {
        vtkDebugMacro(<<"No data.");
        return 1;
    }

    output->ShallowCopy(input);

    // allocate mem for output data arrays
    int type = scalarArray->GetDataType();
    std::string name = std::string(scalarArray->GetName())+"_smooth";
    vtkDataArray* smoothed = vtkDataArray::CreateDataArray(type);
    smoothed->SetName( name.c_str() );
    smoothed->SetNumberOfComponents(scalarArray->GetNumberOfComponents());
    smoothed->SetNumberOfTuples(scalarArray->GetNumberOfTuples());
    smoothed->DeepCopy(scalarArray); // init with scalar data

    vtkSmartPointer<vtkKdTree> kd = vtkSmartPointer<vtkKdTree>::New();
    kd->BuildLocatorFromPoints(inPts);

    vtkSmartPointer<vtkIdList> neighPoints = vtkSmartPointer<vtkIdList>::New();
    double center[3], pos[3];

    double sig2 = SmoothRadius*SmoothRadius;
    double sig3 = sig2 * SmoothRadius;
    double f    = sqrt(8.0*M_PI*M_PI*M_PI);

    for (vtkIdType point = 0; point < numPts; point++)
    {
        input->GetPoint(point, center);
        kd->FindPointsWithinRadius(SmoothRadius, center, neighPoints);

        double weight = 0, totalweight = 0;
        double distance2;

        int numNodes = neighPoints->GetNumberOfIds();
        double* weights = new double[numNodes];

        // estimate weights of neighbors
        for( vtkIdType neighPtId = 0; neighPtId < numNodes; neighPtId++ )
        {
            vtkIdType ptId = neighPoints->GetId(neighPtId);
            input->GetPoint(ptId, pos);
            distance2 = vtkMath::Distance2BetweenPoints(pos, center);
            weight = (1.0/(sig3 * f))*exp(-distance2/2.*sig2);
            weights[neighPtId] = weight;
            totalweight += weight;
        }        

        double value = 0;
        for( vtkIdType neighPtId = 0; neighPtId < numNodes; neighPtId++ )
        {
            vtkIdType ptId = neighPoints->GetId(neighPtId);
            value += scalarArray->GetTuple1(ptId) * weights[neighPtId] / totalweight;
        }
        smoothed->SetTuple1(point, value);

        delete[] weights;
    }

    output->GetPointData()->AddArray(smoothed);
    smoothed->Delete();

    output->Squeeze();
//
    return 1;
}
