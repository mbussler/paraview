#include "vtkTetSmooth.h"

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

vtkStandardNewMacro(vtkTetSmooth);

//----------------------------------------------------------------------------
vtkTetSmooth::vtkTetSmooth()
{
    this->SetNumberOfOutputPorts(1);

    // by default process active point scalars
    this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
       vtkDataSetAttributes::SCALARS);

    this->SmoothingMethod = 0;
    this->CenterWeight = 1;
    this->NumberOfIterations = 1;
}

//----------------------------------------------------------------------------
vtkTetSmooth::~vtkTetSmooth()
{
}

//----------------------------------------------------------------------------
//
// Clip through data generating surface.
//
int vtkTetSmooth::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and output
    vtkUnstructuredGrid *input = vtkUnstructuredGrid::SafeDownCast(
        inInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
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

    vtkDebugMacro(<< "Calculating Tetrahedral Smoothing");

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
    vtkSmartPointer<vtkDataArray> smoothed = vtkDataArray::CreateDataArray(type);
    smoothed->SetName( name.c_str() );
    smoothed->SetNumberOfComponents(scalarArray->GetNumberOfComponents());
    smoothed->SetNumberOfTuples(scalarArray->GetNumberOfTuples());
    smoothed->DeepCopy(scalarArray); // init with scalar data
    output->GetPointData()->AddArray(smoothed);


    if( this->SmoothingMethod == 0) // stencil smoothing
    {
        vtkSmartPointer<vtkDataArray> proc = vtkDataArray::CreateDataArray(type);
        proc->SetName( "proc" );
        proc->SetNumberOfComponents(scalarArray->GetNumberOfComponents());
        proc->SetNumberOfTuples(scalarArray->GetNumberOfTuples());

        vtkSmartPointer<vtkIdList> currentPoint = vtkSmartPointer<vtkIdList>::New();
        currentPoint->SetNumberOfIds(1);
        vtkSmartPointer<vtkIdList> cellsOnPoint = vtkSmartPointer<vtkIdList>::New();

        for( int i = 0; i<NumberOfIterations; i++) 
        {
            for (vtkIdType point = 0; point < numPts; point++)
            {
                currentPoint->SetId(0, point);
                // Get all cells touching this point.
                input->GetCellNeighbors(-1, currentPoint, cellsOnPoint);
                vtkIdType numCellNeighbors = cellsOnPoint->GetNumberOfIds();
                vtkIdType numValidCellNeighbors = 0;

                double value = smoothed->GetTuple1(point)*this->CenterWeight;

                // Iterate on all cells and find all points connected to current point
                // by an edge.
                for (vtkIdType neighbor = 0; neighbor < numCellNeighbors; neighbor++)
                {
                    vtkCell *cell = input->GetCell(cellsOnPoint->GetId(neighbor));

                    int NumberOfCellPoints = cell->GetNumberOfPoints();
                    for (int i = 0; i < NumberOfCellPoints; i++)
                    {
                        vtkIdType ptId = cell->GetPointId(i);
                        if( ptId != point ){
                            numValidCellNeighbors++;
                            value += smoothed->GetTuple1(ptId);
                        }
                    }
                } // iterating over neighbors

                value /= (numValidCellNeighbors+this->CenterWeight);

                proc->SetTuple1(point, value);
            }

            smoothed->DeepCopy(proc); // write back processed values
        }
        proc->Delete();
    }
    else if( this->SmoothingMethod == 2) // directional smoothing
    {
        vtkSmartPointer<vtkKdTree> kd = vtkSmartPointer<vtkKdTree>::New();
        kd->BuildLocatorFromPoints(inPts);

        vtkSmartPointer<vtkDoubleArray> covArray = vtkSmartPointer<vtkDoubleArray>::New();
        covArray->SetName( "Covariance" );
        covArray->SetNumberOfComponents(9);
        covArray->SetNumberOfTuples(numPts);
        output->GetPointData()->AddArray(covArray);       

        vtkSmartPointer<vtkIdList> neighPoints = vtkSmartPointer<vtkIdList>::New();
        double center[3], pos[3];
        double cov[9];

        for (vtkIdType point = 0; point < numPts; point++)
        {
            input->GetPoint(point, center);
            kd->FindPointsWithinRadius(SmoothRadius, center, neighPoints);

            // build (weighted) covariance matrix from neighboring points
            vtkIdType numNodes = neighPoints->GetNumberOfIds();
            double* pointSet[3];
            pointSet[0] = new double[numNodes];
            pointSet[1] = new double[numNodes];
            pointSet[2] = new double[numNodes];

            for( vtkIdType neighPtId = 0; neighPtId < numNodes; neighPtId++ )
            {
                vtkIdType ptId = neighPoints->GetId(neighPtId);
                input->GetPoint(ptId, pos);
                double weight = scalarArray->GetTuple1(ptId);
                pointSet[0][neighPtId] =  (pos[0]-center[0]) * weight;
                pointSet[1][neighPtId] =  (pos[1]-center[1]) * weight;
                pointSet[2][neighPtId] =  (pos[2]-center[2]) * weight;
            }        

            //Covariance Matrix
            //mat3 m;
            for(int j=0; j<3; ++j) {
                for(int i=0; i<3; ++i) {
                    double sum = 0.0;
                    for(int n=0; n<numNodes; ++n){
                        sum += pointSet[i][n] * pointSet[j][n];
                    }
                    //m[i][j] = sum/(double)numNodes;
                    cov[3*j+i] = sum/(double)numNodes;
                }
            }
            covArray->SetTuple(point, cov);

            delete[] pointSet[0];
            delete[] pointSet[1];
            delete[] pointSet[2];
        }
    }
    else if( this->SmoothingMethod == 1) // distance weighted
    {
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
    }
    
    smoothed->Delete();
    output->Squeeze();
//
    return 1;
}

int vtkTetSmooth::FillInputPortInformation(int port, vtkInformation* info)
{
    if ( port == 0 ) {
        info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid" );
        return 1;
    }
    return 0;
}

int vtkTetSmooth::FillOutputPortInformation(int port, vtkInformation* info)
{
    if ( port == 0 ) {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
        return 1;
    }
    return 0;
}
