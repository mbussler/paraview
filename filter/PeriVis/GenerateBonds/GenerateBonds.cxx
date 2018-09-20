#include "GenerateBonds.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkGenericCell.h"
#include "vtkImplicitFunction.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkLine.h"
#include "vtkMergePoints.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkTriangle.h"
#include "vtkIncrementalPointLocator.h"
#include <vtkKdTreePointLocator.h>

#include <math.h>
#include "vtkSmartPointer.h"
#include <algorithm>

vtkStandardNewMacro(GenerateBonds);

//----------------------------------------------------------------------------
GenerateBonds::GenerateBonds()
{
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
GenerateBonds::~GenerateBonds()
{
}

//----------------------------------------------------------------------------
//
// Clip through data generating surface.
//
int GenerateBonds::RequestData(
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
    vtkPolyData *output = vtkPolyData::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkIdType cellId, i, updateTime;
    vtkPoints *cellPts;
    vtkDataArray *clipScalars;
    vtkFloatArray *cellScalars;
    vtkCellArray *newVerts, *newLines, *newPolys, *connList=NULL;
    vtkCellArray *clippedVerts=NULL, *clippedLines=NULL;
    vtkCellArray *clippedPolys=NULL, *clippedList=NULL;
    vtkPoints *newPoints;
    vtkIdList *cellIds;
    double s;
    vtkIdType estimatedSize, numCells=input->GetNumberOfCells();
    vtkIdType numPts=input->GetNumberOfPoints();
    vtkPoints *inPts=input->GetPoints();
    int numberOfPoints;
    vtkPointData *inPD=input->GetPointData(), *outPD = output->GetPointData();
    vtkCellData *inCD=input->GetCellData(), *outCD = output->GetCellData();
    vtkCellData *outClippedCD = NULL;

    vtkDebugMacro(<< "Clipping polygonal data");

    // Initialize self; create output objects
    //
    if ( numPts < 1 || inPts == NULL )
    {
        vtkDebugMacro(<<"No data to clip");
        return 1;
    }

    // Generate bonds
    vtkSmartPointer<vtkCellArray> bonds = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkKdTreePointLocator> kd = vtkSmartPointer<vtkKdTreePointLocator>::New();

    kd->SetDataSet(input);
    kd->BuildLocator();

    vtkSmartPointer<vtkIdList> nnPoints = vtkSmartPointer<vtkIdList>::New();

    for( vtkIdType ptId=0; ptId<numPts; ptId++) 
    {
        double pos[3];
        inPts->GetPoint(ptId, pos);

        kd->FindPointsWithinRadius( this->NeighborhoodSize, pos, nnPoints);

        for( int i=0; i<nnPoints->GetNumberOfIds(); i++)
        {
            vtkIdType nnId = nnPoints->GetId(i);
            bonds->InsertNextCell(2);
            bonds->InsertCellPoint(ptId);
            bonds->InsertCellPoint(nnId);
        }

        this->UpdateProgress(ptId/numPts);
    }

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->DeepCopy(input->GetPoints());
    output->SetPoints(points);
    output->SetLines(bonds);

    output->Squeeze();

    return 1;
}


//----------------------------------------------------------------------------
void GenerateBonds::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}
