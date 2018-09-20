#include "vtkTetTrimSurface.h"

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
#include "vtkDataSetSurfaceFilter.h"

#include <math.h>
#include <algorithm>

vtkStandardNewMacro(vtkTetTrimSurface);

//----------------------------------------------------------------------------
vtkTetTrimSurface::vtkTetTrimSurface()
{
    this->TrimDepth = 1;
}

//----------------------------------------------------------------------------
vtkTetTrimSurface::~vtkTetTrimSurface()
{
}

//----------------------------------------------------------------------------
//
// Clip through data generating surface.
//
int vtkTetTrimSurface::RequestData(
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

    vtkDebugMacro(<< "Trimming Tetrahedral Surface");

    // Initialize self; create output objects
    //
    if ( !input || numPts < 1 || inPts == NULL )
    {
        vtkDebugMacro(<<"No data.");
        return 1;
    }

    output->ShallowCopy(input);

    for( int i=0; i<TrimDepth; i++) {

        // get surface with respective cell ids
        vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = 
            vtkSmartPointer<vtkDataSetSurfaceFilter>::New();

        surfaceFilter->SetInputData(output);
        surfaceFilter->PassThroughCellIdsOn();
        surfaceFilter->Update(); 

        vtkPolyData* surface = surfaceFilter->GetOutput();
        vtkIdTypeArray* surfaceCells = vtkIdTypeArray::SafeDownCast(
            surface->GetCellData()->GetArray( surfaceFilter->GetOriginalCellIdsName()));

        vtkSmartPointer<vtkCellArray> tets = vtkSmartPointer<vtkCellArray>::New();

        // iterate over all cells, insert cells that do not belong to surface
        for( cellId=0; cellId<output->GetNumberOfCells(); cellId++) 
        {
            if( surfaceCells->LookupValue(cellId) < 0){
                tets->InsertNextCell(output->GetCell(cellId));
            }
        }

        output->SetCells(VTK_TETRA, tets);
    }

    output->Squeeze();

    return 1;
}
