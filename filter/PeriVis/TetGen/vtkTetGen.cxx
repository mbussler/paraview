#include "vtkTetGen.h"

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

#include <math.h>
#include <algorithm>
#include "vtkUnstructuredGrid.h"
#include "vtkDataSet.h"
#include "vtkIdList.h"

#include "tetgen.h"

/* performance measure */
#include "timer.h"
#include <QElapsedTimer>

vtkStandardNewMacro(vtkTetGen);

//----------------------------------------------------------------------------
vtkTetGen::vtkTetGen()
{
    this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
vtkTetGen::~vtkTetGen()
{
}

//----------------------------------------------------------------------------
//
// Clip through data generating surface.
//
int vtkTetGen::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
    QElapsedTimer timer;

    // get the info objects
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and output
    vtkDataSet *input = vtkDataSet::SafeDownCast(
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
    int numberOfPoints;
    vtkPointData *inPD=input->GetPointData(), *outPD = output->GetPointData();
    vtkCellData *inCD=input->GetCellData(), *outCD = output->GetCellData();
    vtkCellData *outClippedCD = NULL;

    vtkDebugMacro(<< "Calculating Delaunay Triangulation");

    // Initialize self; create output objects
    //
    if ( !input || numPts < 1 || input->GetNumberOfPoints() == 0 )
    {
        vtkDebugMacro(<<"No data.");
        return 1;
    }

    tetgenio tetgen_in;
    tetgen_in.numberofpoints = numPts;
    tetgen_in.pointlist = new REAL[3*numPts];

    for( vtkIdType ptId = 0; ptId<numPts; ptId++)
    {
        input->GetPoint(ptId, &tetgen_in.pointlist[3*ptId]);
    }
    
    tetgenio tetgen_out;

    timer.start();
    tetrahedralize("", &tetgen_in, &tetgen_out);
    write_timer("tetgen", "tetrahedralize", timer.elapsed());


    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for( vtkIdType ptId = 0; ptId<tetgen_out.numberofpoints; ptId++)
    {
        points->InsertNextPoint(&tetgen_out.pointlist[3*ptId]);
    }


    output->SetPoints(points);

    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    for( vtkIdType cellId = 0; cellId < tetgen_out.numberoftetrahedra; cellId++)
    {
        vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
        ids->InsertNextId(tetgen_out.tetrahedronlist[4*cellId+0]);
        ids->InsertNextId(tetgen_out.tetrahedronlist[4*cellId+1]);
        ids->InsertNextId(tetgen_out.tetrahedronlist[4*cellId+2]);
        ids->InsertNextId(tetgen_out.tetrahedronlist[4*cellId+3]);
        cells->InsertNextCell(ids);
    }
    output->SetCells(VTK_TETRA, cells);
    
    // copy point data arrays
    int numArrays = inPD->GetNumberOfArrays();
    for( int i=0; i<numArrays; i++ ){
        vtkDataArray* arr = vtkDataArray::CreateDataArray( inPD->GetArray(i)->GetDataType());
        arr->SetNumberOfComponents( inPD->GetArray(i)->GetNumberOfComponents());
        arr->SetName( inPD->GetArray(i)->GetName());
        arr->DeepCopy( inPD->GetArray(i));
        outPD->AddArray(arr);
        arr->Delete();
    }

    output->Squeeze();
    return 1;
}

int vtkTetGen::FillInputPortInformation(int port, vtkInformation* info)
{
    if ( port == 0 ) {
        info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet" );
        return 1;
    }
    return 0;
}

int vtkTetGen::FillOutputPortInformation(int port, vtkInformation* info)
{
    if ( port == 0 ) {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
        return 1;
    }
    return 0;
}
