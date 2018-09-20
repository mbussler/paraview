#include "PolyPointSet.h"

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
#include "vtkUnstructuredGrid.h"
#include "vtkTriangle.h"
#include "vtkIncrementalPointLocator.h"
#include <vtkImplicitPolyDataDistance.h>
#include <vtkPolygon.h>

#include <math.h>
#include "vtkSmartPointer.h"
#include <algorithm>

vtkStandardNewMacro(PolyPointSet);

//----------------------------------------------------------------------------
PolyPointSet::PolyPointSet(vtkImplicitFunction *cf)
{
    this->MaxClipDistance = 0.0;
    this->MinClipDistance = 0.0;

	this->SetNumberOfInputPorts(2);
	this->SetNumberOfOutputPorts(1);
}

int PolyPointSet::FillInputPortInformation(int port, vtkInformation* info)
{
    if ( port == 0 ) {
        info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid" );
        return 1;
    }
    if ( port == 1 ) {
        info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData" );
        return 1;
    }
    return 0;

}

//----------------------------------------------------------------------------
PolyPointSet::~PolyPointSet()
{
}

//----------------------------------------------------------------------------
//
// Clip through data generating surface.
//
int PolyPointSet::RequestData(
	vtkInformation *vtkNotUsed(request),
	vtkInformationVector **inputVector,
	vtkInformationVector *outputVector)
{
	// get the info objects
	vtkInformation *inDataInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *inClipPolyInfo = inputVector[1]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the input and output
	vtkUnstructuredGrid *input = vtkUnstructuredGrid::SafeDownCast(
		inDataInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *inClipPD = vtkPolyData::SafeDownCast(
		inClipPolyInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
		outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkIdType cellId, i, updateTime;
	vtkPoints *cellPts;
	vtkDataArray *clipScalars=NULL;
	vtkFloatArray *cellScalars;
	vtkCellArray *newVerts=NULL, *newLines=NULL, *newPolys=NULL, *connList=NULL;
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

	vtkDebugMacro(<< "Clipping unstructured grid data");

	// Initialize self; create output objects
	//
	if ( numPts < 1 || inPts == NULL )
	{
		vtkDebugMacro(<<"No data to clip");
		return 1;
	}

    if( inClipPD->GetNumberOfCells() == 0)
    {
        vtkDebugMacro(<<"Empty clip polygon");
        return 1;
    }

	// generate implicit function from poly data 
	vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance =
		vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
	implicitPolyDataDistance->SetInput(inClipPD);

	// copy point data arrays
	int numArrays = inPD->GetNumberOfArrays();
	for( int i=0; i<numArrays; i++ ){
		vtkDataArray *arr = vtkDataArray::CreateDataArray(inPD->GetArray(i)->GetDataType());
		arr->SetNumberOfComponents(inPD->GetArray(i)->GetNumberOfComponents());
		arr->SetName(inPD->GetArray(i)->GetName());
		outPD->AddArray(arr);
        arr->Delete();
	}
    if( outPD->HasArray("Normals"))
        outPD->SetActiveNormals("Normals");
    if( outPD->HasArray("TextureCoordinates"))
        outPD->SetActiveTCoords("TextureCoordinates");

	// clip point data
	newPoints =  vtkPoints::New();
	//newVerts = vtkCellArray::New();
	std::vector<vtkIdType> pointIds;
    pointIds.resize(numPts, -1);

	double pos[3];
	double data[9];
	vtkIdType vertexId = 0;

	for( vtkIdType ptId=0; ptId<numPts; ptId++) {

		inPts->GetPoint(ptId, pos);
		double distance = abs(implicitPolyDataDistance->EvaluateFunction(pos));

		if( distance <= this->MaxClipDistance &&
            distance >= this->MinClipDistance )
		{
			vtkIdType newPtId = newPoints->InsertNextPoint(pos);

			// copy point data
			for( int i=0; i<numArrays; i++ ){
				inPD->GetArray(i)->GetTuple(ptId, data);
				outPD->GetArray(i)->InsertNextTuple(data);
			}

			pointIds[ptId] = newPtId;
		}

		this->UpdateProgress(0.5*ptId/numPts);
	}
	output->SetPoints(newPoints);
	newPoints->Delete();

 	// copy cell data arrays
    numArrays = inCD->GetNumberOfArrays();
	for( int i=0; i<numArrays; i++ ){
		vtkDataArray *arr = vtkDataArray::CreateDataArray(inCD->GetArray(i)->GetDataType());;
		arr->SetNumberOfComponents(inCD->GetArray(i)->GetNumberOfComponents());
		arr->SetName(inCD->GetArray(i)->GetName());
        outCD->AddArray(arr);
        arr->Delete();
	}

    cellId = 0;

    //// Process Tetrahedra
    vtkCellArray* newTets = vtkCellArray::New();
    vtkSmartPointer<vtkIdList> newCellPtIds;

    // iterate over all cells, insert cells that do not belong to surface
    for( cellId=0; cellId<input->GetNumberOfCells(); cellId++) 
    {
        vtkCell* cell = input->GetCell(cellId);
        if( cell && cell->GetCellType() == VTK_TETRA ) 
        {
            vtkIdList* cellPtIds = cell->GetPointIds();

            newCellPtIds = vtkSmartPointer<vtkIdList>::New();
            newCellPtIds->InsertNextId( pointIds[ cellPtIds->GetId(0) ]);
            newCellPtIds->InsertNextId( pointIds[ cellPtIds->GetId(1) ]);
            newCellPtIds->InsertNextId( pointIds[ cellPtIds->GetId(2) ]);
            newCellPtIds->InsertNextId( pointIds[ cellPtIds->GetId(3) ]);

            //check if cell contains deleted point ids
            if( newCellPtIds->GetId(0) > -1 &&
                newCellPtIds->GetId(1) > -1 &&
                newCellPtIds->GetId(2) > -1 &&
                newCellPtIds->GetId(3) > -1 )
            {
                newTets->InsertNextCell(newCellPtIds);

                // copy cell data
                for( int i=0; i<numArrays; i++ ){
                    inCD->GetArray(i)->GetTuple(cellId, data);
                    outCD->GetArray(i)->InsertNextTuple(data);
                }
            }
        }
        this->UpdateProgress(0.5+0.5*cellId/numCells);
    }

    output->SetCells(VTK_TETRA, newTets);
    newTets->Delete();
    
    output->Squeeze();

	return 1;
}
