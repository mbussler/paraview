#include "ImplicitPolydataClip.h"

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
#include <vtkImplicitPolyDataDistance.h>

#include <math.h>
#include "vtkSmartPointer.h"
#include <algorithm>

vtkStandardNewMacro(ImplicitPolydataClip);

//----------------------------------------------------------------------------
ImplicitPolydataClip::ImplicitPolydataClip(vtkImplicitFunction *cf)
{
	this->SetInsideOut(false);

	this->SetNumberOfInputPorts(2);
	this->SetNumberOfOutputPorts(1);
}

//----------------------------------------------------------------------------
ImplicitPolydataClip::~ImplicitPolydataClip()
{
}

//----------------------------------------------------------------------------
//
// Clip through data generating surface.
//
int ImplicitPolydataClip::RequestData(
	vtkInformation *vtkNotUsed(request),
	vtkInformationVector **inputVector,
	vtkInformationVector *outputVector)
{
	// get the info objects
	vtkInformation *inDataInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *inClipPolyInfo = inputVector[1]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// get the input and output
	vtkPolyData *input = vtkPolyData::SafeDownCast(
		inDataInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *inClipPD = vtkPolyData::SafeDownCast(
		inClipPolyInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *output = vtkPolyData::SafeDownCast(
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

	vtkDebugMacro(<< "Clipping polygonal data");

	// Initialize self; create output objects
	//
	if ( numPts < 1 || inPts == NULL )
	{
		vtkDebugMacro(<<"No data to clip");
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
	}
    if( outPD->HasArray("Normals"))
        outPD->SetActiveNormals("Normals");
    if( outPD->HasArray("TextureCoordinates"))
        outPD->SetActiveTCoords("TextureCoordinates");

	// clip point data
	newPoints =  vtkPoints::New();
	//newVerts = vtkCellArray::New();
	std::vector<vtkIdType> pointIds;

	double pos[3];
	double data[9];
	vtkIdType vertexId = 0;

	for( vtkIdType ptId=0; ptId<numPts; ptId++) {

		inPts->GetPoint(ptId, pos);
		float signedDistance = implicitPolyDataDistance->EvaluateFunction(pos);

		if(( this->InsideOut && signedDistance >= 0.0) || // point inside polygon
		   (!this->InsideOut && signedDistance <= 0.0))
		{
			vtkIdType newPtId = newPoints->InsertNextPoint(pos);
			//newVerts->InsertNextCell(1, &newPtId);

			// copy point data
			for( int i=0; i<numArrays; i++ ){
				inPD->GetArray(i)->GetTuple(ptId, data);
				outPD->GetArray(i)->InsertNextTuple(data);
			}

			pointIds.push_back(ptId);
		}

		this->UpdateProgress(0.5*ptId/numPts);
	}
	output->SetPoints(newPoints);
	//output->SetVerts(newVerts);

	newPoints->Delete();
	//newVerts->Delete();

 	// copy cell data arrays
    numArrays = inCD->GetNumberOfArrays();
	for( int i=0; i<numArrays; i++ ){
		vtkDataArray *arr = vtkDataArray::CreateDataArray(inCD->GetArray(i)->GetDataType());;
		arr->SetNumberOfComponents(inCD->GetArray(i)->GetNumberOfComponents());
		arr->SetName(inCD->GetArray(i)->GetName());
        outCD->AddArray(arr);
	}

    cellId = 0;

    //// Process Vertices
    vtkCellArray* verts = input->GetVerts();

    //// Process Lines
    cellId = input->GetNumberOfVerts();
    vtkCellArray* lines = input->GetLines();
    if( lines && lines->GetNumberOfCells() > 0 ) {
        vtkSmartPointer<vtkCellArray> newLines = vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();
        lines->InitTraversal();
        while( lines->GetNextCell(pts) ) 
        {
            int i=0;
            while( i < pts->GetNumberOfIds() )
            {
                vtkSmartPointer<vtkIdList> newids = vtkSmartPointer<vtkIdList>::New();
                std::vector<vtkIdType>::iterator it;

                for( ;i<pts->GetNumberOfIds(); i++)
                {
                    int a = pts->GetId(i);
                    it = std::find(pointIds.begin(), pointIds.end(), a);
                    if( it != pointIds.end())
                    {
                        newids->InsertNextId(it-pointIds.begin());
                    } 
                    else 
                    {
                        i++;
                        break;
                    }
                }

                if( newids->GetNumberOfIds() > 1) 
                {
                    newLines->InsertNextCell(newids);

                    // copy cell data
			        for( int i=0; i<numArrays; i++ ){
                        inCD->GetArray(i)->GetTuple(cellId, data);
				        outCD->GetArray(i)->InsertNextTuple(data);
			        }
                }
            }
            
            /*
            int a0 = pts->GetId(0), 
                a1 = pts->GetId(1);
            
			std::vector<vtkIdType>::iterator it0, it1;
            it0 = std::find (pointIds.begin(), pointIds.end(), a0);
			if( it0 != pointIds.end() )
			{
				it1 = std::find (pointIds.begin(), pointIds.end(), a1);
				if( it1 != pointIds.end()) 
				{
					newLines->InsertNextCell(2);
					newLines->InsertCellPoint(it0-pointIds.begin());
					newLines->InsertCellPoint(it1-pointIds.begin());

                    // copy cell data
			        for( int i=0; i<numArrays; i++ ){
                        inCD->GetArray(i)->GetTuple(cellId, data);
				        outCD->GetArray(i)->InsertNextTuple(data);
			        }
				}
			}
            */
            cellId++;
        }
        output->SetLines(newLines);
        this->UpdateProgress(0.5+0.5*cellId/numCells);
    }

    //// Process Polygons
    cellId = input->GetNumberOfVerts() + input->GetNumberOfLines();

    vtkCellArray* polys = input->GetPolys();
    if( polys && polys->GetNumberOfCells() > 0  ) 
    {
        vtkSmartPointer<vtkCellArray> newPolys = vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();
        polys->InitTraversal();
        while( polys->GetNextCell(pts)) 
        {
            int numPoints = pts->GetNumberOfIds();

            std::vector<vtkIdType>::iterator it;
            vtkSmartPointer<vtkIdList> newIds = vtkSmartPointer<vtkIdList>::New();

            bool valid = true;
            for( int i=0; i<numPoints; i++) 
            {
                it = std::find (pointIds.begin(), pointIds.end(), pts->GetId(i));
                valid &= (it != pointIds.end());
                newIds->InsertNextId( it-pointIds.begin() );
            }

            if( valid )
            {
                newPolys->InsertNextCell( newIds);

			    // copy cell data
			    for( int i=0; i<numArrays; i++ ){
                    inCD->GetArray(i)->GetTuple(cellId, data);
				    outCD->GetArray(i)->InsertNextTuple(data);
			    }
            }
            
            cellId++;
        }
        output->SetPolys(newPolys);
        this->UpdateProgress(0.5+0.5*cellId/numCells);
    }

	//vtkDebugMacro(<<"Created: "
	//	<< newPoints->GetNumberOfPoints() << " points, "
	//	<< newLines->GetNumberOfCells() << " lines, " );

	output->Squeeze();

	return 1;
}
