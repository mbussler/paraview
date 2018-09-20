#include "ImplicitTetClip.h"

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

/* performance measure */
#include "timer.h"
#include <QElapsedTimer>

vtkStandardNewMacro(ImplicitTetClip);

//----------------------------------------------------------------------------
ImplicitTetClip::ImplicitTetClip(vtkImplicitFunction *cf)
{
	this->SetInsideOut(false);
    this->ClipDistance = 0.0;

	this->SetNumberOfInputPorts(2);
	this->SetNumberOfOutputPorts(1);
}

int ImplicitTetClip::FillInputPortInformation(int port, vtkInformation* info)
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
ImplicitTetClip::~ImplicitTetClip()
{
}

//----------------------------------------------------------------------------
//
// Clip through data generating surface.
//
int ImplicitTetClip::RequestData(
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
		outInfo->Get(vtkUnstructuredGrid::DATA_OBJECT()));

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

    vtkCellArray* newTets = vtkCellArray::New();

    QElapsedTimer timer;
    timer.start();

    if( this->ClipByCellEdges ) {

        // copy all points and data
        output->DeepCopy(input);
        
        // clear cell data arrays
        int numArrays = inCD->GetNumberOfArrays();
        for( int i=0; i<numArrays; i++ ){
            outCD->GetArray(i)->SetNumberOfTuples(0);
        }

        cellId = 0;

        //// Process Tetrahedra
        vtkSmartPointer<vtkIdList> newCellPtIds;

        // generate polygons from input poly data mesh
        std::vector<vtkSmartPointer<vtkPolygon> > polys;
        vtkCellArray* tris_ref = inClipPD->GetPolys();
        tris_ref->InitTraversal();
        vtkSmartPointer<vtkIdList> cellPts = vtkSmartPointer<vtkIdList>::New();
        vtkIdType cellId = 0;
        while( tris_ref->GetNextCell(cellPts) ) 
        {
            vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
            vtkSmartPointer<vtkPoints> points = polygon->GetPoints();
            polygon->GetPointIds()->SetNumberOfIds(cellPts->GetNumberOfIds());

            double p[3];
            for( int ptId=0; ptId<cellPts->GetNumberOfIds(); ptId++)
            {
                inClipPD->GetPoint(cellPts->GetId(ptId), p);
                points->InsertNextPoint(p);
                polygon->GetPointIds()->SetId(ptId, ptId);
            }
            polys.push_back(polygon);
        }

        std::vector<std::pair<int,int> > edges;
        edges.push_back(std::pair<int,int>(0,1));
        edges.push_back(std::pair<int,int>(0,2));
        edges.push_back(std::pair<int,int>(0,3));
        edges.push_back(std::pair<int,int>(1,2));
        edges.push_back(std::pair<int,int>(1,3));
        edges.push_back(std::pair<int,int>(2,3));

        // generate implicit function from poly data 
        vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance =
            vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
        implicitPolyDataDistance->SetInput(inClipPD);

        // process point data
        std::vector<bool> boundaryPoints;
        boundaryPoints.resize(numPts, false);

        outPD = output->GetPointData();
        vtkSmartPointer<vtkIntArray> boundaryPointArray = vtkSmartPointer<vtkIntArray>::New();
        boundaryPointArray->SetName("Boundary Point");
        boundaryPointArray->SetNumberOfComponents(1);
        boundaryPointArray->SetNumberOfTuples(numPts);
        outPD->AddArray(boundaryPointArray);

        double pos[3];
        for( vtkIdType ptId=0; ptId<numPts; ptId++) {

            inPts->GetPoint(ptId, pos);
            double signedDistance = implicitPolyDataDistance->EvaluateFunction(pos);

            if(( this->InsideOut && signedDistance <= this->ClipDistance) || // point inside polygon
                (!this->InsideOut && signedDistance >= -this->ClipDistance))
            {
                boundaryPoints[ptId] = true;
                boundaryPointArray->SetValue(ptId, true);
            }
            else
            {
                boundaryPointArray->SetValue(ptId, false);
            }

            this->UpdateProgress(0.5*ptId/numPts);
        }

        std::cout << std::endl;
        std::cout << std::setprecision(3);

        // iterate over all cells, insert cells which edges don't intersect polygon
        for( cellId=0; cellId<input->GetNumberOfCells(); cellId++) 
        {
            vtkCell* cell = input->GetCell(cellId);
            if( cell && cell->GetCellType() == VTK_TETRA ) 
            {
                vtkIdList* cellPtIds = cell->GetPointIds();

                // check if cell has only boundary point
                bool hasBoundaryPoint = false;
                bool allBoundaryPoint = true;
                for( vtkIdType ptId = 0; ptId < cellPtIds->GetNumberOfIds(); ptId++)
                {
                    vtkIdType cellPtId = cellPtIds->GetId(ptId);
                    hasBoundaryPoint |= boundaryPoints[cellPtId];
                    allBoundaryPoint &= boundaryPoints[cellPtId];
                }
                
                bool intersect = false;

                //if ( hasBoundaryPoint )
                //{
                //    std::vector<std::pair<int,int> >::iterator edgesIter;
                //    for( edgesIter = edges.begin(); edgesIter != edges.end(); edgesIter++)
                //    {
                //        int ptId1 = edgesIter->first;
                //        int ptId2 = edgesIter->second;
                //        if( boundaryPoints[ptId1] || boundaryPoints[ptId2]) {
                //            intersect |= IntersectWithPolys(output, cellPtIds, ptId1, ptId2, polys);
                //            if( intersect) {
                //                break;
                //            }
                //        }
                //    }
                //}

                //if( !intersect )
                if( !allBoundaryPoint )
                {
                    newTets->InsertNextCell(cellPtIds);

                    // copy cell data
                    double data[9];
                    for( int i=0; i<numArrays; i++ ){
                        inCD->GetArray(i)->GetTuple(cellId, data);
                        outCD->GetArray(i)->InsertNextTuple(data);
                    }
                }
            }

            this->UpdateProgress(0.5+0.5*cellId/numCells);
            //if( cellId % 1000 == 0 )
            //    std::cout << "\r\t\t\t\t\r" << cellId/(float)numCells << "%";
        }

    } else {

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
		    double signedDistance = implicitPolyDataDistance->EvaluateFunction(pos);

		    if(( this->InsideOut && signedDistance >= this->ClipDistance) || // point inside polygon
		       (!this->InsideOut && signedDistance <= -this->ClipDistance))
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
    }
    
    write_timer("ImplicitTetClip", "clip", timer.elapsed());

    output->SetCells(VTK_TETRA, newTets);
    newTets->Delete();

    output->Squeeze();

	return 1;
}

bool ImplicitTetClip::IntersectWithPolys(vtkUnstructuredGrid * output, vtkIdList* cellPtIds, 
                                         int ptId1, int ptId2, 
                                         std::vector<vtkSmartPointer<vtkPolygon> > &polys)
{
    double p1[3], p2[3];
    double tolerance = this->Tolerance;
    double t;
    double x[3];
    double pcoords[3];
    int subId;

    output->GetPoint(cellPtIds->GetId(ptId1), p1);
    output->GetPoint(cellPtIds->GetId(ptId2), p2);

    std::vector<vtkSmartPointer<vtkPolygon> >::iterator polysIter;
    for( polysIter = polys.begin(); polysIter != polys.end(); polysIter++)
    {
        vtkSmartPointer<vtkPolygon> poly = *polysIter;
        if( 0 != poly->IntersectWithLine(p1, p2, tolerance, t, x, pcoords, subId))
            return true;
    }
    return false;
}
