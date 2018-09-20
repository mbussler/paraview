#include "ClipLines.h"

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

#include <math.h>
#include "vtkSmartPointer.h"
#include <algorithm>

vtkStandardNewMacro(ClipLines);
vtkCxxSetObjectMacro(ClipLines,ClipFunction,vtkImplicitFunction);

//----------------------------------------------------------------------------
// Construct with user-specified implicit function; InsideOut turned off; value
// set to 0.0; and generate clip scalars turned off.
ClipLines::ClipLines(vtkImplicitFunction *cf)
{
  this->ClipFunction = cf;
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  this->CrincleClip = 0;
}

//----------------------------------------------------------------------------
ClipLines::~ClipLines()
{
  this->SetClipFunction(NULL);
}

//----------------------------------------------------------------------------
// some stl magic to sort pair by second argument
typedef std::pair<vtkIdType, vtkIdType> vtkIdLink;

struct sort_by_first {
	bool operator()(const vtkIdLink &left, const vtkIdLink &right) {
		return left.first < right.first;
	}
};
struct sort_by_second {
	bool operator()(const vtkIdLink &left, const vtkIdLink &right) {
		return left.second < right.second;
	}
};
struct compare_by_first
{
	compare_by_first(vtkIdType const& val) : _val(val) { }
	bool operator () (vtkIdLink const& lnk) { return (lnk.first == _val); }
	vtkIdType _val;
};
struct compare_pair
{
	bool operator()(vtkIdLink l, size_t s) const { return l.first < s;}
	bool operator()(size_t s, vtkIdLink l) const { return s < l.first;}
};

bool compFunc (vtkIdLink i,vtkIdLink j) { return (i.first<j.first); }

//----------------------------------------------------------------------------
// Overload standard modified time function. If Clip functions is modified,
// then this object is modified as well.
unsigned long ClipLines::GetMTime()
{
  unsigned long mTime=this->Superclass::GetMTime();
  unsigned long time;

  if ( this->ClipFunction != NULL )
    {
    time = this->ClipFunction->GetMTime();
    mTime = ( time > mTime ? time : mTime );
    }
  return mTime;
}

//----------------------------------------------------------------------------
//
// Clip through data generating surface.
//
int ClipLines::RequestData(
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

  if ( !this->ClipFunction )
    {
    vtkErrorMacro(<<"Cannot generate clip scalars if no clip function defined");
    return 1;
    }


  // copy point data arrays
  int numArrays = inPD->GetNumberOfArrays();
  for( int i=0; i<numArrays; i++ ){
	  int type = inPD->GetArray(i)->GetDataType();
	  vtkDataArray *arr = vtkDataArray::CreateDataArray(type);
	  arr->SetNumberOfComponents(inPD->GetArray(i)->GetNumberOfComponents());
	  arr->SetName(inPD->GetArray(i)->GetName());
	  outPD->AddArray(arr);
	  arr->Delete();
  }

  // copy cell data arrays
  int numCellArrays = inCD->GetNumberOfArrays();
  for( int i=0; i<numCellArrays; i++ ){
	  int type = inCD->GetArray(i)->GetDataType();
	  vtkDataArray *arr = vtkDataArray::CreateDataArray(type);
	  arr->SetNumberOfComponents(inCD->GetArray(i)->GetNumberOfComponents());
	  arr->SetName(inCD->GetArray(i)->GetName());
	  outCD->AddArray(arr);
	  arr->Delete();
  }

  // clip point data
  newPoints =  vtkPoints::New();
  //newVerts = vtkCellArray::New();

  // links new point ids to previous ones
  //std::vector<vtkIdLink> ptIdsLink;
  std::vector<vtkIdLink> refPtIdsLink;
  
  vtkIdType vertexId = 0;

  for( vtkIdType ptId=0; ptId<numPts; ptId++) 
  {
	  double pos[3];
		inPts->GetPoint(ptId, pos);
		double v = ClipFunction->FunctionValue(pos);

		if( (  this->InsideOut && (v >  0.0 )) ||
			( !this->InsideOut && (v <= 0.0 )) )
		{
			vtkIdType newId = newPoints->InsertNextPoint(pos);

			// copy point data
			for( int i=0; i<numArrays; i++ ){
				double data[9];
				inPD->GetArray(i)->GetTuple(ptId, data);
				outPD->GetArray(i)->InsertNextTuple(data);
			}

			// store old id as reference
			//ptIdsLink.push_back(    vtkIdLink(newId, ptId));
			refPtIdsLink.push_back( vtkIdLink(ptId, newId));
		}

	  this->UpdateProgress(0.5*ptId/numPts);
  }

  // calculate new cell data
  input->BuildLinks();

  // currently only lines (bonds) are supported
  newLines = vtkCellArray::New();

  // sort point Ids
  std::sort(refPtIdsLink.begin(), refPtIdsLink.end(), sort_by_first());

  if( this->CrincleClip ){

	  size_t numRefPoints = refPtIdsLink.size();
	  for( size_t i=0; i<numRefPoints; i++)
	  {
		  vtkIdType refPtId = refPtIdsLink[i].first;
		  vtkSmartPointer<vtkIdList> refNeighIds = getConnectedPoints(input, refPtId);

		  for( size_t j = 0; j<refNeighIds->GetNumberOfIds(); j++) 
		  {
			  vtkIdType refNeighPtId = refNeighIds->GetId(j);

			  if( !std::binary_search(refPtIdsLink.begin(), refPtIdsLink.end(), refNeighPtId, compare_pair()))
			  {
				  // insert new point
				  double pos[3];
				  inPts->GetPoint(refNeighPtId, pos);
				  vtkIdType newId = newPoints->InsertNextPoint(pos);

				  // copy point data
				  for( size_t k=0; k<numArrays; k++ ){
					  double data[9];
					  inPD->GetArray(k)->GetTuple(refNeighPtId, data);
					  outPD->GetArray(k)->InsertNextTuple(data);
				  }

				  // store old id as reference
				  //ptIdsLink.push_back(    vtkIdLink(newId, ptId));
				  refPtIdsLink.push_back( vtkIdLink(refNeighPtId, newId));
			  }
		  }
	  }

	  // sort point Ids
	  std::sort(refPtIdsLink.begin(), refPtIdsLink.end(), sort_by_first());
  }

  for( size_t i=0; i<refPtIdsLink.size(); i++)
  {
	  vtkIdType refPtId = refPtIdsLink[i].first;
	  vtkSmartPointer<vtkIdList> refCellIds  = vtkSmartPointer<vtkIdList>::New();
	  vtkSmartPointer<vtkIdList> refNeighIds = getConnectedPoints(input, refPtId, refCellIds);

	  // remove bonds smaller than ptId, prevents double bonds, requires sorted point ids
	  vtkSmartPointer<vtkIdList> refNeighIdsCleaned = vtkSmartPointer<vtkIdList>::New();
	  vtkSmartPointer<vtkIdList> refCellIdsCleaned = vtkSmartPointer<vtkIdList>::New();
	  for( size_t j=0; j<refNeighIds->GetNumberOfIds(); j++) 
	  {
		  vtkIdType refNeighPtId = refNeighIds->GetId(j);
		  if( refNeighPtId > refPtId ){
			  refNeighIdsCleaned->InsertNextId(refNeighPtId);
			  refCellIdsCleaned->InsertNextId(refCellIds->GetId(j));
		  }
	  }

	  // iterate over (cleaned) neighbors of current point 
	  // and check if neighbor is in clipped dataset
	  for( size_t j=0; j<refNeighIdsCleaned->GetNumberOfIds(); j++) 
	  {
		  vtkIdType refNeighPtId = refNeighIdsCleaned->GetId(j);

		  for( size_t k=i+1; k<refPtIdsLink.size(); k++) // requires sorted links
		  {
			  if( refNeighPtId == refPtIdsLink[k].first ) 
			  {
				  vtkIdType newPtId1, newPtId2;
				  newPtId1 = refPtIdsLink[i].second;
				  newPtId2 = refPtIdsLink[k].second;

				  newLines->InsertNextCell(2);
				  newLines->InsertCellPoint(newPtId1);
				  newLines->InsertCellPoint(newPtId2);

				  // copy cell data
				  for( int i=0; i<numCellArrays; i++ ){
					  double data[9];
					  inCD->GetArray(i)->GetTuple(refCellIdsCleaned->GetId(j), data);
					  outCD->GetArray(i)->InsertNextTuple(data);
				  }

				  break;
			  }
		  }
	  }
  }  
  
  vtkDebugMacro(<<"Created: "
               << newPoints->GetNumberOfPoints() << " points, "
               << newLines->GetNumberOfCells() << " lines, " );


  output->SetPoints(newPoints);
  output->SetLines(newLines);
  //output->SetVerts(newVerts);
  newPoints->Delete();
  newLines->Delete();
  //newVerts->Delete();
  output->Squeeze();

  return 1;
}


//----------------------------------------------------------------------------
void ClipLines::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  if ( this->ClipFunction )
    {
    os << indent << "Clip Function: " << this->ClipFunction << "\n";
    }
  else
    {
    os << indent << "Clip Function: (none)\n";
    }
}


vtkSmartPointer<vtkIdList> ClipLines::getConnectedPoints( vtkPolyData * pd, vtkIdType ptId, vtkSmartPointer<vtkIdList> bondIds /*= NULL*/)
{
	vtkSmartPointer<vtkIdList> cellIds  = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> neighPts = vtkSmartPointer<vtkIdList>::New();

	pd->GetPointCells( ptId, cellIds);

	for(vtkIdType i = 0; i < cellIds->GetNumberOfIds(); i++)
	{
		vtkIdType cellId = cellIds->GetId(i);

		vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
		pd->GetCellPoints( cellId, pointIdList);

		vtkIdType ptIdN = pointIdList->GetId(0);
		if( ptIdN == ptId)
			ptIdN = pointIdList->GetId(1);

		neighPts->InsertNextId(ptIdN);

		if( bondIds )
			bondIds->InsertNextId(cellId);
	}

	return neighPts;
}