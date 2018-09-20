#include "ClipPolyData.h"

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

vtkStandardNewMacro(ClipPolyData);
vtkCxxSetObjectMacro(ClipPolyData,ClipFunction,vtkImplicitFunction);

//----------------------------------------------------------------------------
// Construct with user-specified implicit function; InsideOut turned off; value
// set to 0.0; and generate clip scalars turned off.
ClipPolyData::ClipPolyData(vtkImplicitFunction *cf)
{
    this->ClipFunction = cf;
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
    this->CrincleClip = 0;
}

//----------------------------------------------------------------------------
ClipPolyData::~ClipPolyData()
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
unsigned long ClipPolyData::GetMTime()
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
int ClipPolyData::RequestData(
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
    vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
    std::vector<vtkIdType> newPtIds;
    newPtIds.resize(numPts, -1);

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

            newPtIds[ptId] = newId;
        }

        this->UpdateProgress(0.5*ptId/numPts);
    }

    output->SetPoints(newPoints);
    output->Allocate();

    vtkSmartPointer<vtkIdList> newCellPtIds;

    // iterate over all cells, insert cells that do not belong to surface
    for( cellId=0; cellId<input->GetNumberOfCells(); cellId++) 
    {
        vtkCell* cell = input->GetCell(cellId);

        if( cell ){
            vtkIdList* cellPtIds = cell->GetPointIds();
            int numPtIds = cellPtIds->GetNumberOfIds();
            newCellPtIds = vtkSmartPointer<vtkIdList>::New();

            bool active = true;
            for( int i=0; i<numPtIds; i++) {
                int newPtId = newPtIds[cellPtIds->GetId(i)];
                newCellPtIds->InsertNextId( newPtId );
                active &= (newPtId > 0);
            }

            if( active ) {
                output->InsertNextCell(cell->GetCellType(), newCellPtIds);

                // copy cell data
                for( int i=0; i<numCellArrays; i++ ){
                    double data[9];
                    inCD->GetArray(i)->GetTuple(cellId, data);
                    outCD->GetArray(i)->InsertNextTuple(data);
                }
            }

        }
   }

    output->Squeeze();

    return 1;
}


//----------------------------------------------------------------------------
void ClipPolyData::PrintSelf(ostream& os, vtkIndent indent)
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


vtkSmartPointer<vtkIdList> ClipPolyData::getConnectedPoints( vtkPolyData * pd, vtkIdType ptId, vtkSmartPointer<vtkIdList> bondIds /*= NULL*/)
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