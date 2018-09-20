#include "TetClip.h"

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
#include "vtkUnstructuredGrid.h"
#include "vtkTriangle.h"
#include "vtkIncrementalPointLocator.h"

#include <math.h>
#include "vtkSmartPointer.h"
#include <algorithm>

/* performance measure */
#include "timer.h"
#include <QElapsedTimer>

vtkStandardNewMacro(TetClip);
vtkCxxSetObjectMacro(TetClip,ClipFunction,vtkImplicitFunction);

//----------------------------------------------------------------------------
// Construct with user-specified implicit function; InsideOut turned off; value
// set to 0.0; and generate clip scalars turned off.
TetClip::TetClip(vtkImplicitFunction *cf)
{
    this->ClipFunction = cf;
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
    this->CrincleClip = 0;
}

//----------------------------------------------------------------------------
TetClip::~TetClip()
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
unsigned long TetClip::GetMTime()
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
int TetClip::RequestData(
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
    vtkCellArray *newTets, *connList=NULL;
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
    
    // list to store new point ids
    vtkSmartPointer<vtkIdList> newPtIds = vtkSmartPointer<vtkIdList>::New();
    newPtIds->Resize(numPts);

    QElapsedTimer timer;
    timer.start();

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

            // store new point id as reference
            newPtIds->SetId(ptId, newId);

        } else {
            // mark point as deleted
            newPtIds->SetId(ptId, -1);
        }

        this->UpdateProgress(0.5*ptId/numPts);
    }

    output->SetPoints(newPoints);

    // new tet cells
    newTets = vtkCellArray::New();
    vtkSmartPointer<vtkIdList> newCellPtIds;

    // iterate over all cells, insert cells that do not belong to surface
    for( cellId=0; cellId<input->GetNumberOfCells(); cellId++) 
    {
        vtkCell* cell = input->GetCell(cellId);
        if( cell && cell->GetCellType() == VTK_TETRA ) 
        {
            vtkIdList* cellPtIds = cell->GetPointIds();

            newCellPtIds = vtkSmartPointer<vtkIdList>::New();
            newCellPtIds->InsertNextId( newPtIds->GetId( cellPtIds->GetId(0) ));
            newCellPtIds->InsertNextId( newPtIds->GetId( cellPtIds->GetId(1) ));
            newCellPtIds->InsertNextId( newPtIds->GetId( cellPtIds->GetId(2) ));
            newCellPtIds->InsertNextId( newPtIds->GetId( cellPtIds->GetId(3) ));

            //check if cell contains deleted point ids
            if( newCellPtIds->GetId(0) > -1 &&
                newCellPtIds->GetId(1) > -1 &&
                newCellPtIds->GetId(2) > -1 &&
                newCellPtIds->GetId(3) > -1 )
            {
                newTets->InsertNextCell(newCellPtIds);
            }
        }
    }

    write_timer("TetClip", "clip", timer.elapsed());

    output->SetCells(VTK_TETRA, newTets);
    output->Squeeze();

    vtkDebugMacro(<<"Created: "
        << newPoints->GetNumberOfPoints() << " points, "
        << newTets->GetNumberOfCells() << " tets, " );


    newPoints->Delete();
    newTets->Delete();
    output->Squeeze();

    return 1;
}


//----------------------------------------------------------------------------
void TetClip::PrintSelf(ostream& os, vtkIndent indent)
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

