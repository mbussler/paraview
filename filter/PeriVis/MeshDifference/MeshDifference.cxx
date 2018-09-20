#include "MeshDifference.h"

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
#include "vtkDoubleArray.h"

#include <math.h>
#include "vtkSmartPointer.h"
#include <algorithm>

vtkStandardNewMacro(MeshDifference);

//----------------------------------------------------------------------------
// Construct with user-specified implicit function; InsideOut turned off; value
// set to 0.0; and generate clip scalars turned off.
MeshDifference::MeshDifference()
{
    this->SetNumberOfInputPorts(2);

    // by default, process active point tensors
    this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
        vtkDataSetAttributes::SCALARS);
}

//----------------------------------------------------------------------------
MeshDifference::~MeshDifference()
{
}

//----------------------------------------------------------------------------
//
//
int MeshDifference::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
    // get the info objects
    vtkInformation *inInfo1 = inputVector[0]->GetInformationObject(0);
    vtkInformation *inInfo2 = inputVector[1]->GetInformationObject(0);
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    // get the input and output
    vtkPolyData *input1 = vtkPolyData::SafeDownCast(
        inInfo1->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *input2 = vtkPolyData::SafeDownCast(
        inInfo2->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output = vtkPolyData::SafeDownCast(
        outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkDataArray* inData1 = this->GetInputArrayToProcess(0, inputVector);

    if( !inData1 ) {
        vtkErrorMacro(<<"No data.");
        return 1;
    }

    vtkDataArray* inData2 = input2->GetPointData()->GetArray( inData1->GetName());

    if( !inData2 ) {
        vtkErrorMacro(<<"Data array "<<inData1->GetName()<<" not found in second input.");
        return 1;
    }

    vtkIdType numPts=input1->GetNumberOfPoints();
    vtkPoints *inPts=input1->GetPoints();

    if ( numPts < 1 || !inData2 )
    {
        vtkDebugMacro(<<"No data");
        return 1;
    }

    // allocate mem for output data arrays
    vtkSmartPointer<vtkFloatArray> diffArray = vtkSmartPointer<vtkFloatArray>::New();
    diffArray->SetName("difference");
    diffArray->SetNumberOfComponents(inData1->GetNumberOfComponents());
    diffArray->SetNumberOfTuples(numPts);

    double pos[3];
    double data1[9], data2[9], diff[9];
    for( vtkIdType ptId=0; ptId<numPts; ptId++) 
    {
        inPts->GetPoint(ptId, pos);
        inData1->GetTuple(ptId, data1);
        //std::cout << "Search position: "<<pos[0]<<","<<pos[1]<<","<<pos[2]<<std::endl;

        double dataSum[9]={0,0,0,0,0,0,0,0,0};
        diffArray->SetTuple(ptId, dataSum);

        int subId;
        double pcoords[3];
        double weights[9];

        vtkCell* cell = input2->FindAndGetCell(pos, 0, 0, 0.0001, subId, pcoords, weights);
        if( cell ) {

            //std::cout << "Cell found!" << std::endl;

            vtkIdList* pts = cell->GetPointIds();
            for( int i=0; i<pts->GetNumberOfIds(); i++)
            {
                input2->GetPoint(pts->GetId(i), pos);
                //std::cout << "cell point "<<i<<": "<<pos[0]<<","<<pos[1]<<","<<pos[2]<<"; w="<<weights[i]<<std::endl;
                
                inData2->GetTuple( pts->GetId(i), data2);
                for( int j=0; j<inData2->GetNumberOfComponents(); j++) {
                    dataSum[j] += weights[i]*data2[j];
                }
            }
            
            for( int i=0; i<inData1->GetNumberOfComponents(); i++) {
                diff[i] = data1[i]-dataSum[i];
                //std::cout << diff[i] << ",";
            }
            //std::cout << std::endl;

            diffArray->SetTuple(ptId, diff);
        }

    }

    output->ShallowCopy(input1);
    output->GetPointData()->AddArray(diffArray);

    return 1;
}


//----------------------------------------------------------------------------
void MeshDifference::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}

