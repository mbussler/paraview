/*=========================================================================

Program:   Visualization Toolkit
Module:    PeriStress.h

Copyright (c) Michael Bußler

=========================================================================*/

#include "PeriStress.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkDataSet.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkTransform.h"
#include "vtkSmartPointer.h"
#include "vtkMath.h"

#include "linalg.h"

vtkStandardNewMacro(PeriStress);

inline bool SameSign(float a, float b) {
	return a*b >= 0.0f;
}

inline void swap( int&a, int&b ){
	int swp=a; a=b; b=swp;
};

//==============================================================================
PeriStress::PeriStress()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);
}
//==============================================================================
PeriStress::~PeriStress()
{
}

//==============================================================================
int PeriStress::RequestData(
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

	output->ShallowCopy(input);

	vtkPointData *inPD = input->GetPointData();

	vtkIdType numPts, inPtId, inCellId;

	vtkDebugMacro(<<"Calculating Stress per Atom");

	vtkPointData *outPD = output->GetPointData();
	numPts = input->GetNumberOfPoints();

	if ( numPts < 1 )
	{
		vtkErrorMacro(<<"No data!");
		return 1;
	}

	vtkDebugMacro(<<"Initialise..");

    // get input arrays
    vtkSmartPointer<vtkFloatArray> v[3];
    v[0] = vtkFloatArray::SafeDownCast(inPD->GetArray("vx"));
    v[1] = vtkFloatArray::SafeDownCast(inPD->GetArray("vy"));
    v[2] = vtkFloatArray::SafeDownCast(inPD->GetArray("vz"));

    vtkSmartPointer<vtkFloatArray> f[3];
    f[0] = vtkFloatArray::SafeDownCast(inPD->GetArray("fx"));
    f[1] = vtkFloatArray::SafeDownCast(inPD->GetArray("fy"));
    f[2] = vtkFloatArray::SafeDownCast(inPD->GetArray("fz"));


    input->BuildLinks(); // build up connectivity lookup table

    vtkSmartPointer<vtkFloatArray> stress = vtkSmartPointer<vtkFloatArray>::New();
    stress->SetName("stress");
    stress->SetNumberOfComponents(9);
    stress->SetNumberOfTuples(numPts);

	// calculate and store per atom stress values
	
	vtkDebugMacro(<<"Calculate stress per atom..");

	// calculate PCA per particle
	for (inPtId=0; inPtId < numPts; inPtId++) {

		double p1[3];
		input->GetPoint(inPtId, p1);
        double m = 1;

        double v0 = v[0]->GetValue(inPtId);
        double v1 = v[1]->GetValue(inPtId);
        double v2 = v[2]->GetValue(inPtId);

        double f1[3];
        f1[0] = f[0]->GetValue(inPtId);
        f1[1] = f[1]->GetValue(inPtId);
        f1[2] = f[2]->GetValue(inPtId);

        float S[9];
        S[0] = -m*v0*v0; // S_xx
        S[4] = -m*v1*v1; // S_yy
        S[8] = -m*v2*v2; // S_zz
        S[1] = -m*v0*v1; // S_xy
        S[2] = -m*v0*v2; // S_xz
        S[5] = -m*v1*v2; // S_yz

		const vtkSmartPointer<vtkIdList>& bondedPoints = GetConnectedVertices(input, inPtId);
		for( vtkIdType i=0; i<bondedPoints->GetNumberOfIds(); i++ )
		{
            vtkIdType bondedPointId = bondedPoints->GetId(i);

            double p2[3];
            input->GetPoint( bondedPointId, p2);

            double f2[3];
            f2[0] = f[0]->GetValue(bondedPointId);
            f2[1] = f[1]->GetValue(bondedPointId);
            f2[2] = f[2]->GetValue(bondedPointId);

            S[0] -= p1[0]*f1[0] + p2[0]*f2[0]; // S_xx
            S[4] -= p1[1]*f1[1] + p2[1]*f2[1]; // S_yy
            S[8] -= p1[2]*f1[2] + p2[2]*f2[2]; // S_zz
                                  
            S[1] -= p1[0]*f1[1] + p2[0]*f2[1]; // S_xy
            S[2] -= p1[0]*f1[2] + p2[0]*f2[2]; // S_xz
            S[5] -= p1[1]*f1[2] + p2[1]*f2[2]; // S_yz

        }

        S[3] = S[1];
        S[6] = S[2];
        S[7] = S[5];

        stress->SetTuple(inPtId, S);
	}

	outPD->AddArray(stress);
    
	output->Squeeze();

	return 1;
}

//==============================================================================
vtkSmartPointer<vtkIdList> 
	PeriStress::GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int id)
{
	vtkSmartPointer<vtkIdList> connectedVertices =
		vtkSmartPointer<vtkIdList>::New();

	//get all cells that vertex 'id' is a part of
	vtkSmartPointer<vtkIdList> cellIdList =
		vtkSmartPointer<vtkIdList>::New();
	mesh->GetPointCells(id, cellIdList);

	/*
	cout << "Vertex 0 is used in cells ";
	for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
	{
	cout << cellIdList->GetId(i) << ", ";
	}
	cout << endl;
	*/

	for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
	{
		//cout << "id " << i << " : " << cellIdList->GetId(i) << endl;

		vtkSmartPointer<vtkIdList> pointIdList =
			vtkSmartPointer<vtkIdList>::New();
		mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);

		//cout << "End points are " << pointIdList->GetId(0) << " and " << pointIdList->GetId(1) << endl;

		if(pointIdList->GetId(0) != id)
		{
			//cout << "Connected to " << pointIdList->GetId(0) << endl;
			connectedVertices->InsertNextId(pointIdList->GetId(0));
		}
		else
		{
			//cout << "Connected to " << pointIdList->GetId(1) << endl;
			connectedVertices->InsertNextId(pointIdList->GetId(1));
		}
	}

	return connectedVertices;
}

//==============================================================================
void PeriStress::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
}

int PeriStress::ProcessRequest(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	return vtkPolyDataAlgorithm::ProcessRequest(request, inputVector, outputVector);
}

//==============================================================================
