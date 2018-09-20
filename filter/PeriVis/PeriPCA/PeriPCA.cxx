/*=========================================================================

Program:   Visualization Toolkit
Module:    PeriPCA.h

Copyright (c) Michael Bußler

=========================================================================*/

#include "PeriPCA.h"
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

vtkStandardNewMacro(PeriPCA);

inline bool SameSign(float a, float b) {
	return a*b >= 0.0f;
}

inline void swap( int&a, int&b ){
	int swp=a; a=b; b=swp;
};

//==============================================================================
PeriPCA::PeriPCA()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);
}
//==============================================================================
PeriPCA::~PeriPCA()
{
}

//==============================================================================
int PeriPCA::RequestData(
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

	vtkPointData *pd = output->GetPointData();

	vtkIdType numPts, inPtId, inCellId;

	vtkDebugMacro(<<"Calculating Light Scattering");

	vtkPointData *outPD = output->GetPointData();
	numPts = input->GetNumberOfPoints();

	if ( numPts < 1 )
	{
		vtkErrorMacro(<<"No data!");
		return 1;
	}

	vtkDebugMacro(<<"Initialise..");

	// allocate mem for output data arrays
	// array stores intensity per point
	vtkSmartPointer<vtkDoubleArray> eigenValues, covTensor, 
                                    ev0, ev1, ev2,
                                    cl, cp, cs,
                                    normals;

	eigenValues = vtkDoubleArray::New();
	eigenValues->SetName("EigenValues");
	eigenValues->SetNumberOfComponents(3);
	eigenValues->SetNumberOfTuples(numPts);

	covTensor = vtkDoubleArray::New();
	covTensor->SetName("Covariance");
	covTensor->SetNumberOfComponents(9);
	covTensor->SetNumberOfTuples(numPts);

    ev0 = vtkDoubleArray::New();
    ev0->SetName("ev0");
    ev0->SetNumberOfComponents(3);
    ev0->SetNumberOfTuples(numPts);

    ev1 = vtkDoubleArray::New();
    ev1->SetName("ev1");
    ev1->SetNumberOfComponents(3);
    ev1->SetNumberOfTuples(numPts);

    ev2 = vtkDoubleArray::New();
    ev2->SetName("ev2");
    ev2->SetNumberOfComponents(3);
    ev2->SetNumberOfTuples(numPts);

    cl = vtkDoubleArray::New();
    cl->SetName("c_l");
    cl->SetNumberOfComponents(1);
    cl->SetNumberOfTuples(numPts);
    
    cp = vtkDoubleArray::New();
    cp->SetName("c_p");
    cp->SetNumberOfComponents(1);
    cp->SetNumberOfTuples(numPts);

    cs = vtkDoubleArray::New();
    cs->SetName("c_s");
    cs->SetNumberOfComponents(1);
    cs->SetNumberOfTuples(numPts);

    normals = vtkDoubleArray::New();
    normals->SetName("normals");
    normals->SetNumberOfComponents(3);
    normals->SetNumberOfTuples(numPts);

    //for (inPtId=0; inPtId < numPts; inPtId++)
	//	eigenValues->SetValue(inPtId, 0.0f);  // x0	

	input->BuildLinks(); // build up connectivity lookup table
	std::vector<vtkSmartPointer<vtkIdList> > conns;
	conns.reserve(numPts);

	// calculate and store bond connections per particle
	for (inPtId=0; inPtId < numPts; inPtId++)
	{
		// get neighbors for diffusion kernel calculations
		vtkSmartPointer<vtkIdList> connectedVertices =
			GetConnectedVertices(input, inPtId);
		conns.push_back(connectedVertices);
	}
	
	vtkDebugMacro(<<"Calculate PCA..");

	// calculate PCA per particle
	for (inPtId=0; inPtId < numPts; inPtId++) {

		double pos[3];
		input->GetPoint(inPtId, pos);

		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		points->InsertNextPoint(pos);

		const vtkSmartPointer<vtkIdList>& bondedPoints = conns[inPtId];

		for( vtkIdType i=0; i<bondedPoints->GetNumberOfIds(); i++ )
		{
			vtkIdType bondedPointId = bondedPoints->GetId(i);
			input->GetPoint( bondedPointId, pos);
			points->InsertNextPoint(pos);
		}

		double ew[3];
		double cov[9];
        double ev[9];
		
        calculatePCA(points, ew, ev, cov);
		
        eigenValues->SetTuple(inPtId, ew);
		covTensor->SetTuple(inPtId, cov);
        ev0->SetTuple(inPtId, &ev[0]);
        ev1->SetTuple(inPtId, &ev[3]);
        ev2->SetTuple(inPtId, &ev[6]);

        double c_l =     (ew[0]-ew[1]) / (ew[0]+ew[1]+ew[2]);
        double c_p = 2.0*(ew[1]-ew[2]) / (ew[0]+ew[1]+ew[2]);
        double c_s = 3.0*(ew[2])       / (ew[0]+ew[1]+ew[2]);

        cl->SetValue( inPtId, c_l);
        cp->SetValue( inPtId, c_p);
        cs->SetValue( inPtId, c_s);

        double normal[3];
        if (c_l>0.2) {
            normal[0] = ev[0];
            normal[1] = ev[1];
            normal[2] = ev[2];
        } else if (c_p>0.3) {
            normal[0] = ev[3];
            normal[1] = ev[4];
            normal[2] = ev[5];
        } else {
            normal[0] = ev[6];
            normal[1] = ev[7];
            normal[2] = ev[8];
        }

        normals->SetTuple(inPtId, normal);
	}
	
	//vtkPointData *pd = output->GetPointData();
	//if( pd->HasArray("Intensity"))
	//	pd->RemoveArray("Intensity");
	//pd->AddArray(intensityIn);

	pd->AddArray(eigenValues);
	pd->AddArray(covTensor);
    pd->AddArray(ev0);
    pd->AddArray(ev1);
    pd->AddArray(ev2);
    pd->AddArray(cl);
    pd->AddArray(cp);
    pd->AddArray(cs);
    pd->AddArray(normals);
    
	output->Squeeze();

	return 1;
}

//==============================================================================
vtkSmartPointer<vtkIdList> 
	PeriPCA::GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int id)
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

void PeriPCA::calculatePCA( vtkPoints* points, double eigenValues[3], double* ev/*=NULL */, double* covariance/*=NULL */ )
{
	vtkIdType numPoints = points->GetNumberOfPoints()-1; // skip p0 (center)
	std::vector<double> evSet[3];

	evSet[0].resize(numPoints);
	evSet[1].resize(numPoints);
	evSet[2].resize(numPoints);

	// estimate mean position
	//double pos[3], u[3] = {0.0, 0.0, 0.0};
	//for( vtkIdType ptId = 0; ptId < numPoints; ptId++ ) {
	//	points->GetPoint(ptId, pos);
	//	u[0] += pos[0]; u[1] += pos[1]; u[2] += pos[2];
	//}	
	//u[0] /= numPoints; u[1] /= numPoints; u[2] /= numPoints;

	// use particle position as center..
	double pos[3], u[3];
	points->GetPoint(0, u);
	for( vtkIdType ptId = 1; ptId < numPoints; ptId++ )
	{
		points->GetPoint(ptId, pos);
		evSet[0][ptId] = pos[0]-u[0];
		evSet[1][ptId] = pos[1]-u[1];
		evSet[2][ptId] = pos[2]-u[2];
	}	

	//Covariance Matrix
	mat3 m;
	for(int j=0; j<3; ++j) {
		for(int i=0; i<3; ++i) {
			double sum = 0.0;
			for(int n=0; n<numPoints; ++n){
				sum += evSet[i][n] * evSet[j][n];
			}
			m[i][j] = sum/(double)(numPoints-1);
		}
	}

	vec3 lambda;
	mat3eigenvalues(m, lambda);

	// find largest EW
	int n=3;
	int sorted[3];
	for( int i=0; i<n; i++) {
		sorted[i] = i;
	}

	// sort by descending order
	bool swapped=true;
	while(swapped) 
	{
		swapped=false;
		for( int i=0; i<(n-1); i++) {
			if( lambda[sorted[i]] < lambda[sorted[i+1]]) {
				swap( sorted[i], sorted[i+1]);
				swapped = true;
			}
		}
	}

	eigenValues[0] = lambda[sorted[0]];
	eigenValues[1] = lambda[sorted[1]];
	eigenValues[2] = lambda[sorted[2]];

	if( covariance ){
		for(int j=0; j<3; ++j) {
			for(int i=0; i<3; ++i) {
				covariance[3*j+i] = m[i][j];
			}
		}
	}

    if( ev ){
        vec3 v;
        for( int i=0; i<3; i++){
            mat3realEigenvector(m, lambda[sorted[i]], v);
            ev[3*i+0] = v[0];
            ev[3*i+1] = v[1];
            ev[3*i+2] = v[2];
        }
    }
}
//==============================================================================
void PeriPCA::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
}

int PeriPCA::ProcessRequest(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	return vtkPolyDataAlgorithm::ProcessRequest(request, inputVector, outputVector);
}

//==============================================================================
void PeriPCA::GetInvariants(double* m, double* pqr) 
{
	pqr[0] = -m[0*3+0]*m[1*3+1]*m[2*3+2]
	-m[0*3+1]*m[1*3+2]*m[2*3+0]
	-m[0*3+2]*m[1*3+0]*m[2*3+1]
	+m[2*3+0]*m[1*3+1]*m[0*3+2]
	+m[2*3+1]*m[1*3+2]*m[0*3+0]
	+m[2*3+2]*m[1*3+0]*m[0*3+1];


	pqr[1] = m[1*3+1]*m[2*3+2] - m[2*3+1]*m[1*3+2]
	+ m[0*3+0]*m[1*3+1] - m[1*3+0]*m[0*3+1]
	+ m[0*3+0]*m[2*3+2] - m[0*3+2]*m[2*3+0];

	pqr[2] = -(m[0*3+0]+m[1*3+1]+m[2*3+2]);
}


int PeriPCA::CalculateRoots(double* a, double* r) 
{
	double c1 = a[1] - a[2]*a[2]/3.;
	double c0 = a[0] - a[1]*a[2]/3. + 2./27.*a[2]*a[2]*a[2];
	// Make cubic coefficient 4 and linear coefficient +- 3
	// by substituting y = z*k and multiplying with 4/k^3
	if (c1 == 0) {
		if (c0 == 0)     r[0] =  0;
		else if (c0 > 0) r[0] = -pow(c0, 1./3.);
		else             r[0] =  pow(-c0, 1./3.);
	}
	else {
		bool negc1 = c1 < 0;
		double absc1 = negc1 ? -c1 : c1;
		double k = sqrt(4./3.*absc1);
		double d0 = c0 * 4./(k*k*k);
		// Find the first solution
		if (negc1) {
			if (d0 > 1)       r[0] = -cosh(acosh(d0)/3);
			else if (d0 > -1) r[0] = -cos(acos(d0)/3);
			else              r[0] =  cosh(acosh(-d0)/3);
		}
		else {
			r[0] = -sinh(asinh(d0)/3);
		}
		// Transform back
		r[0] *= k;
	}
	r[0] -= a[2]/3;
	// Other two solutions
	double p = r[0] + a[2];
	double q = r[0]*p + a[1];
	double discrim = p*p - 4*q;
	//    if (forceReal && discrim < 0.0) discrim = 0.0;
	if (discrim >= 0) {
		double root = sqrt(discrim);
		r[1] = (-p - root)/2.;
		r[2] = (-p + root)/2.;
		return 3;
	}
	else {
		double root = sqrt(-discrim);
		r[1] = -p/2;
		r[2] = root/2.;
		return 1;
	}    
}


void PeriPCA::CalculateEigenvector(double* m, double lambda, double* evect) {
	double norm[3];
	double cross[9];
	double red[9];

	red[0] =  m[0]-lambda;
	red[1] =  m[1];
	red[2] =  m[2];
	red[3] =  m[3];
	red[4] =  m[4]-lambda;
	red[5] =  m[5];
	red[6] =  m[6];
	red[7] =  m[7];
	red[8] =  m[8]-lambda;

	CalculateCross(&red[3], &red[6], &cross[0]);
	CalculateCross(&red[6], &red[0], &cross[3]);
	CalculateCross(&red[0], &red[3], &cross[6]);

	norm[0] = NormSquared(&cross[0]);
	norm[1] = NormSquared(&cross[3]);
	norm[2] = NormSquared(&cross[6]);

	int best = getLargest(norm);
	double len = sqrt(norm[best]);

	if (len > 0) {
		evect[0] = cross[best*3+0] / len;
		evect[1] = cross[best*3+1] / len;
		evect[2] = cross[best*3+2] / len;
	}
	else {
		evect[0] = 0.0;
		evect[1] = 0.0;
		evect[2] = 0.0;
	}
}


void PeriPCA::CalculateCross(double* v1, double* v2, double* res) {
	res[0] = v1[1]*v2[2] - v1[2]*v2[1];
	res[1] = v1[2]*v2[0] - v1[0]*v2[2];
	res[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

double PeriPCA::NormSquared(double* v) {
	return (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}


int PeriPCA::getLargest(double* v) 
{
	if (v[0]>v[1]) {
		if (v[0]>v[2]) {
			return 0;
		} else {
			return 2;
		}
	} else {
		if (v[1]>v[2]) {
			return 1;
		} else {
			return 2;
		}
	}
}
