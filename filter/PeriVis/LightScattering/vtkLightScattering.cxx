/*=========================================================================

Program:   Visualization Toolkit
Module:    vtkLightScattering.h

Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkLightScattering.h"
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

vtkStandardNewMacro(vtkLightScattering);

//==============================================================================
vtkLightScattering::vtkLightScattering()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->m_conns = NULL;
	this->Diffusion = 0.1;
	this->Stepsize = 0.01;
	this->Iterations = 5;
	this->SeedEvery = 0;
	this->SeedParticleStride = 100;
	this->ResetIntensity = false;
	this->InitialSeed = true;
}
//==============================================================================
vtkLightScattering::~vtkLightScattering()
{
	delete[] m_conns;
}

//==============================================================================
int vtkLightScattering::RequestData(
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
	vtkSmartPointer<vtkDoubleArray> intensityIn;

	if( pd->HasArray("Intensity")) {
		intensityIn = dynamic_cast<vtkDoubleArray*>(pd->GetArray("Intensity"));
	
		if( this->ResetIntensity)
		{
			// reset intensity source term
			for (inPtId=0; inPtId < numPts; inPtId++)
				intensityIn->SetValue(inPtId, 0.0f);  // x0	
		}
	} else {
		 intensityIn = vtkDoubleArray::New();
		 intensityIn->SetName("Intensity");
		 intensityIn->SetNumberOfComponents(1);
		 intensityIn->SetNumberOfTuples(numPts);
		 for (inPtId=0; inPtId < numPts; inPtId++)
			 intensityIn->SetValue(inPtId, 0.0f);  // x0	

		 pd->AddArray(intensityIn);
	}

	if( !m_conns ){

		input->BuildLinks(); // build up connectivity lookup table
		m_conns = new vtkSmartPointer<vtkIdList>[numPts];

		for (inPtId=0; inPtId < numPts; inPtId++) // iterate over points and store connections
		{
			// get neighbors for diffusion kernel calculations
			vtkSmartPointer<vtkIdList> connectedVertices =
				GetConnectedVertices(input, inPtId);
			m_conns[inPtId] = connectedVertices;
		}
	}
	
	// init the resulting intensity
	vtkSmartPointer<vtkDoubleArray> intensity;
	intensity = vtkDoubleArray::New();
	intensity->SetName("IntensitySrc");
	intensity->SetNumberOfComponents(1);
	intensity->SetNumberOfTuples(numPts);

	vtkDebugMacro(<<"Calculate Diffusion..");

	// diffusion constants @TODO: get from UI
	const double tstep = this->Stepsize;
	const double diff = this->Diffusion;
	const double a = tstep*diff*numPts;

	// iterate diffusion
	for( int iteration = 0; iteration < this->Iterations; iteration++) 
	{		
		vtkDebugMacro(<<"Iteration #" << iteration);

		// reset result array
		for (inPtId=0; inPtId < numPts; inPtId++)
			intensity->SetValue(inPtId, 0.0f);

		// Seed
		if(( iteration == 0 && this->InitialSeed) ||			
		   ( this->SeedEvery && ((iteration % this->SeedEvery) == 0)))
		{
			for (inPtId=0; inPtId < numPts; inPtId += this->SeedParticleStride)
			{
				//if( m_conns[inPtId]->GetNumberOfIds() > 5 ) // only particles with at least 5 bonds
				{
					float seedIntensity = 1.0f;
					intensity->SetValue( inPtId, seedIntensity);
					intensityIn->SetValue( inPtId, 
						intensityIn->GetValue(inPtId) + seedIntensity*this->Stepsize ); // add source
				}
			}
		}

		// SOLVE: iterate over points and diffuse light
		for( int k=0; k<20; k++ ) 
		{
			for (inPtId=0; inPtId < numPts; inPtId++)
			{
				// get neighbors for diffusion kernel calculations
				vtkSmartPointer<vtkIdList> connectedVertices = m_conns[inPtId];

				// initial intensity value
				const int nConns = connectedVertices->GetNumberOfIds();
				double val = intensityIn->GetValue(inPtId); // x[IX(i,j)]

				//x0[IX(i,j)] = (x[IX(i,j)] + a*( x0[IX(i-1,j)] +
				//								  x0[IX(i+1,j)] +
				//								  x0[IX(i,j-1)] +
				//								  x0[IX(i,j+1)]))/c;

				for(vtkIdType i = 0; i < connectedVertices->GetNumberOfIds(); i++)
				{
					vtkIdType neighPtId = connectedVertices->GetId(i);
					val += a*intensity->GetValue(neighPtId); // a * x0[IX(i-1,j)]
				}
				val = val / (1.0+nConns*a); // c
				intensity->SetValue(inPtId, val); // set x0[IX(i,j)]
			}

			this->UpdateProgress( (20*iteration+k)/(float) 20*this->Iterations);
			if (this->GetAbortExecute())
			{
				break;
			}
		}

		// copy back x = x0 ...
		for (inPtId=0; inPtId < numPts; inPtId++) {
			intensityIn->SetValue(inPtId, intensity->GetValue(inPtId));
		}
	}
	
	//vtkPointData *pd = output->GetPointData();
	//if( pd->HasArray("Intensity"))
	//	pd->RemoveArray("Intensity");
	//pd->AddArray(intensityIn);

	output->Squeeze();

	return 1;
}

//==============================================================================
vtkSmartPointer<vtkIdList> 
	vtkLightScattering::GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int id)
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
void vtkLightScattering::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
}

int vtkLightScattering::ProcessRequest(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
	return vtkPolyDataAlgorithm::ProcessRequest(request, inputVector, outputVector);
}
