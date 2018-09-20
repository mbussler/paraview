#pragma once

#ifndef TIMESTEP_CONTROLLER_CUDA_H
#define TIMESTEP_CONTROLLER_CUDA_H

#include "InlineFunctions.h"
#include <boost/shared_ptr.hpp>
#include <boost/dynamic_bitset.hpp>

#include "TimestepCuda.h"
#include "Typedefs.h"
#include "Vector.h"

#include <list>
#include <vector>

#include "common.cuh"
#include "Particles.cuh"

using namespace std;

class TimeStepControllerCuda {
public:

	// Create a TimeStepController for a certain number of timesteps
	TimeStepControllerCuda();
	~TimeStepControllerCuda();

	/// Add a Timestep, returns the index of the TS in the TSC
	void addTimestep(TimestepCudaPtr ts, int& timestepIndex);

	/// Register the Timestep with index i to be used at a certain point-in-time
	void registerTimestep(int ts_index, double ts_time);

	// Add a new timestep and immediately register it for a certain point-in-time
	void addAndRegisterTimestep(TimestepCudaPtr ts, double ts_time);

	void locateParticlesCuda(float4* points, int numPoints, int* dBlocks);

	void advectParticles(float4* dPos, float timestep, int numParticles);

	void getOccupiedBlocks(int& numBlocks, boost::dynamic_bitset<uchar>& b);

	/// call the init function after all timesteps were loaded
	void initTimesteps();

	// Set time to 0 and reset
	void resetTime();

	void switchIntegration();

	void setIntegrationScheme( IntegrationScheme is ) { integration = is; };

	void printHops( bool b) { m_printHops = b; }
	void printStepsizes( bool b) { m_printStepSizes = b; }

	TetMeshCudaPtr getCurrentCudaMesh() {
		return getCurrentTimestep()->getCudaMesh();
	}
	;
	TetMeshCudaPtr getNextCudaMesh() {
		return getNextTimestep()->getCudaMesh();
	}
	;

	void setATSParams(ATSParams params);
	void setParticleLifetime(int num_steps);

	void printStepSizes( int numParticles);
	void printNumHops( int numParticles, float* hops=0 );
	void setPause( bool pause ) { m_pause = pause; };

	void togglePause()  { m_pause = !m_pause; };
	void toggleLimitLifetime()
	{
		setParticleLifetime( m_limit_lifetime ? 100 : 999999999 );
		m_limit_lifetime = !m_limit_lifetime;
	};

	cudaStream_t getExecutionStream() { return streams[ executionStream ]; };

private:

	void resetTimepoints();
	void incrementTimepoints();

	/// get the current timestep from the timestep controller
	TimestepCudaPtr getCurrentTimestep() {
		return timesteps[currentTimestep];
	}
	TimestepCudaPtr getNextTimestep() {
		return timesteps[nextTimestep];
	}
	TimestepCudaPtr getOvernextTimestep() {
		return timesteps[overnextTimestep];
	}

	void _init();
	void _finalize();

	void resetStepSizes();
	void resetOccupiedBlocks();

	void calculateAdvection(float4* dPos, int numParticles, float timestep,
			int ts_a, int ts_b, int sc_a, int sc_b);

	// The current integration scheme that is used for the particle advection
	IntegrationScheme integration;

	// The parameters for the Dopri-5 Adaptive Timestepping
	ATSParams m_ats_params;

	/// The current time
	double time;

	/// The current \delta t
	double _h;

	vistime_t* vistime;

	/// The Maximum Time, time will be reset to 0 if it reaches timeMax
	double timeMax;

	// The Array holding the different time steps
	vector<TimestepCudaPtr> timesteps;

	/// the current, next and overnext timestep
	int currentTimestep, nextTimestep, overnextTimestep;

	// index of start-cell arrays for current timesteps
	int currScIdx;
	int nextScIdx;

	// References to the loaded timesteps
	list<TimePoint> timepoints;
	list<TimePoint>::iterator currentTp;
	list<TimePoint>::iterator nextTp;
	list<TimePoint>::iterator overnextTp;

	// device pointers
	float4* dVelocities;
	int* dStartCells[3];
	uchar* dOccupiedBlocks;
	float* dStepSizes;

	// two cuda-streams are alternated to perform asynchronous data transfer to the GPU
	cudaStream_t streams[2];
	int transferStream;
	int executionStream;

	// some flags..
	bool m_printHops;
	bool m_printStepSizes;
	bool m_pause;
	bool m_limit_lifetime;

};
typedef boost::shared_ptr<TimeStepControllerCuda> TimeStepControllerCudaPtr;

#endif
