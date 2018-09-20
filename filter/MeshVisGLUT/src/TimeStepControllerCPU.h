#pragma once

#ifndef TIMESTEP_CONTROLLER_CPU_H
#define TIMESTEP_CONTROLLER_CPU_H

#include <boost/shared_ptr.hpp>

#include "TimestepCPU.h"
#include "Typedefs.h"
#include "Vector.h"
#include "TetMesh.h"
#include <list>
#include <vector>

#include "common.cuh"
#include "timer.h"

using namespace std;

class TimeStepControllerCPU
{
public:

	// Create a TimeStepController for a certain number of timesteps
	TimeStepControllerCPU();
	~TimeStepControllerCPU();

	/// Add a Timestep, returns the index of the TS in the TSC
	void addTimestep( TimestepCPUPtr ts, int& timestepIndex );

	/// Register the Timestep with index i to be used at a certain point-in-time
	void registerTimestep( int ts_index, double ts_time );

	/// get a timestep from the timestep controller
	TimestepCPUPtr getTimestep( int timestepIndex );

	/// get the current timestep from the timestep controller
	TimestepCPUPtr getCurrentTimestep();

	/// get the next timestep from the timestep controller
	TimestepCPUPtr getNextTimestep();

	TimestepCPUPtr getOvernextTimestep();

	TetMeshGLPtr getCurrentMeshGL();
	TetMeshGLPtr getNextMeshGL();

	/// call the init function after all timesteps were loaded
	void initTimesteps();

	  // Set time to 0 and reset
	  void resetTime();

	void advectParticles( ParticleList& particles, float timestep );

	int timestepCount;

	void switchIntegration();

	void setATSParams( ATSParams params )
	{
		m_ats_params = params;
	};

  double getTime() { return timer1.time(); };
  
private:

	void calculateAdvection( ParticleList& particles, float timestep, int ts_a, int ts_b );

	void CalculateIntegration(		Particle& p, Vector& v, float h, float a, float ah, int ts_a, int ts_b);
	void CalculateEulerIntegration(	Particle& p, Vector& v, float h, float a, float ah, int ts_a, int ts_b);
	void CalculateRK3Integration(	Particle& p, Vector& v, float h, float a, float ah, int ts_a, int ts_b);
	void CalculateRK4Integration(	Particle& p, Vector& v, float h, float a, float ah, int ts_a, int ts_b);
	void CalculateDopri5Integration( Particle& p, Vector& v_rk4, Vector& v_rk5, float h, float a, float ah, int ts_a, int ts_b);
	void CalculateDopri5Integration_ATS ( Particle& p, Vector& v, float timestep, float a, float ah, int ts_a, int ts_b	);
	/// The Array holding the different time steps
	vector<TimestepCPUPtr> timesteps;

	list<TimePoint> timepoints;
	list<TimePoint>::iterator currentTp;
	list<TimePoint>::iterator nextTp;
	list<TimePoint>::iterator overnextTp;

	/// the current and next timestep
	int currentTimestep, nextTimestep, overnextTimestep;

	void resetTimepoints();
	void incrementTimepoints();

	/// The current time 
	double time;

	/// The current \delta t
	double _h;

	/// The Maximum Time, time will be reset to 0 if it reaches timeMax
	double timeMax;

	ATSParams m_ats_params;

	IntegrationScheme integration;

  timer timer1;

};
typedef boost::shared_ptr<TimeStepControllerCPU> TimeStepControllerCPUPtr;


#endif
