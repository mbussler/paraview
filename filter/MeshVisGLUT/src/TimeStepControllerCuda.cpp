#include "TimeStepControllerCuda.h"

//#include <boost/thread.hpp>
//#include <boost/bind.hpp>

#include <string.h>

TimeStepControllerCuda::TimeStepControllerCuda()
{
	time=0.0;
	timeMax=0.0;
	currentTimestep=0;
	nextTimestep=0;
	overnextTimestep=0;
	integration=Dopri5_ATS;

	transferStream=0;
	executionStream=1;

	currScIdx = 0;

	m_ats_params.tol = 10e-2f;
	m_ats_params.h_min = 0.001f;
	m_ats_params.rho = 0.9f;
	m_ats_params.eta = 2.0f;

	m_printHops = false;
	m_printStepSizes = false;
	m_pause = false;
	m_limit_lifetime = false;

};

void TimeStepControllerCuda::addTimestep( TimestepCudaPtr ts, int& timestepIndex )
{
	timesteps.push_back(ts);
	timestepIndex = currentTimestep;
	currentTimestep++;
};

/// Register a Timestep for a certain point-in-time
void TimeStepControllerCuda::registerTimestep( int ts_index, double ts_time  )
{
	TimePoint tp;
	tp.index = ts_index;
	tp.time = ts_time;
	timepoints.push_back(tp);
	timepoints.sort(TimePointSortPredicate);

	if( ts_time > timeMax )
		timeMax = ts_time;
};

void TimeStepControllerCuda::addAndRegisterTimestep(TimestepCudaPtr ts, double time)
{
	int tsi;
	addTimestep(ts, tsi);
	registerTimestep(tsi, time);
}

TimeStepControllerCuda::~TimeStepControllerCuda()
{
	_finalize();
};

void TimeStepControllerCuda::_init()
{
	size_t size_vel = sizeof(float4) * MAX_PARTICLES;
	allocateArray( (void**)&dVelocities, size_vel);
	setArray( dVelocities, 0, size_vel );

	size_t size_idx = sizeof(int) * MAX_PARTICLES;
	allocateArray( (void**)&dStartCells[0], size_idx);
	allocateArray( (void**)&dStartCells[1], size_idx);
	allocateArray( (void**)&dStartCells[2], size_idx);

	setArray( dStartCells[0], 0, size_idx);
	setArray( dStartCells[1], 0, size_idx);
	setArray( dStartCells[2], 0, size_idx);

	allocateArray( (void**)&dOccupiedBlocks, (MAX_PARTICLES * sizeof(uchar)/32) );
	allocateArray( (void**)&dStepSizes, MAX_PARTICLES * sizeof( float) );

	resetStepSizes();
	resetOccupiedBlocks();

	setATSParameters( &m_ats_params );

	createStreams(2, streams);

	allocatePageLockedArray((void**)&vistime, sizeof( vistime_t), false );
};

void TimeStepControllerCuda::resetStepSizes()
{
	setArray( dStepSizes, 0, MAX_PARTICLES * sizeof( float));
}

void TimeStepControllerCuda::resetOccupiedBlocks()
{
	setArray( dOccupiedBlocks, 0, (MAX_PARTICLES * sizeof(uchar)/32) );
}

void TimeStepControllerCuda::setATSParams(ATSParams params)
{
	m_ats_params = params;
	setATSParameters( &m_ats_params );
};

void TimeStepControllerCuda::setParticleLifetime(int num_steps)
{
	setMaxParticleSteps( num_steps );
}

void TimeStepControllerCuda::_finalize()
{
	freeArray( dVelocities );

	freeArray( dStartCells[0]);
	freeArray( dStartCells[1]);
	freeArray( dStartCells[2]);

	freeArray( dOccupiedBlocks );
	freeArray( dStepSizes );

	freePageLockedHostMemory(vistime);
};

void TimeStepControllerCuda::initTimesteps()
{
	_init();

	for( int i=0; i< timesteps.size(); i++)
	{
		timesteps[i]->init(i);
	}

	resetTimepoints();
};

void TimeStepControllerCuda::resetTime()
{
	// reset time to the beginning 
	time = 0.0f;

	// reset timepoints
	resetTimepoints();
};

void TimeStepControllerCuda::resetTimepoints()
{

	transferStream = 0;
	executionStream = 1;

	// set indices for current, next and overnext timestep
	currentTp = timepoints.begin();
	currentTimestep = currentTp->index;

	nextTp = currentTp;
	safeIncrement(timepoints, nextTp);
	nextTimestep = nextTp->index;

	overnextTp = nextTp;
	safeIncrement(timepoints, overnextTp);
	overnextTimestep = overnextTp->index;

#ifdef GPU_STREAMING
	// copy data to GPU for first three time steps
	getCurrentTimestep()->copyDataToGPU();
	getNextTimestep()->copyDataToGPU();
	getOvernextTimestep()->copyDataToGPUAsync(streams[transferStream]);
#endif

	// set start-cell field for current ts to 0
	currScIdx = 0;
	nextScIdx = 1;
};

void TimeStepControllerCuda::incrementTimepoints()
{
	// switch transfer stream in order to asure asynchronous copy of former increment operations has completed
	executionStream = transferStream;
	++transferStream %= 2;

#ifdef GPU_STREAMING
	// free gpu-mem of the current timestep
	getCurrentTimestep()->freeGPUMem();
#endif

	// increment all indices
	currentTp = nextTp;
	currentTimestep = nextTimestep;

	nextTp = overnextTp;
	nextTimestep = overnextTimestep;

	safeIncrement( timepoints, overnextTp );
	overnextTimestep = overnextTp->index;

	currScIdx = nextScIdx;
	++nextScIdx %= 3;

#ifdef GPU_STREAMING
	// pre-fetch gpu data for the over-next timestep
	getOvernextTimestep()->copyDataToGPUAsync( streams[transferStream] );
#endif

};

void TimeStepControllerCuda::locateParticlesCuda( float4* points, int numPoints, int* dBlocks )
{
	//getNextTimestep()->getCudaMesh()->resetTraversedCells();

	// locate particles in current and next ts + update cell indices
	getCurrentTimestep()->locateParticlesCuda( points, numPoints, dStartCells[ currScIdx], dBlocks, streams[executionStream] );
	getNextTimestep()->locateParticlesCuda( points, numPoints, dStartCells[ nextScIdx], dBlocks, streams[executionStream] );

	// HOLODEMO: reset traversed cells after point location, looks ugly otherwise
	//getCurrentTimestep()->getCudaMesh()->resetTraversedCells();
};

void TimeStepControllerCuda::advectParticles( float4* dPos, float h, int numParticles )
{

	if( time+h < nextTp->time )
	{
		// calculate advection for current and next timestep
		calculateAdvection( dPos, numParticles, h, currentTimestep, nextTimestep, currScIdx, nextScIdx );

		if( !m_pause ) time += h;
	}
	else 
	{
		// calculate Advection for the first half of the timestep
		double h1 = nextTp->time - time;
		double h2 = h-h1;

		if( h1 > 0.1 * h)
			calculateAdvection( dPos, numParticles, (float)h1, currentTimestep, nextTimestep, currScIdx, nextScIdx );

		time += h1;

		// check if time has not yet crossed the maximum time barrier
		if( time < timeMax )
		{
			// increment references and copy new data to GPU
			incrementTimepoints();

			// locate all particles in the new timestep
			getNextTimestep()->getCudaMesh()->resetTraversedCells();
			getNextTimestep()->locateParticlesCuda( dPos, numParticles, dStartCells[(currScIdx+1)%3], 0, streams[executionStream] );

			// perform advection for the rest of the timestep

			// dismiss very small steps
			if( h2 > 0.1 * h)
				calculateAdvection( dPos, numParticles, h2, currentTimestep, nextTimestep, currScIdx, nextScIdx );
			//resetStepSizes();
		}
		else
		{
			// time has crossed the maximum time -> reset everything
			#ifdef GPU_STREAMING
				getCurrentTimestep()->freeGPUMem();
				getNextTimestep()->freeGPUMem();
				getOvernextTimestep()->freeGPUMem();
			#endif

			resetTimepoints();
			time = time - timeMax;

			// locate particles
			getCurrentTimestep()->locateParticlesCuda( dPos, numParticles, dStartCells[0], 0, streams[executionStream] );
			getNextTimestep()->locateParticlesCuda( dPos, numParticles, dStartCells[1], 0, streams[executionStream] );

			// perform advection for the rest of the timestep

			// dismiss very small steps
			if( h2 > 0.1 * h)
				calculateAdvection( dPos, numParticles, h2, currentTimestep, nextTimestep, currScIdx, nextScIdx );
			//resetStepSizes();
		}

		time += h2;
	}

};

void TimeStepControllerCuda::calculateAdvection( float4* dPos, int numParticles, float timestep, int ts_a, int ts_b, int sc_a, int sc_b )
{
	// The normalized time for the interpolation
	double a = (time - currentTp->time) / (nextTp->time - currentTp->time);

	// The normalized timestep
	double ah = timestep / (nextTp->time - currentTp->time);

	// Clear traversed cells
	timesteps[ts_a]->getCudaMesh()->resetTraversedCells();
	timesteps[ts_b]->getCudaMesh()->resetTraversedCells();

	MeshGPU mesh1 = timesteps[ts_a]->getCudaMesh()->mapVBOtoCuda();
	MeshGPU mesh2 = timesteps[ts_b]->getCudaMesh()->mapVBOtoCuda();

	vistime->t = (float) time;
	vistime->h = (float) timestep;
	vistime->a = (float) a;
	vistime->ah = (float) ah;

	IntegrateVelocitiesCuda( integration, dPos, dVelocities, numParticles, 
		vistime,
		mesh1, mesh2,
		dStartCells[sc_a], dStartCells[sc_b],
		dStepSizes, streams[executionStream] );

	timesteps[ts_a]->getCudaMesh()->unmapVBOtoCuda();
	timesteps[ts_b]->getCudaMesh()->unmapVBOtoCuda();

	AdvectParticlesCuda( dPos, dVelocities, numParticles, timestep, dOccupiedBlocks, streams[executionStream] );

	if( m_printHops ) 
	{
		printNumHops( numParticles );
		printStepSizes( numParticles );
	}
};

void TimeStepControllerCuda::getOccupiedBlocks( int& numBlocks, boost::dynamic_bitset<uchar>& b )
{
	// create host array for blocks and copy from gpu
	uchar* blocks;
	size_t size = numBlocks * sizeof(uchar);
	allocatePageLockedArray( (void**)&blocks, size );

	copyArrayFromDeviceAsync( blocks, dOccupiedBlocks, size, getExecutionStream() );

	// prepare bit array
	b.resize( numBlocks, false);

	// lastblock store index of the last occupied block
	int lastBlock = numBlocks;

	// assure async transfer has finished
	cudaStreamSynchronize(getExecutionStream());

	// iterate over host array and copy values to bit array
	for( int i=0; i<numBlocks; i++ )
	{
		b[i] = blocks[i];
		if( b[i] == 0 ) lastBlock = i+1;
	}

	// set new value for number of blocks according to the last occupied block
	numBlocks = lastBlock;

	// clean up
	if( size > 0)
	{
		freePageLockedHostMemory(blocks);
	}
};

void TimeStepControllerCuda::switchIntegration()
{
	switch( integration)
	{
	case Euler:
		integration = RK3;
		printf("switched to RK-3 Integration Scheme.\n");
		break;
	case RK3:
		integration = RK4;
		printf("switched to RK-4 Integration Scheme.\n");
		break;
	case RK4:
		integration = Dopri5;
		printf("switched to Dopri-5 Integration Scheme.\n");
		break;
	case Dopri5:
		integration = Dopri5_ATS;
		printf("switched to Dopri-5 with adaptive timestepping.\n");
		break;
	case Dopri5_ATS:
		integration = Euler;
		printf("switched to Euler Integration Scheme.\n");
		break;
	}
};

void TimeStepControllerCuda::printStepSizes(int numParticles )
{
	float* hStepSize;
	size_t size = numParticles * sizeof(float);
	allocatePageLockedArray((void**)&hStepSize, size);

	copyArrayFromDevice( hStepSize, dStepSizes, numParticles*sizeof(float));

	float sum = 0.0f;
	float min = 9999.0f;
	float max = 0.0f;
	int c = 0;

	for( int i=0; i<numParticles; i++ )
	{
		float stepsize = hStepSize[i];

		if( stepsize  > 0 )
		{
			sum += stepsize;
			min = ( stepsize < min ) ? stepsize : min;
			max = ( stepsize > max ) ? stepsize : max;
			c++;
		}
	}

	float avg = sum / (float) c;

	//printf("Current Step sizes:\n");
	//printf("min: %.4f max: %.4f avg: %.4f\n", min, max, avg);
	printf(", %.4f,%.4f, %.4f\n", min, max, avg);

	freePageLockedHostMemory(hStepSize);
};

void TimeStepControllerCuda::printNumHops( int numParticles, float* hops )
{
	size_t size = numParticles * sizeof(int);
	int* hNumHops1;
	int* hNumHops2;
	int* hIterations;

	// allocate host memory to store the number of traversed cells to
	allocatePageLockedArray((void**)&hNumHops1, size);
	allocatePageLockedArray((void**)&hNumHops2, size);
	allocatePageLockedArray((void**)&hIterations, size);

	// download values from the GPU
	getNumberOfHops( numParticles, hNumHops1, hNumHops2, hIterations );

	// initialize the result values 
	float resSumHops = 0.0f;
	float resMinHops = 1000000.0f;
	float resMaxHops = 0.0f;

	int resSumIterations = 0;
	int resMaxIterations = 0;
	int resMinIterations = 1000000;

	// number of "valid" results
	int c = 0;

	// create a histogram of the different number of hops
	int histogram[20] = {0};

	vector<float> vHops;
	vector<int> iters;

	for (int i = 0; i < numParticles; i++)
	{
		int numHops1 = hNumHops1[i];
		int numHops2 = hNumHops2[i];
		int numIterations = hIterations[i];

		if( numIterations > 0 )
		{
			iters.push_back( numIterations );

			resSumIterations += numIterations;
			resMinIterations = std::min<int>( resMinIterations, numIterations);
			resMaxIterations = std::max<int>( resMaxIterations, numIterations);

			float numHops = ((numHops1+numHops2)/2.0f) / (float)numIterations;
			vHops.push_back( numHops );

			resSumHops += numHops;
			resMinHops = std::min<int>( resMinHops, numHops);
			resMaxHops = std::max<int>( resMaxHops, numHops);

			if( numHops > 0 && numHops < 19 ) { histogram[(int)(numHops+0.5)]++; }

			c++;
		}
	}

	float avg_hops = resSumHops / c;
	float avg_iters = resSumIterations / (float) c;

	vector<float>::iterator medHops = vHops.begin() + (vHops.end() - vHops.begin())/2;
	vector<int>::iterator medIters = iters.begin() + (iters.end() - iters.begin())/2;

	nth_element(vHops.begin(), medHops, vHops.end());
	nth_element(iters.begin(), medIters, iters.end());

	if( hops )
	{
		hops[0] = resMinHops;
		hops[1] = resMaxHops;
		hops[2] = avg_hops;
	}
	else
	{
		printf("%d, %.3f,%.3f,%.3f,%.3f, %d,%d,%d,%.3f", integration, resMinHops, resMaxHops, *medHops, avg_hops, resMinIterations, resMaxIterations, *medIters, avg_iters );
	}

	//	for( int i=0; i<20; i++)
	//	{
	//		printf(",%d", histogram[i]);
	//	}
	//
	//	printf("\n");

	freePageLockedHostMemory(hNumHops1);
	freePageLockedHostMemory(hNumHops2);
	freePageLockedHostMemory(hIterations);
};
