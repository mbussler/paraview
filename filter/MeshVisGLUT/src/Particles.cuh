
#ifndef PARTICLES_CUH
#define PARTICLES_CUH

// includes
#include "common.cuh"

extern "C"
{
	void setATSParameters( ATSParams *hostParams );
	void setMaxParticleSteps( int num_steps );
	
	void AdvectParticlesCuda( float4* dPositions, float4* dVelocities, int numParticles, float timestep, uchar* occupiedBlocks, cudaStream_t stream);
    
	void IntegrateVelocitiesCuda
    (
        IntegrationScheme integration,
        float4* dPositions, float4* dVelocities, int numParticles,
		vistime_t* vistime,
		MeshGPU mesh1, MeshGPU mesh2,
        int* dStartCells1, int* dStartCells2,
		float* dStepSizes, cudaStream_t stream
	);

	void getNumberOfHops( int count, int* hNumHops1, int* hNumHops2, int* hIterations );
}

#endif
