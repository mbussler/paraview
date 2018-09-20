
#ifndef PARTICLES_KERNEL_CU
#define PARTICLES_KERNEL_CU

#include "common.cuh"

// A particle is deleted after this number of steps
__device__ __constant__ int max_particle_steps = 99999999;

__device__ void
writeOccupiedFlag( uchar blockOccupied, uchar* occupiedBlocks)
{
	// id of the current warp where warpSize is 32
	// For a block-size of 256, this is in 0..7
	const unsigned int warp_id = threadIdx.x / warpSize;

	// offset for the current block, given as number of warps per block multiplied by the block-id
	// For a block-size of 256, this is in 0, 8, 16, ...
	const unsigned int block_off = ( blockDim.x / warpSize ) * blockIdx.x;

	// write "flag"
	occupiedBlocks[ block_off + warp_id] = blockOccupied;
}

__device__ float4 outOfFieldPosition( unsigned int tid )
{
	//return make_float4(	-1.0f + ((tid % 2000) / 1000.0f ), 1.5f + (tid / 2000)*0.05 , 0.0f, -10.0f);
  return make_float4(-1000.0f, -1000.0f, -1000.0f, -1000.0f  );
}

__global__ void
AdvectParticlesKernel
( 
	float4* positions,
	float4* velocities,
	float h,
	int numParticles,
	uchar* occupiedBlocks
)
{
	const unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;

	// this flag is one, if there are active particles in the current block a 32 particles
	// and 0 otherwise
	__shared__ char blockOccupied;
	blockOccupied = 1;

	if( tid < numParticles )
	{
		// read current particle position from global mem
		float4 pos = positions[ tid ];

		// check if particle is inside the flow field
		if(( pos.w > 0 ) && ( pos.w < max_particle_steps ))
		{
			// advect particle
			pos += h*velocities[tid];

			// increase life counter
			pos.w += 1.0f;

			// write new position back to global mem
			positions[tid] = pos;

			// current block still has active particles, set occupied flag to 0
			blockOccupied = 0;
		}
		else if( pos.w > -10.0f ) 
			//positions[tid].w = -1000.0f;
			positions[tid] = outOfFieldPosition( tid );
	}

	// first thread of each warp writes blockOccupied flag
	if( threadIdx.x % warpSize == 0 )
	{
		writeOccupiedFlag( blockOccupied, occupiedBlocks);
	}
};


#endif
