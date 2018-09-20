

#ifndef KD_TRAVERSE_KERNEL_CU
#define KD_TRAVERSE_KERNEL_CU

#include "common.cuh"

// store KD-Tree as texture
texture<float, 1, cudaReadModeElementType>  STex;
texture<int,   1, cudaReadModeElementType>  ITex;
texture<int,   1, cudaReadModeElementType>  LTex;

// store mesh nodes as texture
texture<float4, 1, cudaReadModeElementType> nodesTex;

__global__ void
traverseTreeSPKernel
( 
     float4* points,
	 int numPoints,
	 int* cells, 
     int* blocks,
     KdTreeGPU kdTree
)
{
	unsigned int tid = 0;
	unsigned int max_tid = 0;

	if( blocks ) 
	{
		// id of the current warp where warpSize=32
		// For a block-size of 256, this is in 0..7
		const unsigned int warp_id = threadIdx.x / 32;

		// offset for the current block, given as number of warps per block multiplied by the block-id
		// For a block-size of 256, this is in 0, 8, 16, ...
		const unsigned int block_off = ( blockDim.x / 32 ) * blockIdx.x;

		const unsigned int off = block_off + warp_id;

		// blocks[0] holds the number of blocks
		if( off < blocks[0] )
		{
			unsigned int base_id = blocks[off+1];
			
			// read current block from offset-array and add thread-id
			tid = base_id + (threadIdx.x % 32);

			// maximum thread-id
			max_tid = base_id + 32;
		}
	} 
	else 
	{
		tid = blockDim.x * blockIdx.x + threadIdx.x;
		max_tid = numPoints;
	}

	if( tid < max_tid )
	{
        float4 p = points[tid];

		if( p.w > 0 ) // test, if position is oof
		{
			float val;
			char splitDim;
			float splitValue;

			int index=0;

			for( int l=0; l<kdTree.levels; l++ )
			{
				// Get split dimension from global mem
				splitDim = kdTree.dD[index];

				// Get actual split value from texture mem
   				splitValue = tex1Dfetch(STex, index);

				// Choose component from p according to splitdim
	        	val = *(&p.x + splitDim);

				// Calculate successor node
				int succ = (( val < splitValue) ? 1 : 2);

				// Calculate new index
				index = 2*index + succ;
			}
			
			// unshift index
			index = index - (1<<kdTree.levels) +1;

			// write cell-index from texture to global mem
			cells[ tid ] = tex1Dfetch(LTex, index);
		}
    }
}

__device__ int 
nearestNeighborEstimate( float4 p, int levels, char* D )
{
	int l=0;
	int i=0;

	while( l<levels )
	{
		// get current split from kd tree
		int	  splitDim   = D[i];
		float splitValue = tex1Dfetch(STex, i);
    
		// get component of point according to split dimension
		float point_component = *(&p.x + splitDim);

		//if( splitDim == 0 ) val = p.x;
		//if( splitDim == 1 ) val = p.y;
		//if( splitDim == 2 ) val = p.z;

		// calculate successor node offset
		int succ = (( point_component < splitValue) ? 1 : 2);

		i = 2*i+succ;
		l++;
	}

	// calculate offset and return
	return i - (1<<levels) +1;
}

__device__ float
calculateDistance( float4 p1, float4 p2 )
{
	return length( p1-p2 );
}

__device__ float 
Random() 
{ 
	return 0.5; // (rand() / (float) RAND_MAX); 
}

__device__ float 
Random( float lower_bound, float upper_bound)
{
	return Random() * (upper_bound-lower_bound)+lower_bound;
}

__device__ float4
RandomPointOnSphere(float4 p, float r)
{
	float x,y,z,s;
	do {
		x = Random(-1.0f,1.0f);
		y = Random(-1.0f,1.0f);
		s = x*x+y*y;
	} while (s > 1);

	const float t = 2*sqrt(1-s);

	x *= t;
	y *= t;
	z = 2*s - 1;

	return make_float4( r*x+p.x, r*y+p.y, r*z+p.z, 0.0f);
};

__global__ void
traverseTreeKernel
( 
     float4* points,
	 int numPoints,
	 int* cells, 
     int* blocks, 
     KdTreeGPU kdTree,
	 float4* nodes
	 )
{
	unsigned int tid;
	unsigned int max_tid;

	if( blocks ) 
	{
		// read current block from offset-array and add thread-id
		tid = blocks[ blockIdx.x ] + threadIdx.x;
		max_tid = blocks[ blockIdx.x ] + 31;
	} 
	else 
	{
		tid = blockDim.x * blockIdx.x + threadIdx.x;
		max_tid = numPoints;
	}

	if( tid < max_tid )
	{

		float4 p = points[tid];

		if( p.w > 0 ) // test, if position is oof
		{
			// first kd tree traversal to estimate distance
			int nn_off = nearestNeighborEstimate( p, kdTree.levels, kdTree.dD  );

			// get node index
			int nn_idx = tex1Dfetch(ITex, nn_off);

			// get coordinates for node index
			float4 nn = tex1Dfetch(nodesTex, nn_idx);

			// calculate current distance
			float dist = calculateDistance( p, nn );

			// traverse with random points
			for( int i=0; i<10; i++) 
			{
				// calculate next query point
				float4 rnd = RandomPointOnSphere(p, dist/2.0f);

				// get offset to access I and L
				int rnd_off = nearestNeighborEstimate( rnd, kdTree.levels, kdTree.dD );

				// get node index by offset from I
				int rnd_idx = tex1Dfetch(ITex, rnd_off);

				// get coordinates from index
				nn = tex1Dfetch(nodesTex, rnd_idx);

				float rnd_dist = calculateDistance( p, nn );

				if( rnd_dist < dist) 
				{
					// store new best distance
					dist=rnd_dist;
					// store offset of current best nearest neighbor
					nn_off = rnd_off;
				}
			}

			// store cell
			cells[tid] = tex1Dfetch( LTex, nn_off );
		}
	}
};

#endif
