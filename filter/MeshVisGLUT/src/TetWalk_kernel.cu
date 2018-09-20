

#ifndef TETWALK_KERNEL_CU
#define TETWALK_KERNEL_CU

// USE_TEX
#include "common.cuh"

#include "MeshFunctions.cu"

__global__ void
TetWalkKernel
(
	float4* points,
	int numPoints,
	int* startCells,
	int* blocks,
	char* traversedCells
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
		float4 search_point = points[tid];

		if( search_point.w > 0 ) // test, if position is oof
		{
			int current_cell = startCells[tid];
			int stepCount = 0;

			bool oof = ( current_cell < 0 );
			bool crdsValid = false;

			// natural coodinates
			float4 crds;

			// points indices of current cell
			int4 indices;

			while( !( crdsValid || oof ))
			{
				if( current_cell > 0 ) { traversedCells[ current_cell ] = 1; }

				// read cell node indices from mem
				indices = tex1Dfetch(cells1Tex, current_cell);

				// read physical coordinates from mem
				float4 p1 = tex1Dfetch(nodes1Tex, indices.x);
				float4 p2 = tex1Dfetch(nodes1Tex, indices.y);
				float4 p3 = tex1Dfetch(nodes1Tex, indices.z);
				float4 p4 = tex1Dfetch(nodes1Tex, indices.w);

				// calculate natural coordinates of search point in current cell
				crds = calculateNaturalCoordinates( search_point, p1, p2, p3, p4 );

				int4 ns = tex1Dfetch( neighbors1Tex, current_cell);
				int next_cell = calculateNextCell( crds, ns, current_cell);

				crdsValid = ( next_cell == current_cell || stepCount > tw_max_hops );

				// otherwise, use next cell and iterate
				if( !crdsValid )
				{
					current_cell = next_cell;
					oof = ( current_cell < 0 );
					stepCount++;
				}
			}

			// update startcell index
			startCells[tid] = current_cell;
		}
	}
};


#endif
