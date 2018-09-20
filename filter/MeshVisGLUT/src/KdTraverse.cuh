
#ifndef KD_TRAVERSE_CUH
#define KD_TRAVERSE_CUH

#include <vector_types.h>

// number of levels to traverse in the kdtree
extern __constant__ int cLevels;

extern "C"
{
	void traverseTreeCuda
    ( 
    // input: (initialized) device pointer to searchpoints, number of points
        float4* dPoints, int numPoints, int* dCells,
        int* dBlocks,

    // input: Kd-Tree on GPU
		KdTreeGPU kdTree,

    // input: Mesh on GPU
        MeshGPU mesh,

        cudaStream_t stream
	);
}

#endif
