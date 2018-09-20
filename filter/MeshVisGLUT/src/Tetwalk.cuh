
#ifndef TETWALK_CUH
#define TETWALK_CUH

// includes, project
#include "common.cuh"

extern "C"
{
    void PerformTetWalkCuda
    ( 
    // input: device pointer (DP) to search points, number of points, DP to cell indices to start walk at 
        float4* dPoints, int numPoints, int* startCells, int* dBlocks,
        MeshGPU mesh, cudaStream_t stream
    );

	void BindTextureCuda( int textureRef, MeshGPU mesh );
	void UnbindTextureCuda( int textureRef );

	void SetTraversedCellsCuda( int meshIdx, char* dTraversedCells );

}

#endif
