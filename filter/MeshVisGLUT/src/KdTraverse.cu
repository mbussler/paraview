
#ifndef KD_TRAVERSE_CU
#define KD_TRAVERSE_CU

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include <helper_cuda.h>

// includes, kernels
#include "common.cuh" //USE_TEX
#include "KdTraverse_kernel.cu"


////////////////////////////////////////////////////////////////////////////////
//! Entry point for Cuda functionality on host side
//! @param argc  command line argument count
//! @param argv  command line arguments
//! @param data  data to process on the device
//! @param len   len of \a data
////////////////////////////////////////////////////////////////////////////////

extern "C" 
{
	void traverseTreeCuda
    ( 
        float4* dPoints, int numPoints, int* dCells,
        int* dBlocks,
		KdTreeGPU kdTree,
        MeshGPU mesh,
        cudaStream_t stream
	)
    {

		if( numPoints < 1 ) return;

		// Bind textures
		checkCudaErrors(cudaBindTexture(0, STex, kdTree.dS, kdTree.sS));
		checkCudaErrors(cudaBindTexture(0, LTex, kdTree.dL, kdTree.sL));
		//cutilSafeCall(cudaBindTexture(0, ITex, kdTree.dI, kdTree.sI));
		//cutilSafeCall(cudaBindTexture(0, nodesTex, mesh.nodes, mesh.s_nodes));

        // setup execution parameters
        int threadsPerBlock = 256; // Evaluate 128 positions per block
        int numBlocks = ( numPoints + threadsPerBlock-1 ) / threadsPerBlock;

        CREATE_TIMER
        START_TIMER(stream)

        //// execute the kernel
        traverseTreeSPKernel<<< numBlocks, threadsPerBlock, 0, stream >>> ((float4*) dPoints, numPoints, dCells, dBlocks, kdTree );
        
        getLastCudaError("Kernel execution failed");
        
        STOP_TIMER(stream)
        PRINT_TIMER("Tree-Traverse: ", "\n");
        STORE_TIMER(kdTreeTimer)
        DESTROY_TIMER

        checkCudaErrors(cudaUnbindTexture(STex));
        checkCudaErrors(cudaUnbindTexture(LTex));
        //cutilSafeCall(cudaUnbindTexture(ITex));
        //cutilSafeCall(cudaUnbindTexture(nodesTex));

    };
}
#endif



