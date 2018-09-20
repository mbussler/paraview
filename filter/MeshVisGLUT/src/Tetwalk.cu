
#ifndef TETWALK_CU
#define TETWALK_CU

#include "Tetwalk.cuh"

// includes, kernels
#include "TetWalk_kernel.cu"

extern "C"
{

    void PerformTetWalkCuda
    ( 
        float4* dPoints, int numPoints, int* startCells,  int* dBlocks,
        MeshGPU mesh,
        cudaStream_t stream
    )
    {
		if( numPoints < 1 ) return; // dont allow empty blocks

		checkCudaErrors(cudaBindTexture(0, nodes1Tex, mesh.nodes, mesh.s_nodes));
		checkCudaErrors(cudaBindTexture(0, cells1Tex, mesh.cells, mesh.s_cells));
		checkCudaErrors(cudaBindTexture(0, neighbors1Tex, mesh.neighbors, mesh.s_cells));

		// setup execution parameters
		int threadsPerBlock = 256; // number of threads per block
		int numBlocks = ( numPoints + threadsPerBlock-1) / threadsPerBlock;

		CREATE_TIMER
		START_TIMER(stream)

		//// execute the kernel
		TetWalkKernel<<< numBlocks, threadsPerBlock, 0, stream >>> ((float4*)dPoints, numPoints, startCells, dBlocks, mesh.traversedCells );
	    
		getLastCudaError("Kernel execution failed");

        STOP_TIMER(stream)
        PRINT_TIMER("TetWalk: ", "\n");
		STORE_TIMER(tetwalkTimer)
        DESTROY_TIMER


		checkCudaErrors(cudaUnbindTexture(nodes1Tex));
		checkCudaErrors(cudaUnbindTexture(cells1Tex));
		checkCudaErrors(cudaUnbindTexture(neighbors1Tex));

	};

}

#endif
