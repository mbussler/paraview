
#include "KdTreeTest.h"

#include "common.cuh"
#include "KdTraverse.cuh"

#include <vector_types.h>
#include <boost/shared_ptr.hpp>

#ifndef KDTREE_CUDA_H
#define KDTREE_CUDA_H

class KdTreeCuda : public KdTreeTest
{
public:
    KdTreeCuda();
    ~KdTreeCuda();

    // Asynchronous host to device copy using page-locked memory
    void copyDataToGPUAsync(cudaStream_t stream);

    void copyDataToGPU();
    void freeGPUMem();

    /// Traverse Kd-Tree with CUDA
	/*!
	 *  [in] dPoints device pointer to search points
     *  [in] num_points number of search points
     *  [out] dIndices device pointer to store found cell indices at
	 */    
	void traverseTree( float4* dPoints, int numPoints, int* dCells, int* dBlocks, MeshGPU meshGPU, cudaStream_t stream );

private:
	KdTreeGPU m_kdTreeGPU;

	// host pointers to page-locked memory
	float* hS;
	char*  hD;
	int*   hL;
	//int* hI;
};

typedef boost::shared_ptr<KdTreeCuda> KdTreeCudaPtr;

#endif

