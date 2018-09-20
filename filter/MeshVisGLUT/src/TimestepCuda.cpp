#include "TimestepCuda.h"

//#include <boost/thread.hpp>

TimestepCuda::TimestepCuda()
{
	m_mesh.reset( new TetMeshCuda() );
	bDataCopied = false;
};

TimestepCuda::~TimestepCuda()
{
};

void TimestepCuda::init(int timestepIndex)
{
	printf("\n");
	printf("Initializing [Timestep %d]:\n", timestepIndex);

	//m_mesh->init();

	printf("Processing Kd-Tree..\n");
		m_kdtree.reset( new KdTreeCuda() );
		m_kdtree->createKdTree( m_mesh );

#ifndef GPU_STREAMING
	printf("Transferring Data to GPU..\n");
	copyDataToGPU();
#endif

	printf("done.\n");
};

void TimestepCuda::locateParticlesCuda( float4* dPoints, int numPoints, int* dStartCells, int* dBlocks, cudaStream_t stream)
{
	// get reference to mesh on GPU
	MeshGPU meshGPU = m_mesh->mapVBOtoCuda();

	// localize, global
	m_kdtree->traverseTree( dPoints, numPoints, dStartCells, dBlocks, meshGPU, stream );

	// release reference
	m_mesh->unmapVBOtoCuda();

    // localize, local
    m_mesh->PerformTetWalk( dPoints, numPoints, dStartCells, dBlocks, stream );

};

void TimestepCuda::copyDataToGPU()
{
	if (!bDataCopied)
	{
		m_mesh->copyDataToGPU();
		m_kdtree->copyDataToGPU();

		bDataCopied = true;
	}
};


void TimestepCuda::copyDataToGPUAsync(cudaStream_t stream)
{
#ifdef ASYNC_TRANSFER
	if (!bDataCopied)
	{
		m_kdtree->copyDataToGPUAsync(stream);
		//boost::thread copyTreeThread(&KdTreeCuda::copyDataToGPUAsync, m_kdtree, stream );

		m_mesh->copyDataToGPUAsync(stream);
		//boost::thread copyMeshThread(&TetMeshCuda::copyDataToGPUAsync, m_mesh, stream );

		bDataCopied = true;
	}
#else
	copyDataToGPU();
#endif
};

void TimestepCuda::freeGPUMem()
{
	if( bDataCopied )
	{
		m_mesh->freeGPUMem();
		m_kdtree->freeGPUMem();

		bDataCopied = false;
	}
};
