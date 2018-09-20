
#include "KdTreeCuda.h"

KdTreeCuda::KdTreeCuda()
{
	hS=0;
	hD=0;
	hL=0;
	//hI=0;
};

KdTreeCuda::~KdTreeCuda()
{
   //freeGPUMem();
};

void KdTreeCuda::freeGPUMem()
{
	freeArray( m_kdTreeGPU.dS );
    freeArray( m_kdTreeGPU.dD );
    freeArray( m_kdTreeGPU.dL );
    //freeArray( m_kdTreeGPU.dI );

	if( hS ) freePageLockedHostMemory(hS);
    if( hD ) freePageLockedHostMemory(hD);
    if( hL ) freePageLockedHostMemory(hL);
    //if( hI ) freePageLockedHostMemory(hI);

	hS=0;
	hD=0;
	hL=0;
	//hI=0;
};

void KdTreeCuda::copyDataToGPU()
{
    // calculate memsize of data
	m_kdTreeGPU.sS = sizeof(float) * point_count;
	m_kdTreeGPU.sD = sizeof(char) * point_count;
	m_kdTreeGPU.sL = sizeof(int) * point_count;
	//m_kdTreeGPU.sI = sizeof(int) * point_count;

    // allocate memory on GPU
    allocateArray( (void**) &m_kdTreeGPU.dS, m_kdTreeGPU.sS );
    allocateArray( (void**) &m_kdTreeGPU.dD, m_kdTreeGPU.sD );
    allocateArray( (void**) &m_kdTreeGPU.dL, m_kdTreeGPU.sL );
    //allocateArray( (void**) &m_kdTreeGPU.dI, m_kdTreeGPU.sI );

    copyArrayToDevice( m_kdTreeGPU.dS, S, m_kdTreeGPU.sS );
    copyArrayToDevice( m_kdTreeGPU.dD, D, m_kdTreeGPU.sD );
    copyArrayToDevice( m_kdTreeGPU.dL, L, m_kdTreeGPU.sL );
    //copyArrayToDevice( m_kdTreeGPU.dI, I, m_kdTreeGPU.sI );

	m_kdTreeGPU.levels = levels;
};

// Asynchronous host to device copy using page-locked memory
void KdTreeCuda::copyDataToGPUAsync(cudaStream_t stream)
{
	KdTreeGPU kd;

    // calculate memsize of data
	kd.sS = sizeof(float) * point_count;
	kd.sD = sizeof(char) * point_count;
	kd.sL = sizeof(int) * point_count;
	//kd.sI = sizeof(int) * point_count;

    // allocate memory on GPU
    allocateArray( (void**) &kd.dS, kd.sS );
    allocateArray( (void**) &kd.dD, kd.sD );
    allocateArray( (void**) &kd.dL, kd.sL );
    //allocateArray( (void**) &kd.dI, kd.sI );

	bool writeCombined = true;
    allocatePageLockedArrayPortable((void**)&hS, kd.sS, writeCombined);
    allocatePageLockedArrayPortable((void**)&hD, kd.sD, writeCombined);
    allocatePageLockedArrayPortable((void**)&hL, kd.sL, writeCombined);
	//allocatePageLockedArray((void**)&hI, kd.sI, writeCombined);

    // allocate page-locked host memory and copy asynchronously to device
	copyArrayToPageLockedHostMemory(hS, this->S, kd.sS);
    copyArrayToPageLockedHostMemory(hD, this->D, kd.sD);
    copyArrayToPageLockedHostMemory(hL, this->L, kd.sL);
	//copyArrayToPageLockedHostMemory(hI, I, kd.sI);

    printf("Copy Kd-Tree: ");
    startTimer(stream);

    copyArrayToDeviceAsync( kd.dS, hS, kd.sS, stream );
    copyArrayToDeviceAsync( kd.dD, hD, kd.sD, stream );
    copyArrayToDeviceAsync( kd.dL, hL, kd.sL, stream );
    //copyArrayToDeviceAsync( kd.dI, hI, kd.sI );

    stopTimer(stream);
    printTimer();
    printf("\n");
    destroyTimer(stream);

	kd.levels = levels;
	m_kdTreeGPU = kd;

}


void KdTreeCuda::traverseTree( float4* dPoints, int numPoints, int* dCells, int* dBlocks, MeshGPU meshGPU, cudaStream_t stream )
{
    traverseTreeCuda( dPoints, numPoints, dCells, dBlocks, m_kdTreeGPU, meshGPU, stream);
};
