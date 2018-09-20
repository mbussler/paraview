
#pragma once

#ifndef TIMESTEP_H
#define TIMESTEP_H

#include "common.cuh"

#include "TetMeshCuda.h"
#include "KdTreeCuda.h"

/**
 *	\brief Timestep class
 *  This class represents one timestep, including the mesh and the search structure
 *
 *	\author Michael Bu�ler
*/

class TimestepCuda
{
public:
	TimestepCuda();
	~TimestepCuda();

	//void setMesh( TetMeshGL* mesh );
	void init(int timestepIndex);

    void copyDataToGPU();
    void freeGPUMem();

    void copyDataToGPUAsync(cudaStream_t stream);

	/// locate cells for a set of particles and store for each Particle.
	/*!
	 *  This function needs to be called once to locate the particles cells
	 *  The Cell indices for this Timestep are stored in the particles Cell field given by timestepIndex
	 */
	void locateParticlesCuda( float4* dPoints, int numPoints, int* dStartCells, int* dBlocks, cudaStream_t stream=0);

	TetMeshCudaPtr getCudaMesh() { return m_mesh; };
	TetMeshPtr   getTetMesh() { return m_mesh; };
	
	char* filename;

private:

	TetMeshCudaPtr  m_mesh;
	KdTreeCudaPtr   m_kdtree;

	bool bDataCopied;

};

typedef boost::shared_ptr<TimestepCuda> TimestepCudaPtr;

#endif
