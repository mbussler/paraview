
#pragma once

#ifndef TIMESTEP_CPU_H
#define TIMESTEP_CPU_H


#include "TetMeshGL.h"
#include "KdTreeTest.h"

#include "TetrahedronCreator.h"
#include "TetWalker.h"
#include "Particle.h"

/**
 *	\brief TimestepCPU class
 *  This class represents one timestep, including the mesh and the search structure
 *
 *	\author Michael Bu�ler
*/

class TimestepCPU
{
public:
	TimestepCPU();
	~TimestepCPU();

	//void setMesh( TetMeshGL* mesh );
	void init(int _timestepIndex);

	/// locate cells for a set of particles and store for each Particle.
	/*!
	 *  This function needs to be called once to locate the particles cells
	 *  The Cell indices for this Timestep are stored in the particles Cell field given by timestepIndex
	 */
	void locateParticles( ParticleList& pl );
	void locateParticle( Particle& p );

	void interpolateVelocityAt( const Point& position, int& startCell, Vector& velocity);

	TetMeshPtr getTetMesh()
	{ 
		return m_meshgl; 
	};
	TetMeshGLPtr getMeshGL()
	{ 
		return m_meshgl; 
	};
	TetWalkerPtr getTetWalker()		
	{ 
		return m_tetWalker; 
	};

	void storeTraversedCells( bool store );
	
	double pointInTime;
	int	timestepIndex;

	char* filename;

private:

	TetMeshGLPtr  m_meshgl;
	TetWalkerPtr  m_tetWalker;
	KdTreeTestPtr  m_kdtree;
	TetrahedronCreatorPtr  m_tetCreator;

};

typedef boost::shared_ptr<TimestepCPU> TimestepCPUPtr;

#endif
