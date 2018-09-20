#pragma once

#ifndef PARTICLE_IMPLEMENTATION_CUDA_H
#define PARTICLE_IMPLEMENTATION_CUDA_H

#include "TimeStepControllerCuda.h"
#include "ParticleSeeder.h"

#include "Typedefs.h"
#include "common.cuh"
#include "Particles.cuh"

#include <boost/shared_ptr.hpp>
#include <boost/dynamic_bitset.hpp>

#include <vector>

using namespace std;

#define		ParticlesPerBlock	32

class ParticleImplementationCuda
{
public:
	ParticleImplementationCuda();
	~ParticleImplementationCuda();

	void createParticles();
	void createSingleParticle();
	void createParticlesTest( int count );
    void ParticleStep();
	void clearParticles();

    void switchSeeder();

    int getTotalParticleCount() { return m_numParticles; };
	int getActiveParticleCount() { return m_activeParticles; };

	void setTimeStepController( TimeStepControllerCudaPtr tsc)
		{ m_tsc = tsc;}

	const float timestep() const { return m_timestep; };
	void set_timestep( const float timestep) 
	{ 
		m_timestep = timestep; 
		printf("timestep: %f\n", timestep);
	};

	const bool clear_particle() const { return false; }
	void set_clear_particle( const bool clear_particle ) { ; }

	ParticleSeeder* seeder() { if (!m_seeder) m_seeder = new SphereParticleSeeder(); return m_seeder; };

protected:

    void _init();
    void _finalize();

	void updatePositionsVBO( float4* newParticles, int count, int* hOffsets=0);
	void updateColorVBO( int count, int* hOffsets=0);

	void calculateFreeBlocks( int* hOffsets, int numBlocks, int& totalBlocks );

    unsigned int m_numParticles;
    unsigned int m_activeParticles;

    unsigned int m_particleCounter;

    uint createVBO(uint size);

    uint   m_posVbo;        // vertex buffer object for particle positions
    uint   m_colorVBO;      // vertex buffer object for colors

    struct cudaGraphicsResource *m_cuda_posvbo_resource; // handles OpenGL-CUDA exchange
    
	float4* dParticles;

	// The TimeStepController is used to advect the particles
	TimeStepControllerCudaPtr 	m_tsc;

	// The timestep by which the particles are moved
	float 					m_timestep;// = 0.1f;

	// The Particle Seeder
	ParticleSeeder* 		m_seeder; //=0;

    GLUquadric* m_refGeometry;

    void dumpParticles(uint start, uint count);

	int particleArrayPos;

};
typedef boost::shared_ptr<ParticleImplementationCuda> ParticleImplementationCudaPtr;


#endif
