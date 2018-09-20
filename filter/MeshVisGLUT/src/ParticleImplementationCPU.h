#pragma once

#ifndef PARTICLE_IMPLEMENTATION_GPU_H
#define PARTICLE_IMPLEMENTATION_GPU_H

#include "TimeStepControllerCPU.h"
#include "Particle.h"
#include "ParticleSeeder.h"

#include "Typedefs.h"

#include <boost/shared_ptr.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>

class ParticleImplementationCPU
{
public:
	ParticleImplementationCPU();
	~ParticleImplementationCPU();

	void createParticles();
	void createParticlesTest(int count);
	void createParticlesVBO();
	void createTracesVBO();
	void updateParticlesVBO();
	void ParticleStep();
	void clearParticles();
	void deleteParticleVBO();

	void drawParticles();
	void drawParticleSpheres( float radius = 0.1f);
	void drawTeapots( float radius = 0.1f);
	void drawTraces();
	void drawSeeder();

	void drawReferenceGeometry();

	void switchSeeder();

	int getTotalParticleCount() { return m_particleList.size(); };
	int getActiveParticleCount() { return m_activeParticles; };

	void setTimeStepController( TimeStepControllerCPUPtr tsc)
		{ m_tsc = tsc;}

	const float timestep() const { return m_timestep; };
	void set_timestep( const float timestep) 
	{ 
		m_timestep = timestep; 
		printf("timestep: %4.2f\n", timestep);
	};

	const bool clear_particle() const { return m_clearParticles; }
	void set_clear_particle( const bool clear_particle ) { m_clearParticles = clear_particle; }

	ParticleSeeder* seeder() { if (!m_seeder) m_seeder = new ParticleSeeder(); return m_seeder; };

protected:

	// The Particles are Stored in a ParticleList
	ParticleList			m_particleList;

	// This flag indicates whether particles, the leave the flow field are deleted
	bool m_clearParticles;

	// Iterators to traverse the ParticleList and PointLists
	ParticleListIterator	p_i;
	PointListIterator		pl_i;

	// The TimeStepController is used to calculate the velocity for the particles
	TimeStepControllerCPUPtr 	m_tsc;

	// The timestep by which the particles are moved
	float 					m_timestep;// = 0.1f;

	// The Particle Seeder
	ParticleSeeder* 		m_seeder; //=0;

	// OpenGL Vertex Buffer Objects
	GLuint					vbo_traces;
	GLuint					vbo_particles;

	// array of offsets to draw the traces from the traces-VBO
	vector<int>				m_traces_offset;
	int						m_traces_count;

	GLUquadricObj 			*m_refGeometry;

	int m_activeParticles;

};
typedef boost::shared_ptr<ParticleImplementationCPU> ParticleImplementationCPUPtr;


#endif
