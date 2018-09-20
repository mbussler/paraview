#pragma once

#ifndef PARTICLES_GL_H
#define PARTICLES_GL_H

#include "Typedefs.h"
#include "common.cuh"
#include "ParticleImplementationCuda.h"

#include <boost/shared_ptr.hpp>
#include <boost/dynamic_bitset.hpp>
#include <vector>

using namespace std;

class ParticlesGL : public ParticleImplementationCuda
{
public:
	ParticlesGL();
	~ParticlesGL();

	void drawParticles();
	void drawParticleSpheres( float radius );
	void drawParticlesTrace();

	void drawSeeder();
	void drawReferenceGeometry();

    void setParticleRadius(float r) { m_particleRadius = r; }
    void setFOV(float fov) { m_fov = fov; }
    void setWindowSize(int w, int h) { m_window_w = w; m_window_h = h; }

	void switchTraceColor();

protected:

    void _initGL();
	GLuint _compileProgram(const char *vsource, const char *fsource);

	void _drawPoints();

    GLUquadric* m_refGeometry;

    float m_particleRadius;
    float m_fov;
    int m_window_w, m_window_h;

    GLuint m_program;

	int trace_color;
};

typedef boost::shared_ptr<ParticlesGL> ParticlesGLPtr;

#endif
