#pragma once

#ifndef MESHVIS_H
#define MESHVIS_H

/**
 *
 *	\brief Mesh Visualization
 *  
 *	\author Michael Bu�ler
 */

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <cstdio>
#include <string.h>
#include <math.h>

// includes, GL
#include <GL/glew.h>
#include <GL/freeglut.h>
#define BUFFER_OFFSET(bytes) ((GLubyte*) NULL + (bytes))

#include "InlineFunctions.h"
#include "BoundingBox.h"
#include <list>
#include "Camera.h"

#include "ParticlesGL.h"
//#include "ParticleImplementationCPU.h"

// fps counter
#include "DisplayFPS.h"
#include "TetMeshGL.h"

// recodr
#include "RecordingTools.h"

const GLfloat LightAmbient[]= { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat LightDiffuse[]= { 1.0f, 0.576f, 0.172f, 1.0f };
const GLfloat LightPosition[]= { 5.0f, -5.0f, -5.0f};

const GLfloat HeadLightAmbient[]= { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat HeadLightDiffuse[]= { 1.0f, 0.576f, 0.172f, 1.0f };
const GLfloat HeadLightPosition[] = {0.0f, 0.0f, -2.0f};

// Testing parameters for framecount
static bool run_test = false;

// particle populations for testing
static int test_particles[] = {10000, 20000, 50000, 100000, 200000, 500000};

// Number of particle counts
static int test_numTests = 6;
static int test_currentTest = 2;
static int test_startAt = 2;


// Number of repetitions per particle count
static int test_numIterations = 1;
static int test_currentIteration = 0;

// Number of steps per particle count
static int test_numSteps = 300;
static int test_currentStep = 0;

void register_callbacks();

class MeshVis 
{

public:
	MeshVis();
	~MeshVis();
	
	
	bool Init( int argc, char* argv[] );
	bool Do();    

	void runMeshVis();

	void setTimeStepController(TimeStepControllerCudaPtr tsc) { m_tsc = tsc;}
	//void setTimeStepController(TimeStepControllerCPUPtr tsc) { m_tsc = tsc;}

	void setParticles( ParticlesGLPtr particles)
	{
		m_particles = particles;
		m_particles->setFOV( 45.0 );
		m_particles->setWindowSize(window_width, window_height);
	}
	//void setParticles( ParticleImplementationCPUPtr particles) { m_particles = particles; }

	// Recorder
	RecordingTools rec;
	bool isRecording() {return rec.isRecording();};

public:

	// rendering callbacks
	void display();
	void keyboard(unsigned char key, int x, int y);
	void keyboard_special(int key, int x, int y);
	void mouse(int button, int state, int x, int y);
	void motion(int x, int y);
	void reshape(int w, int h);

private:

	bool initGLUT( int argc, char* argv[]);
	void setup_scene();       //< init scene
	void setup_interaction(); //< init interaction
	void setup_navigation();  //< init scene navigation

	void draw_scene(); //< draw scene

	void render(); //< draw stuff
	void update(); //< update scene

	void setSeederPosition(int x, int y);

	// reference to the TimeStepController to draw the mesh etc
	TimeStepControllerCudaPtr		m_tsc;
	//TimeStepControllerCPUPtr		m_tsc;

	// reference to the ParticleImplementation to draw the Particles etc
	ParticlesGLPtr				m_particles;
	//ParticleImplementationCPUPtr	m_particles;

	bool m_createParticles;

	bool m_pause;

	int m_last_time;   //< GLUT_ELAPSED_TIME is stored here
	int m_frame_time;  //< also GLUT_ELAPSED_TIME, but relative
	int m_update_time;  //< accumulate frame time

	float m_stepsPerSecond;

	bool m_update_scene;

	// Camera
	Camera vis_cam;
	float move_speed, rotate_x, rotate_y;
	Vector translate;

	// Glut Window dimensions
	unsigned int window_width;
	unsigned int window_height;

	float nearplane, farplane;

	// mouse controls
	int mouse_old_x, mouse_old_y;
	int mouse_buttons;
	int glut_modifiers;

	// flags
	bool render_free_cam;
	
	// render state control
	bool m_render_billboards;
	bool m_drawBox;
	bool m_drawVertices;
	bool m_drawVelocities;
	bool m_drawCells;
	bool m_drawTraversedCells;
	bool m_drawKdTreewalk;
	bool m_drawTraces;
    bool m_drawOuterFaces;

	float m_shereRadius;

};
typedef boost::shared_ptr<MeshVis> MeshVisPtr;


#endif
