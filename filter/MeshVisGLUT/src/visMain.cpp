/**
 *	\brief Main class
 *
 *	\author Michael Bu�ler
 */

#ifdef _MEM_LEAK_DETECTION
#define VLD_AGGREGATE_DUPLICATES
#define VLD_MAX_DATA_DUMP 10
#include "vld.h" //visual leak detection
#endif

#include "common.cuh"

#ifdef TETLIBRARY
	#include "TetgenLoader.h"
#endif

#include <stdio.h>
#include <iostream>
#include <vector>

#include "VTKLoader.h"
#include "TimestepCuda.h"
#include "TimeStepControllerCuda.h"

#include "ParticlesGL.h"
#include "MeshVis.h"

// global definition required for GLUT callback
MeshVisPtr vis(new MeshVis());

int main(int argc, char* argv[]) 
{
	srand(time(NULL));

	if (vis->Init( argc, argv))
	{
		cudaGLInit( argc, argv );
		TimeStepControllerCudaPtr TSC(new TimeStepControllerCuda());

		// testing parameters
		float stepsize = 0.02f;
		int count = 1000;
		int num_iterations = 50;

		// data loading parameters
		int start_timestep = 10;
		int end_timestep = 41;
		float time_gap = 0.2f;

		if( checkCmdLineFlag(argc, (const char**)argv, "stepsize"))
			stepsize = getCmdLineArgumentFloat(argc, (const char**)argv, "stepsize");

		ATSParams params;
		params.tol = 10e-3f;
		params.h_min = 0.005f;
		params.rho = 0.9f;
		params.eta = 2.0f;

		if( checkCmdLineFlag(argc, (const char**)argv, "tolerance"))
			params.tol = getCmdLineArgumentFloat(argc, (const char**)argv, "tolerance");
		if( checkCmdLineFlag(argc, (const char**)argv, "hmin"))
			params.h_min = getCmdLineArgumentFloat(argc, (const char**)argv, "hmin" );

		float time = 0.0f;
		int tsi;

		// load .vtk files
		VTKLoader vtk;
        string meshFile = "data/tecplot_0044_3.8e+02.vtu";
        //string meshFile = "/media/visushome/data/Roehrle/test/test.vtu";

		time_gap = 10.0f;

		// set the desired flow parameters
		for (int i=0; i<3; i++)
		{
			/* create Timestep */
			TimestepCudaPtr ts(new TimestepCuda());

			vtk.loadFromFile(meshFile.c_str(), ts->getTetMesh());

			// register Timesteps within the TSC
			int timestepIndex;
			TSC->addTimestep(ts, timestepIndex);
			TSC->registerTimestep( timestepIndex, time);

			time += time_gap;
		}

		// initialize Timesteps
		TSC->initTimesteps();
		TSC->setATSParams( params );
		TSC->setIntegrationScheme( RK3 );

		ParticlesGLPtr particles(new ParticlesGL());
		particles->setTimeStepController(TSC);
		particles->set_timestep( stepsize );

		vis->setTimeStepController( TSC );
		vis->setParticles( particles );
		vis->runMeshVis();
	}

	return 0;
}
;

// rendering callbacks
void display() {
	vis->display();
}
;
void keyboard(unsigned char key, int x, int y) {
	vis->keyboard(key, x, y);
}
;
void keyboard_special(int key, int x, int y) {
	vis->keyboard_special(key, x, y);
}
;
void mouse(int button, int state, int x, int y) {
	vis->mouse(button, state, x, y);
}
;
void motion(int x, int y) {
	vis->motion(x, y);
}
;
void reshape(int w, int h) {
	vis->reshape(w, h);
}
;

void timerFunction(int value) {
	vis->Do();
	if (vis->isRecording())
		glutTimerFunc(1000 / 15.0, timerFunction, 1);
	else
		glutTimerFunc(1, timerFunction, 1);
}
;
void register_callbacks() {
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(keyboard_special);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutReshapeFunc(reshape);
	glutTimerFunc(1, timerFunction, 1);
}
;
