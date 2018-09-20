
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

#include <stdio.h>
#include <iostream>
#include <vector>

#include "TetgenLoader.h"
#include "VTKLoader.h"
#include "TimestepCuda.h"
#include "TimeStepControllerCuda.h"

#include "ParticlesGL.h"
#include "KdTree.h"

#include "MeshVis.h"

MeshVisPtr vis( new MeshVis() );

float flow_param1[5][3] = { { -0.2, 3.0, -1.5 }, { -0.6, 6.0, -1.5 }, { -1.2,
		9.0, -1.5 }, { -0.6, 6.0, -1.5 }, { -0.2, 3.0, -1.5 } };

float flow_param2[5][3] = { { 0.0, 1.0, 0.1 }, { 0.0, 1.0, 0.1 }, { 0.0, 1.0,
		0.1 }, { 0.0, 1.0, 0.1 }, { 0.0, 1.0, 0.1 } };

float flow_param3[5][3] = { { 0.0, 1.0, -0.3 }, { 0.0, 1.0, 0.0 }, { 0.0, 1.0,
		0.3 }, { 0.0, 1.0, 0.0 }, { 0.0, 1.0, -0.3 } };

const BoundingBox Domain = BoundingBox(Point(-1.0, -1.0, -1.0), Point(1.0, 1.0,
		1.0));

int main(int argc, char* argv[]) {
	srand(time(NULL));

	if (vis->Init( argc, argv))
	{
		cudaGLInit( argc, argv );
		//cudaInit( argc, argv );
		TimeStepControllerCudaPtr TSC(new TimeStepControllerCuda());

		// testing parameters
		float stepsize = 0.02f;
		int count = 100000;
		int num_iterations = 50;
		int num_timesteps = 5;

		cutGetCmdLineArgumentf(argc, (const char**)argv, "stepsize", &stepsize);
		cutGetCmdLineArgumenti(argc, (const char**)argv, "count", &count);
		cutGetCmdLineArgumenti(argc, (const char**)argv, "iterations", &num_iterations);

		ATSParams params;
		params.tol = 10e-4f;
		params.h_min = 0.001f;
		params.rho = 0.9f;
		params.eta = 2.0f;

		cutGetCmdLineArgumentf(argc, (const char**)argv, "tolerance", &params.tol);
		cutGetCmdLineArgumentf(argc, (const char**)argv, "hmin", &params.h_min);
		bool printHops = cutCheckCmdLineFlag(argc, (const char**)argv, "printHops");

		TSC->printHops( printHops );

		// load meshes with tetgen
		TetgenLoaderPtr loader( new TetgenLoader() );
		loader->setDomain( Domain );

		// load .vtk files
		VTKLoader vtk;
		char meshFile[255];
		int fileIdx=20;

		float time = 0.0f;
		int tsi;

		// set the desired flow parameters
		for (int i = 0; i <= num_timesteps; i++)
		{
			/* create Timestep */
			TimestepCudaPtr ts(new TimestepCuda());

			/* create a Tetrahedral Grid with the given flow parameters */
			//loader->setSynthFlowParameter( flow_param3[1] );
			//loader->createTetrahedralGrid(ts->getTetMesh(), 64, 64, 64, 0.01f, true);

			/* load an existing mesh with tetgen. */
			//sprintf(meshFile, "TetMesh_20x20x20_a%4.2f_b%4.2f_c%4.2f", flow_param3[1][0], flow_param3[1][1], flow_param3[1][2]);
			//loader->loadMeshFromFiles( meshFile, ts->getTetMesh());

			// load engine??.vtk-files with the VTKLoader:
		#ifdef _WIN32
			sprintf(meshFile, "d:\\Engine\\engine_%02d.vtk", fileIdx);
		#else
			sprintf(meshFile, "./Engine/engine_%02d.vtk", fileIdx);
		#endif

			vtk.loadFromFile(meshFile, ts->getTetMesh());
			fileIdx += 1;

			// register Timesteps within the TSC
			TSC->addTimestep(ts, tsi);
			TSC->registerTimestep( tsi, time);

			time += (num_iterations * stepsize) / (float) num_timesteps;
		}

		// initialize Timesteps
		TSC->initTimesteps();
		ParticlesGLPtr particles( new ParticlesGL());
		particles->setTimeStepController(TSC);

		TSC->setATSParams( params );
		particles->set_timestep( stepsize );

		//for( int i=0; i<5; i++ )
		//{
		//	particles->createParticlesTest(100000);
		//	TSC->switchIntegration();

		//	for( int j=0; j<100; j++ )
		//	{
		//		particles->ParticleStep();

		//		if( j%10 == 0 ) 
		//		{
		//			getGPUMemoryUsage( &mem_free, &mem_total );
		//			std::cout << particles->getTotalParticleCount() << "\t";
		//			std::cout << "Memory free/total: " << mem_free << "/" << mem_total  << std::endl;
		//		}
		//	}
		//}

		vis->setTimeStepController( TSC );
		vis->setParticles( particles );
		vis->runMeshVis();

		return 0;
	}
};


// rendering callbacks
void display()
{ vis->display(); };
void keyboard(unsigned char key, int x, int y)
{ vis->keyboard( key, x, y); };
void keyboard_special(int key, int x, int y)
{ vis->keyboard_special(key, x, y); };
void mouse(int button, int state, int x, int y)
{ vis->mouse(button, state, x, y); };
void motion(int x, int y)
{ vis->motion(x, y); };
void reshape(int w, int h)
{ vis->reshape(w, h); };

void timerFunction( int value )
{
	vis->Do();
	if( vis->isRecording() )
		glutTimerFunc( 1000/15.0, timerFunction, 1);
	else
		glutTimerFunc( 1, timerFunction, 1);
};
void register_callbacks()
{
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(keyboard_special);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutReshapeFunc(reshape);
	glutTimerFunc( 1, timerFunction, 1);
};
