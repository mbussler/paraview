
#include "MeshVis.h"
#include "glh_linear.h"

MeshVis::MeshVis()
{
};

MeshVis::~MeshVis()
{
};

void MeshVis::runMeshVis()
{
    // create VBO
	glutMainLoop();
};


bool MeshVis::Init( int argc, char* argv[] )
{
	window_width = 1024;
	window_height = 768;

	if( initGLUT( argc, argv ) )
	{
		// setup_navigation
		setup_navigation();

		// setup scene
		setup_scene();

		// setup interaction
		setup_interaction();

		m_last_time = glutGet(GLUT_ELAPSED_TIME);
		m_update_time = 0;

		return true;
	}

	return false;
};

bool MeshVis::Do()
{
	// update timing and animation
	update();

	// drawing stuff
	render();

	return true;

};

void MeshVis::draw_scene()
{

    // save current lighting and meterial state of the OpenGL state machine
	glPushAttrib(GL_LIGHTING );
	glEnable(GL_LIGHTING);


	if( m_render_billboards)
	{
		m_particles->drawParticleSpheres(m_shereRadius);
	}
	else
	{
		m_particles->drawParticles();
	}

	if( m_drawTraces)
	{
		m_particles->drawParticlesTrace();
	}

	m_particles->drawSeeder();
	//m_particles->drawReferenceGeometry();
	
	//SetupRC();
	//RenderJet();

	DWORD renderFlags = 0;
	renderFlags |= //DRAW_BOX |
			(m_drawOuterFaces   ? DRAW_OUTER_FACES : 0) |
			(m_drawCells 		? DRAW_CELLS : 0) |
			(m_drawVertices 	? DRAW_VERTICES : 0) |
			(m_drawVelocities 	? DRAW_VELOCITIES : 0) |
			(m_drawTraversedCells ? DRAW_TRAVERSED_CELLS : 0);

	m_tsc->getCurrentCudaMesh()->draw( renderFlags );
	//m_tsc->getCurrentMeshGL()->draw( renderFlags );
//
//	if( m_drawKdTreewalk )
//	{
//		m_tsc->getNextCudaMesh()->draw( DRAW_TRAVERSED_CELLS );
//	}

    // restore lighting state (which has been pushed via glPushAttrib)
	glPopAttrib();
};

void MeshVis::setup_navigation()
{
	// flags
	render_free_cam	= false;
	vis_cam.Position.z = -3.0f;
	move_speed = 0.05f;

	rotate_x = 0.0;
	rotate_y = 0.0;
	translate = Vector(0.0, -1.0, -3.0);

	// mouse controls
	mouse_buttons = 0;
};

void MeshVis::setup_interaction()
{
	m_createParticles = true;
	m_update_scene = true;
	m_pause = false;

	m_stepsPerSecond = 30.0f;
};

void MeshVis::setup_scene()
{
    // Set Viewport
    glViewport(0, 0, window_width, window_height);

    // Use Depth Test
    glEnable(GL_DEPTH_TEST);

    // Set Projection
	nearplane = 0.1; farplane = 30.0;
    glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
    gluPerspective(45.0, (GLfloat)window_width / (GLfloat) window_height, nearplane, farplane);

	m_drawVelocities	= false;
	m_drawVertices		= false;
	m_drawBox			= true;
	m_drawCells			= false;
	m_drawTraversedCells = false;
	m_drawKdTreewalk	= false;
	m_drawTraces		= false;
	m_drawOuterFaces	= true;

	m_render_billboards = true;
	m_shereRadius = .01;

	/// setup light
	glEnable( GL_LIGHT0);
	glLightfv(GL_LIGHT0,GL_DIFFUSE,LightDiffuse);
	glLightfv(GL_LIGHT0,GL_AMBIENT,LightAmbient);
	glLightfv(GL_LIGHT0,GL_POSITION,LightPosition);

	glEnable( GL_LIGHT1);
	glLightfv(GL_LIGHT1,GL_DIFFUSE,HeadLightDiffuse);
	glLightfv(GL_LIGHT1,GL_AMBIENT,HeadLightAmbient);
	glLightfv(GL_LIGHT1,GL_POSITION,HeadLightPosition);
};


void MeshVis::update()
{
	// update timing
	const int time = glutGet(GLUT_ELAPSED_TIME);
	
	m_frame_time = time - m_last_time;
	m_last_time = time;

	m_update_time += m_frame_time;

	if( m_update_time > 1000.0f / m_stepsPerSecond ) // GLUT_ELAPSED_TIME is in ms!!!
	{
		m_update_time = 0;

		if( run_test )
		{
			if( test_currentTest < test_numTests )
			{
				if( test_currentStep == 0)
				{
					// initialize test settings
					m_update_scene = false;
					m_createParticles = false;
					m_tsc->resetTime();
					m_particles->clearParticles();
					FPS::resetRecorder();
					FPS::record = true;

					m_particles->seeder()->center = Point(0.0f, 0.7f, 0.0f);
					m_particles->createParticlesTest( test_particles[test_currentTest] );
					test_currentStep++;
				}
				else
				{
					m_particles->ParticleStep();
					test_currentStep++;

					if( test_currentStep > test_numSteps )
					{
						printf("NumParticles:%d ", test_particles[test_currentTest]);
						FPS::printMedianFPS();
						printf("\n");

						test_currentStep = 0;
						test_currentIteration++;

						if( test_currentIteration >= test_numIterations)
						{
							test_currentIteration = 0;
							test_currentTest++;
						}
					}
				}
			}
			else
			{
				run_test = false;
				test_currentTest = test_startAt;
				FPS::record = false;
			}
		}
		else
		{
			if( m_update_scene )
			{
				// Do a Particle Step
				m_particles->ParticleStep();
			}

			if( m_createParticles )
			{
				m_particles->createParticles();

				//m_tsc->resetTime();
				//m_particles->clearParticles();
				//m_particles->seeder()->center = Point(0.0f, 0.7f, 0.0f);;
				//m_particles->createParticlesTest(100000);
				//m_createParticles = false;
				//m_update_scene = true;
			}

			if( m_drawTraversedCells)
			{
				m_tsc->getCurrentCudaMesh()->createTraversedCellsVBO();
				//m_tsc->getCurrentMeshGL()->createTraversedCellsVBO();
			}

			if( m_drawKdTreewalk )
			{
				m_tsc->getNextCudaMesh()->createTraversedCellsVBO();
				//m_tsc->getNextMeshGL()->createTraversedCellsVBO();
			}

			//if( m_drawTraces )
			//{
			//	m_particles->createTracesVBO();
			//}
		}

		FPS::iterations++;
		FPS::totalParticles = m_particles->getTotalParticleCount();
		FPS::timestep = m_particles->timestep();

		getGPUMemoryUsage(&FPS::mem_free, &FPS::mem_total, 1 << 20);
		m_tsc->printNumHops( m_particles->getTotalParticleCount(), FPS::hops );
	}
}

void MeshVis::render()
{
	// clear background
	glClearColor(0.4f, 0.4f, 0.4f, 1.0f);
	//glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// set view matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	if( render_free_cam) 
	{
		gluLookAt(vis_cam.Position.x, vis_cam.Position.y  , vis_cam.Position.z ,
				  vis_cam.LookAt.x  , vis_cam.LookAt.y    , vis_cam.LookAt.z   ,
				  vis_cam.Up.x      , vis_cam.Up.y        , vis_cam.Up.z       );
	} 
	else 
	{
		glTranslatef( translate.x, translate.y, translate.z);
	    glRotatef(rotate_x, 1.0, 0.0, 0.0);
		glRotatef(rotate_y, 0.0, 1.0, 0.0);
	}

	GLfloat m[16]; 
	glGetFloatv (GL_MODELVIEW_MATRIX, m); 
	glh::matrix4 mv;
	mv.set_value(m);

    // Update Light0 position according to current MV-Transformation
    glLightfv(GL_LIGHT0, GL_POSITION, LightPosition);

    // Draw the Scene
	draw_scene();

	// Increment frame counter
	FPS::counter++;

    if( rec.isRecording())
	{
        bool push_state = m_update_scene;
        m_update_scene = false;

        rec.saveFrame();

        m_update_scene = push_state;
	}
	else
	{
		FPS::displayFps();
	}

	glutSwapBuffers();
	glutPostRedisplay();
};

bool MeshVis::initGLUT( int argc, char* argv[])
{
	// Create GL context
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(window_width, window_height);
	glutCreateWindow("MeshVis");

	// initialize necessary OpenGL extensions
    glewInit();
    if (! glewIsSupported("GL_VERSION_2_0 ")) {
        fprintf(stderr, "ERROR: Support for necessary OpenGL extensions missing.");
        fflush(stderr);
        return false;
    }

	// register GLUT callbacks
	register_callbacks();

	return true;
};

void MeshVis::display()
{
};

void MeshVis::keyboard(unsigned char key, int x, int y)
{
	bool state;

    switch(key) {
		case '1': TOGGLE(m_drawVertices); break;
		case '2': TOGGLE(m_drawVelocities); break;
		case '3': TOGGLE(m_drawCells); break;
		case '4': TOGGLE(m_drawTraversedCells); break;
		case '5': TOGGLE(m_drawKdTreewalk); break;
		case '6': TOGGLE(m_drawTraces ); break;
		case '7': TOGGLE(m_drawOuterFaces ); break;
		case 'b': TOGGLE(m_render_billboards); break;
		case 'c': m_particles->clearParticles(); break;
		case 'C': m_particles->switchTraceColor(); break;
		case 'z': m_particles->switchSeeder(); break;
		case ' ': TOGGLE(m_update_scene); TOGGLE(m_createParticles); break;
		case 'p': TOGGLE(m_createParticles); break;
		case 'P': TOGGLE(m_pause); m_tsc->setPause(m_pause); break;
		case 'm': m_tsc->switchIntegration(); break;
		case 'n': m_stepsPerSecond+=1.0f; break;
		case 'N': if(m_stepsPerSecond>0) m_stepsPerSecond-=1.0f; break;
		case 'r': rec.toggleRecording(); break;
		case 'f': TOGGLE(render_free_cam); break;
		// free cam
		case 'W': vis_cam.CameraMove( move_speed ); break;
		case 'S': vis_cam.CameraMove( -move_speed ); break;
		case 'A': vis_cam.CameraStrafe( move_speed ); break;
		case 'D': vis_cam.CameraStrafe( -move_speed ); break;
		case 'E': vis_cam.CameraMoveUp( move_speed ); break;
		case 'Q': vis_cam.CameraMoveUp( -move_speed ); break;
		case '+': m_particles->set_timestep( m_particles->timestep() * 1.1); break;
		case '-': m_particles->set_timestep( m_particles->timestep() * 0.9); break;
		case 'o': m_shereRadius+=.001;break;//m_particle_render->set_radius(m_shereRadius);break;
		case 'O': m_shereRadius-=.001;break;//m_particle_render->set_radius(m_shereRadius);break;
		case 'l': m_particles->seeder()->size+=.001;break;
		case 'L': m_particles->seeder()->size-=.001;break;
		
		case 'T': TOGGLE(run_test); break;

		// quit with either 'q' or 'ESC'
		case(27) :
		case 'q':
	        exit(0);
    }
};
void MeshVis::keyboard_special(int key, int x, int y)
{
	switch(key) {
		// Move Seeder
		case GLUT_KEY_PAGE_UP:		m_particles->seeder()->center.z += .1f;break;
		case GLUT_KEY_PAGE_DOWN:	m_particles->seeder()->center.z -= .1f; break;
		case GLUT_KEY_LEFT:			m_particles->seeder()->center.x += .02f; break;
		case GLUT_KEY_RIGHT:		m_particles->seeder()->center.x -= .02f; break;
		case GLUT_KEY_UP:			m_particles->seeder()->center.y += .02f; break;
		case GLUT_KEY_DOWN:			m_particles->seeder()->center.y -= .02f; break;
    }
};

void MeshVis::mouse(int button, int state, int x, int y)
{
	mouse_buttons ^= 1<<button;
	glut_modifiers = glutGetModifiers();

    mouse_old_x = x;
    mouse_old_y = y;
};

void MeshVis::motion(int x, int y)
{
    int dx, dy;
    dx = x - mouse_old_x;
    dy = y - mouse_old_y;

	if(  render_free_cam && mouse_buttons == 1 )
	{
		vis_cam.CameraMove( dx, dy );
	}
	else 
	{
		if( mouse_buttons == 1 )
		{
			rotate_x += dy * 0.2;
			rotate_y += dx * 0.2;
		}
		if( mouse_buttons == 4 )
		{
			float step_y = dy * 0.01;
			float step_x = dx * 0.01;

			if ( glut_modifiers == GLUT_ACTIVE_SHIFT)
			{
				translate.y -= step_y;
				translate.x += step_x;
			}
			else
			{
				translate.z += (translate.z + step_y < -0.1) ? step_y : 0;
			}
		}

		if( mouse_buttons == 2 )
		{
			setSeederPosition(x,y);
		}

		if( mouse_buttons == 6 )
		{
			m_particles->seeder()->center.z += 0.01*dy;
		}
	}

    mouse_old_x = x;
    mouse_old_y = y;
};

void MeshVis::setSeederPosition(int x, int y)
{
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];

	glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
	glGetDoublev( GL_PROJECTION_MATRIX, projection );
	glGetIntegerv( GL_VIEWPORT, viewport );

	Point& c = m_particles->seeder()->center;

	// use gluProject to calculate z-coord
	GLdouble winX, winY, winZ;
	gluProject(c.x, c.y, c.z, modelview, projection, viewport, &winX, &winY, &winZ);

	winX = (float)x;
	winY = (float)viewport[3] - (float)y;
	
	GLdouble posX, posY, posZ;
	gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);

	c.x = posX;
	c.y = posY;
	c.z = posZ;

  printf("%.3f, %.3f, %.3f\n", posX, posY, posZ);
};

void MeshVis::reshape(int w, int h)
{
	if (h==0) h=1;

	glMatrixMode( GL_PROJECTION);

	glLoadIdentity();

	FPS::window_width = window_width = w;
	FPS::window_height = window_height = h;

	glViewport(0, 0, w, h);
    gluPerspective(45.0, (GLfloat)window_width / (GLfloat) window_height, nearplane, farplane);

	m_particles->setFOV( 45.0 );
	m_particles->setWindowSize(w, h);
};


