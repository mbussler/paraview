
#include "ParticlesGL.h"
#include "shaders.h"

ParticlesGL::ParticlesGL() 
:
	m_particleRadius(0.04f),
	m_program(0),
	trace_color(0)
{
	m_refGeometry = gluNewQuadric();
	gluQuadricNormals(m_refGeometry, GLU_SMOOTH);

	_initGL();
};

ParticlesGL::~ParticlesGL()
{
};

void ParticlesGL::_initGL()
{
    m_program = _compileProgram(vertexShader, spherePixelShader);

#if !defined(__APPLE__) && !defined(MACOSX)
    glClampColorARB(GL_CLAMP_VERTEX_COLOR_ARB, GL_FALSE);
    glClampColorARB(GL_CLAMP_FRAGMENT_COLOR_ARB, GL_FALSE);
#endif
}


GLuint
ParticlesGL::_compileProgram(const char *vsource, const char *fsource)
{
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

    glShaderSource(vertexShader, 1, &vsource, 0);
    glShaderSource(fragmentShader, 1, &fsource, 0);
    
    glCompileShader(vertexShader);
    glCompileShader(fragmentShader);

    GLuint program = glCreateProgram();

    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);

    glLinkProgram(program);

    // check if program linked
    GLint success = 0;
    glGetProgramiv(program, GL_LINK_STATUS, &success);

    if (!success) {
        char temp[256];
        glGetProgramInfoLog(program, 256, 0, temp);
        printf("Failed to link program:\n%s\n", temp);
        glDeleteProgram(program);
        program = 0;
    }

    return program;
}

void ParticlesGL::_drawPoints()
{

#ifdef NOGL
	float* hPos = new float[ 4*m_numParticles];
	copyArrayFromDevice(hPos, dParticles, m_numParticles * sizeof(float4));
	glBegin( GL_POINTS );
	{
		int k = 0;
		for (int i = 0; i < m_numParticles; ++i)
		{
			glVertex3fv(&hPos[k]);
			k += 4;
		}
	}
	glEnd();
	delete[] hPos;
#else

	// positions vbo
	glBindBuffer(GL_ARRAY_BUFFER, m_posVbo);
	glVertexPointer(3, GL_FLOAT, sizeof(float4), BUFFER_OFFSET(0));
	glEnableClientState(GL_VERTEX_ARRAY);

	// color vbo
    glBindBuffer(GL_ARRAY_BUFFER, m_colorVBO);
    glColorPointer(4, GL_FLOAT, 0, 0);
    glEnableClientState(GL_COLOR_ARRAY);

	glDrawArrays(GL_POINTS, 0, m_numParticles );

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

#endif
};


void ParticlesGL::drawParticles()
{
	glPushAttrib(GL_DEPTH_BUFFER_BIT
		| GL_ENABLE_BIT
		| GL_VIEWPORT_BIT
		| GL_POLYGON_BIT ) ;

	glDisable( GL_LIGHTING );

	glColor3f( 0.6f, 0.8f, 1.0f);
	
	glPointSize( 2.5f);
	
	_drawPoints();

	glPopAttrib();
};

void ParticlesGL::drawParticleSpheres( float radius)
{
	glPushAttrib(GL_DEPTH_BUFFER_BIT
		| GL_ENABLE_BIT
		| GL_VIEWPORT_BIT
		| GL_POLYGON_BIT ) ;

	glEnable(GL_POINT_SPRITE_ARB);
	glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_NV);
	glDepthMask(GL_TRUE);
	glEnable(GL_DEPTH_TEST);

	glUseProgram(m_program);
	glUniform1f( glGetUniformLocation(m_program, "pointScale"), m_window_h / tanf(m_fov*0.5f*(float)M_PI/180.0f) );
	glUniform1f( glGetUniformLocation(m_program, "pointRadius"), radius );

	glColor3f(1, 1, 1);

	_drawPoints();

	glUseProgram(0);
	glDisable(GL_POINT_SPRITE_ARB);

	glPopAttrib();
};

void ParticlesGL::drawParticlesTrace()
{
	glPushAttrib(GL_DEPTH_BUFFER_BIT
		| GL_ENABLE_BIT
		| GL_VIEWPORT_BIT
		| GL_POLYGON_BIT ) ;

	glDisable( GL_LIGHTING );

	drawColor( trace_colors[trace_color] );
	glLineWidth(2.0f);

	// positions vbo
	glBindBuffer(GL_ARRAY_BUFFER, m_posVbo);
	glVertexPointer(3, GL_FLOAT, sizeof(float4), BUFFER_OFFSET(0));
	glEnableClientState(GL_VERTEX_ARRAY);

	glDrawArrays(GL_LINE_STRIP, 0, m_numParticles );

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glDisableClientState(GL_VERTEX_ARRAY);

	glPopAttrib();
};

void ParticlesGL::switchTraceColor()
{
	++trace_color %= trace_colors_count;
};


void ParticlesGL::drawSeeder()
{
	if( m_seeder)
	{
		glPushAttrib(GL_DEPTH_BUFFER_BIT
			| GL_ENABLE_BIT
			| GL_VIEWPORT_BIT
			| GL_POLYGON_BIT ) ;
		glDisable(GL_LIGHTING);
		m_seeder->draw();
		glPopAttrib();
	}
};

void ParticlesGL::drawReferenceGeometry()
{
	glPushAttrib(GL_DEPTH_BUFFER_BIT
		| GL_ENABLE_BIT
		| GL_VIEWPORT_BIT
		| GL_POLYGON_BIT );

	glEnable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);

	glColor3f(1.0, 1.0, 1.0);

	glPushMatrix();

	// draw closed cylinder
	glTranslatef(0.0, 0.0, -0.2);
	gluQuadricOrientation(m_refGeometry, GLU_OUTSIDE);
	gluCylinder(m_refGeometry, 0.2, 0.05, 0.4, 32, 32);

	gluQuadricOrientation(m_refGeometry, GLU_INSIDE);
	gluDisk( m_refGeometry, 0.0, 0.2, 32, 1);

	glTranslatef(0.0, 0.0, 0.4);
	gluQuadricOrientation(m_refGeometry, GLU_OUTSIDE);
	gluDisk( m_refGeometry, 0.0, 0.05, 32, 1);

	// draw Sphere
	//	gluSphere(m_refGeometry, 0.1, 32, 32);

	glPopMatrix();

	glPopAttrib();
};

