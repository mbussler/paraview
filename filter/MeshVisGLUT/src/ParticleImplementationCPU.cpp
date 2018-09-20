
#include "ParticleImplementationCPU.h"

//#include <boost/thread.hpp>
//#include <boost/bind.hpp>

#include <iostream>

ParticleImplementationCPU::ParticleImplementationCPU()
{
	//m_tsc = 0;
	m_timestep = 0.02f;
	m_seeder = new CubeParticleSeeder();
	m_seeder->center = Point(0.0f, 0.7f, 0.0f);

	m_traces_count = 0;
	m_clearParticles = true;
	
	m_refGeometry = gluNewQuadric();
	gluQuadricNormals(m_refGeometry, GLU_SMOOTH);
};

ParticleImplementationCPU::~ParticleImplementationCPU()
{
#ifndef NOGL
	deleteParticleVBO();
#endif
};


// create some particles
void ParticleImplementationCPU::createParticles()
{
	// set count
	int count=32;

	switch (m_seeder->type)
	{
	default:
	case PointSeederType:
		count = 1;
		break;
	case RingSphereSeederType:
		count = 64;
		break;
	case OnSphereSeederType:
		count = 128;
		break;
	case CubeSeederType:
		count = 512;
		break;
//	case SphereSeederType:
//		break;
	}

	// create new particles
	float4* newParticles = new float4[count];
	m_seeder->createPoints( newParticles, count);

	ParticleList createdParticles;

	for( int i=0; i<count; i++)
	{
		Particle particle( newParticles[i] );
		particle.Cell = new int[ m_tsc->timestepCount ];
		particle.bOutOfField = false;
		createdParticles.push_back( particle );
	}

	m_tsc->getCurrentTimestep()->locateParticles(createdParticles);
	m_tsc->getNextTimestep()->locateParticles(createdParticles);

	m_particleList.insert( m_particleList.end(), createdParticles.begin(), createdParticles.end());

#ifndef NOGL
	createParticlesVBO();
	//createTracesVBO();
#endif

	// clean up
	delete[] newParticles;
};

// create particles for performance testing
void ParticleImplementationCPU::createParticlesTest(int count)
{
	// create new particles
	float4* newParticles = new float4[count];
	m_seeder->createPoints( newParticles, count);

	ParticleList createdParticles;

	for( int i=0; i<count; i++)
	{
		Particle particle( newParticles[i] );
		particle.Cell = new int[ m_tsc->timestepCount ];
		particle.bOutOfField = false;
		createdParticles.push_back( particle );
	}

	m_tsc->getCurrentTimestep()->locateParticles(createdParticles);
	m_tsc->getNextTimestep()->locateParticles(createdParticles);

	m_particleList.insert( m_particleList.end(), createdParticles.begin(), createdParticles.end());

#ifndef NOGL
	createParticlesVBO();
	//createTracesVBO();
#endif

	// clean up
	delete[] newParticles;
}


/* Trace Particles through the flow field */
void ParticleImplementationCPU::ParticleStep()
{

	m_activeParticles = 0;
	// Clear oof particles and refresh active particle counter
	for( p_i = m_particleList.begin(); p_i != m_particleList.end();)
	{
		if( p_i->bOutOfField && m_clearParticles )
		{
			p_i = m_particleList.erase(p_i);
		}
		else
		{
			m_activeParticles++;
			p_i++;
		}
	}

	// advect particles
	m_tsc->advectParticles( m_particleList, m_timestep );

#ifndef NOGL
	updateParticlesVBO();
#endif
};

void ParticleImplementationCPU::clearParticles()
{
	m_particleList.clear();
	m_traces_count = 0;
};

void ParticleImplementationCPU::switchSeeder()
{
	if ( !m_seeder)
		{
			m_seeder = new ParticleSeeder();
			return;
		}

	ParticleSeeder* old_seeder = m_seeder;

	switch (m_seeder->type)
	{
		case PointSeederType:
			m_seeder = new RingParticleSeeder();
		break;
		case RingSphereSeederType:
			m_seeder = new OnSphereParticleSeeder();
		break;
		case OnSphereSeederType:
			m_seeder = new CubeParticleSeeder();
		break;
		case CubeSeederType:
			m_seeder = new ParticleSeeder();
		break;
//		case SphereSeederType:
//			m_seeder = new SphereParticleSeeder());
//		break;
	}

	m_seeder->center = old_seeder->center;
	delete old_seeder;
}

void ParticleImplementationCPU::drawParticleSpheres( float radius)
{
	glPushAttrib(GL_DEPTH_BUFFER_BIT
				| GL_ENABLE_BIT
				| GL_VIEWPORT_BIT
				| GL_POLYGON_BIT );

	glEnable(GL_LIGHTING);

    for( p_i = m_particleList.begin(); p_i != m_particleList.end(); p_i++)
	{
		glPushMatrix();
		glTranslatef(p_i->Position.x, p_i->Position.y, p_i->Position.z);
		glutSolidSphere( radius, 32, 32);
		glPopMatrix();
	}
	glPopAttrib();
};

void ParticleImplementationCPU::drawTeapots( float radius)
{
	glPushAttrib(GL_DEPTH_BUFFER_BIT
				| GL_ENABLE_BIT
				| GL_VIEWPORT_BIT
				| GL_POLYGON_BIT );

	glEnable(GL_LIGHTING);
	glFrontFace(GL_CW);
	for( p_i = m_particleList.begin(); p_i != m_particleList.end(); p_i++)
	{
		glPushMatrix();
		glTranslatef(p_i->Position.x, p_i->Position.y, p_i->Position.z);
		glutSolidTeapot(radius);
		glPopMatrix();
	}
	glFrontFace(GL_CCW);
	glPopAttrib();
};

void ParticleImplementationCPU::drawParticles()
{
	glPushAttrib(GL_DEPTH_BUFFER_BIT
				| GL_ENABLE_BIT
				| GL_VIEWPORT_BIT
				| GL_POLYGON_BIT ) ;

	unsigned int particleCount = m_particleList.size();
	
	glDisable( GL_LIGHTING );
    glEnableClientState(GL_VERTEX_ARRAY);

    glColor3f( 0.6f, 0.8f, 1.0f);
    glPointSize( 2.5f);

	glBindBuffer(GL_ARRAY_BUFFER, vbo_particles);
    glVertexPointer(3, GL_FLOAT, 0, BUFFER_OFFSET(0));

    glDrawArrays(GL_POINTS, 0, particleCount);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glDisableClientState(GL_VERTEX_ARRAY);

	glPopAttrib();
};

void ParticleImplementationCPU::drawSeeder()
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

void ParticleImplementationCPU::drawReferenceGeometry()
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


// VBO Rendering Methods
void ParticleImplementationCPU::createParticlesVBO()
{
	if (!glIsBuffer( vbo_particles ))
		glGenBuffers(1, &vbo_particles );

	int particleCount = m_particleList.size();

	// store 3 floats per particle position
	unsigned int size = (3 * particleCount);

	GLfloat* data = new GLfloat[size];

	int idx = 0;

	for( p_i = m_particleList.begin(); p_i != m_particleList.end(); p_i++)
	{
		Point p = p_i->Position;

		for( int i=0; i<3; i++) {
			data[3*idx+i] = p.c[i];
		}
		idx++;
	}

	glBindBuffer(GL_ARRAY_BUFFER, vbo_particles);
	glBufferData(GL_ARRAY_BUFFER, size * sizeof(GLfloat), data, GL_DYNAMIC_DRAW);

	glBindBuffer( GL_ARRAY_BUFFER, 0 );

	delete data;
};

void ParticleImplementationCPU::drawTraces()
{
//	if (!glIsBuffer( vbo_traces ))
//	{
//		createTracesVBO();
//	}
//
//	glPushAttrib( GL_DEPTH_BUFFER_BIT
//						| GL_ENABLE_BIT
//						| GL_VIEWPORT_BIT
//						| GL_LINE_BIT
//						| GL_POLYGON_BIT ) ;
//
//	glEnable(GL_LIGHTING);
//
//	// load VBOs
//	glBindBuffer(GL_ARRAY_BUFFER, vbo_traces);
//
//	glEnableClientState(GL_VERTEX_ARRAY);
//	glEnableClientState(GL_NORMAL_ARRAY);
//	glEnableClientState(GL_COLOR_ARRAY);
//
//	glVertexPointer(  3, GL_FLOAT, 12*sizeof(GL_FLOAT), BUFFER_OFFSET(0));
//	glNormalPointer(     GL_FLOAT, 12*sizeof(GL_FLOAT), BUFFER_OFFSET(4));
//	glColorPointer(   3, GL_FLOAT, 12*sizeof(GL_FLOAT), BUFFER_OFFSET(8));
//
//	glLineWidth( 4.0 );
//	
//	for( int i=0; i<m_traces_count; i++)
//	{
//		//drawColor( trace_colors[ i%trace_colors_count ] );
//		glDrawArrays( GL_LINE_STRIP, m_traces_offset[i], m_traces_offset[i+1]-m_traces_offset[i]);
//	}
//
//	glDisableClientState(GL_VERTEX_ARRAY);
//	glDisableClientState(GL_NORMAL_ARRAY);
//	glDisableClientState(GL_COLOR_ARRAY);
//
//	glBindBuffer(GL_ARRAY_BUFFER, 0);
//
//	glPopAttrib();
};

void ParticleImplementationCPU::createTracesVBO()
{
//	if (!glIsBuffer( vbo_traces ))
//		glGenBuffers(1, &vbo_traces );
//
//	m_traces_count = m_particleList.size();
//
//	// store a trace offset for each particle
//	m_traces_offset.resize( m_traces_count+1 );
//	m_traces_offset[0] = 0;
//
//	int vertex_count = 0;
//	int idx = 1;
//	for( p_i = m_particleList.begin(); p_i != m_particleList.end(); p_i++)
//	{
//		vertex_count += p_i->Trace.size();
//		m_traces_offset[idx] = vertex_count;
//		idx++;
//	}
//
//	// store all trace vertices in a vbo
//	unsigned int size = ( (4+4+4) * vertex_count );
//
//	GLfloat* data = new GLfloat[size];
//
//	idx=0;
//	
//	// iterate over all particles
//	ParticleListIterator particle = m_particleList.begin();
//	for( ; particle != m_particleList.end(); particle++)
//	{
//		// iterate over the trace
//
//		PointListIterator f = particle->Trace.begin();
//		PointListIterator p = particle->Trace.begin(); // != particle->Trace.end()) ? ++f : f;
//		PointListIterator n = particle->Trace.begin(); // != particle->Trace.end()) ? ++p : p;
//
//		if( particle->Trace.size() < 2)
//		{
//			PointListIterator f = particle->Trace.begin();
//			PointListIterator p = ++f; // != particle->Trace.end()) ? ++f : f;
//			PointListIterator n = ++p; // != particle->Trace.end()) ? ++p : p;
//		}
//		
//		while( n != particle->Trace.end() )
//		{
//			// Write Vertex
//			data[ 12*idx+0 ] = p->x;
//			data[ 12*idx+1 ] = p->y;
//			data[ 12*idx+2 ] = p->z;
//			//data[ 12*idx+3 ] = 1.0;
//
//			// calculate Normal
//			Vector norm;
//
//			if( f!=p && p!=n)
//			{
//				Vector m1 = VectorNormalize(*f - *p);
//				Vector m2 = VectorNormalize(*p - *n);
//				
//				norm = VectorCross(norm, m1);
//			}
//
//			data[ 12*idx+4] = ( norm.x);  // -dy
//			data[ 12*idx+5] = ( norm.y);  //  dx
//			data[ 12*idx+6] = ( norm.z);  // dz (is this correct?)
//			//data[ 12*idx+7 ] = 1.0;
//
//			// Write Color
//			data[ 12*idx+8 ] = trace_colors[idx%trace_colors_count].r;
//			data[ 12*idx+9 ] = trace_colors[idx%trace_colors_count].g;
//			data[ 12*idx+10] = trace_colors[idx%trace_colors_count].b;
//			//data[ 12*idx+11 ] = .6;
//
//			f = p;
//			p = n;
//			n = ++n;
//			
//			idx++;
//		}
//	}
//
//	glBindBuffer(GL_ARRAY_BUFFER, vbo_traces);
//	glBufferData(GL_ARRAY_BUFFER, size * sizeof(GLfloat), data, GL_STATIC_DRAW);
//
//	glBindBuffer( GL_ARRAY_BUFFER, 0 );
//
//	delete[] data;
//
};

void ParticleImplementationCPU::updateParticlesVBO()
{
	glBindBuffer(GL_ARRAY_BUFFER, vbo_particles);

	GLfloat* data;
	data = (GLfloat*) glMapBuffer( GL_ARRAY_BUFFER, GL_WRITE_ONLY);

	if (data != (GLfloat*) NULL)
	{
		int idx = 0;
		for( p_i = m_particleList.begin(); p_i != m_particleList.end(); p_i++)
		{
			Point p = p_i->Position;
			for( int i=0; i<3; i++) {
				data[3*idx+i] = p.c[i];
			}
			idx++;
		}
		glUnmapBuffer(GL_ARRAY_BUFFER);
	} else {
		//printf("MapBuffer error!\n");
	}

	glBindBuffer(GL_ARRAY_BUFFER, 0);
};


void ParticleImplementationCPU::deleteParticleVBO()
{
	glBindBuffer(1, vbo_particles);
	glDeleteBuffers(1, &vbo_particles);
	vbo_particles=0;

	glBindBuffer(1, vbo_traces);
	glDeleteBuffers(1, &vbo_traces);
	vbo_traces=0;

	m_traces_offset.clear();
};


