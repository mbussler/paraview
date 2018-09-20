
#include "ParticleImplementationCuda.h"
#include <iostream>

ParticleImplementationCuda::ParticleImplementationCuda()
{
	_init();
};

ParticleImplementationCuda::~ParticleImplementationCuda()
{
	_finalize();
};

void ParticleImplementationCuda::_init()
{
	m_numParticles = 0;
	m_timestep = 0.02f;
	m_seeder = new OnSphereParticleSeeder();
	m_seeder->center = Point(0.0f, 0.7f, 0.00f);

	m_particleCounter = 0;
	particleArrayPos = 0;

#ifdef NOGL
	allocateArray((void**) &dParticles, MAX_PARTICLES * sizeof( float4 ));
#else
	// allocate GPU data
	unsigned int memSize = MAX_PARTICLES * sizeof( float4 );
	m_posVbo = createVBO(memSize);
	m_colorVBO = createVBO( memSize );

	registerGLBufferObject(m_posVbo, &m_cuda_posvbo_resource);
#endif
};

void ParticleImplementationCuda::_finalize()
{
#ifdef NOGL
	freeArray( dParticles );
#else
	unregisterGLBufferObject(m_cuda_posvbo_resource);
	glDeleteBuffers(1, (const GLuint*)&m_posVbo);
	glDeleteBuffers(1, (const GLuint*)&m_colorVBO);
#endif
};


// create some particles
void ParticleImplementationCuda::createParticles()
{
	if(m_seeder->type == PointSeederType)
	{
		createSingleParticle();
		return;
	}
	else
	{
		particleArrayPos = 0;
	}

	// set count
	int count=32;

	switch (m_seeder->type)
	{
	case PointSeederType:
		count = 32;
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

	// calculate blocks for new particles
	int numBlocks	= ( count + ParticlesPerBlock-1 ) / ParticlesPerBlock;
	int totalBlocks = ( m_numParticles + ParticlesPerBlock-1 ) / ParticlesPerBlock;

	// create new particles on host-side in page-locked memory for faster transfer
	float4* newParticles;
	size_t sp = count * sizeof(float4);
	allocatePageLockedArray((void**)&newParticles, sp, false );
	m_seeder->createPoints( newParticles, count);

	// store block-offsets also in page-locked host memory
	int* hOffsets;
	size_t sb = (numBlocks+1)*sizeof(int);
	allocatePageLockedArray((void**)&hOffsets, sb , false);
	calculateFreeBlocks( hOffsets, numBlocks, totalBlocks );

	// distribute blocks of particles to GPU-Memory
	int* dOffsets;
	allocateArray( (void**)&dOffsets, sb);
	copyArrayToDeviceAsync( dOffsets, hOffsets, sb, m_tsc->getExecutionStream());

	// copy new positions to VBO
#ifdef NOGL
	for( int i=0; i<numBlocks; i++ )
	{
		copyArrayToDeviceAsync( dParticles + hOffsets[i+1], &newParticles[32*i], 32*sizeof(float4), m_tsc->getExecutionStream());
	}

#else
	updatePositionsVBO( newParticles, count, hOffsets );
	updateColorVBO( count, hOffsets );
	dParticles = (float4*) mapGLBufferObject(&m_cuda_posvbo_resource);
#endif

	// update cell indices in tsc
	m_tsc->locateParticlesCuda( dParticles, count, dOffsets );

#ifndef NOGL
	unmapGLBufferObject(m_cuda_posvbo_resource);
#endif

	// increment particle counter
	m_numParticles = totalBlocks * ParticlesPerBlock;

	// clean up
	freeArray( dOffsets );
	freePageLockedHostMemory( hOffsets);
	freePageLockedHostMemory(newParticles);

};

// create some particles
void ParticleImplementationCuda::createSingleParticle()
{
	// create new particle
	float4* newParticle;
	size_t sp = sizeof(float4);
	allocatePageLockedArray((void**)&newParticle, sp, false );
	m_seeder->createPoints(newParticle, 1);

	// copy new positions to VBO
#ifdef NOGL
	copyArrayToDeviceAsync( dParticles + particleArrayPos, newParticle, sizeof(float4), m_tsc->getExecutionStream());
#else
	unregisterGLBufferObject(m_cuda_posvbo_resource);
	glBindBuffer(GL_ARRAY_BUFFER, m_posVbo);
	glBufferSubData(GL_ARRAY_BUFFER, particleArrayPos*sizeof(float4), sizeof(float4), newParticle );
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	registerGLBufferObject(m_posVbo, &m_cuda_posvbo_resource);

	// fill color buffer
	float4 color;
	float t = ( m_particleCounter++ ) / (float) 128000;
	m_particleCounter %= 128000;
	colorRamp(t, &color.x);
	color.w = 1.0f;

	glBindBuffer(GL_ARRAY_BUFFER, m_colorVBO);
	glBufferSubData(GL_ARRAY_BUFFER, particleArrayPos*sizeof(float4), sizeof(float4), &color );
	glBindBuffer(GL_ARRAY_BUFFER, 0);	
	dParticles = (float4*) mapGLBufferObject(&m_cuda_posvbo_resource);
#endif

	particleArrayPos++;
	m_numParticles = particleArrayPos;

	// update cell indices in tsc
	m_tsc->locateParticlesCuda( dParticles, m_numParticles, 0 );

#ifndef NOGL
	unmapGLBufferObject(m_cuda_posvbo_resource);
#endif

	freePageLockedHostMemory(newParticle);
};


// create some particles
void ParticleImplementationCuda::createParticlesTest( int count )
{
	// set count
	//count=1000000;

	// create new particles
	float4* newParticles;
	size_t sp = count * sizeof( float4);
	allocatePageLockedArray((void**)&newParticles, sp, false);
	m_seeder->createPoints( newParticles, count);

	// copy positions
#ifdef NOGL
	copyArrayToDeviceAsync( dParticles, newParticles, sp, m_tsc->getExecutionStream());
#else
	updatePositionsVBO( newParticles, count );
	updateColorVBO( count );
	dParticles = (float4*) mapGLBufferObject(&m_cuda_posvbo_resource);
#endif

	// update cell indices in tsc
	m_tsc->locateParticlesCuda( dParticles, count, 0 );

#ifndef NOGL
	unmapGLBufferObject(m_cuda_posvbo_resource);
#endif

	// increment particle counter
	m_numParticles = count;

	// clean up
	freePageLockedHostMemory(newParticles);
};


/* Trace Particles through the flow field */
void ParticleImplementationCuda::ParticleStep()
{
#ifndef NOGL
	dParticles = (float4*) mapGLBufferObject(&m_cuda_posvbo_resource);
#endif

	m_tsc->advectParticles( dParticles, m_timestep, m_numParticles );

#ifndef NOGL
	unmapGLBufferObject(m_cuda_posvbo_resource);
#endif
};


void ParticleImplementationCuda::clearParticles()
{
	m_numParticles = 0;
	particleArrayPos = 0;
};

void ParticleImplementationCuda::switchSeeder()
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
		//	case SphereSeederType:
		//		m_seeder = new SphereParticleSeeder();
		//		break;
	}

	m_seeder->center = old_seeder->center;
	m_seeder->size = old_seeder->size;


	// HOLODEMO: clear particles, if seeder switched
	clearParticles();

	delete old_seeder;
}

uint
ParticleImplementationCuda::createVBO(uint size)
{
	GLuint vbo;
	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, size, 0, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	return vbo;
}

void
ParticleImplementationCuda::dumpParticles(uint start, uint count)
{
	float* hPos = (float*)malloc( sizeof(float)*4*count);
	copyArrayFromDevice(hPos, 0, sizeof(float)*4*count, &m_cuda_posvbo_resource);

	for(uint i=start; i<start+count; i++) {
		printf("pos: (%.4f, %.4f, %.4f, %.4f)\n", hPos[i*4+0], hPos[i*4+1], hPos[i*4+2], hPos[i*4+3]);
	}
	free( hPos );
}


void ParticleImplementationCuda::calculateFreeBlocks( int* hOffsets, int numBlocks, int& totalBlocks )
{
	hOffsets[0] = numBlocks;

	// get free blocks and distribute new positions
	boost::dynamic_bitset<unsigned char> blocks;
	m_tsc->getOccupiedBlocks( totalBlocks, blocks );

	m_activeParticles = ( blocks.size() - blocks.count()) * ParticlesPerBlock;

	// store offsets to blocks array on host
	int off = blocks.find_first();

	for( int i=0; i<numBlocks; i++)
	{
		if( off < 0 ) { off = totalBlocks++; }

		hOffsets[i+1] = off * ParticlesPerBlock;
		off = blocks.find_next( off );
	}
}

// copy new particle positions to GPU
void ParticleImplementationCuda::updatePositionsVBO( float4* newParticles, int count, int* hOffsets)
{
	unregisterGLBufferObject(m_cuda_posvbo_resource);
	glBindBuffer(GL_ARRAY_BUFFER, m_posVbo);

	if( hOffsets )
	{
		for( int i=0; i<hOffsets[0]; i++ )
		{
			glBufferSubData(GL_ARRAY_BUFFER, hOffsets[i+1]*sizeof(float4), 32*sizeof(float4), &newParticles[32*i] );
		}
	}
	else
	{
		glBufferData(GL_ARRAY_BUFFER, count*sizeof(float4), newParticles, GL_DYNAMIC_DRAW);
	}

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	registerGLBufferObject(m_posVbo, &m_cuda_posvbo_resource);
}

void ParticleImplementationCuda::updateColorVBO( int count, int* hOffsets )
{
	// fill color buffer
	float *colors = new float[4*count];
	float *ptr = colors;

	for( int i=0; i< count; i++)
	{
		//Color3f col = RandomColor();
		float t = ( m_particleCounter++ ) / (float) 128000;
		m_particleCounter %= 128000;

		colorRamp(t, ptr);
		ptr+=3;
		*ptr++ = 1.0f;
	}

	glBindBuffer(GL_ARRAY_BUFFER, m_colorVBO);

	if( hOffsets)
	{
		for( int i=0; i<hOffsets[0]; i++ )
		{
			glBufferSubData(GL_ARRAY_BUFFER, hOffsets[i+1]*sizeof(float4), 32*sizeof(float4), &colors[4*32*i] );
		}
	}
	else
	{
		glBufferData(GL_ARRAY_BUFFER, count*sizeof(float4), colors, GL_DYNAMIC_DRAW );
	}

	glBindBuffer(GL_ARRAY_BUFFER, 0);
}
