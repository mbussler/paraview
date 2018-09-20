
#pragma once

#ifndef PARTICLE_SEEDER_H
#define PARTICLE_SEEDER_H

// includes, GL
#include <GL/glew.h>
#include <GL/freeglut.h>

#include "Typedefs.h"
#include "InlineFunctions.h"
#include "BoundingBox.h"

#include "common.cuh"

typedef enum SeederType
{
	PointSeederType = 0,
	SphereSeederType,
	RingSphereSeederType,
	OnSphereSeederType,
	CubeSeederType
};

class ParticleSeeder
{
public:
	ParticleSeeder()
	{ 
		center = Point(0,0.75,0);
		type = PointSeederType;
		size = 0.1;
	};
	~ParticleSeeder()
	{};

	Point center;
	SeederType type;
	float size;

	virtual void createPoint(Point& p)
	{ p = center; };

	virtual void draw()
	{
		glColor3f(1.0, 0.6f, 0.8f);
		glPointSize(2.0f);
		glBegin(GL_POINTS);
			drawPoint( center );
		glEnd();
	};

	virtual void createPoints( float4* ps, int number )
	{
		assert( number > 0 );
		assert( ps );

		Point p;
		
		for( int i=0; i<number; i++)
		{
			createPoint(p);
			ps[i] = p.toFloat4();
		}
	};

};

class SphereParticleSeeder : public ParticleSeeder
{
public:
	SphereParticleSeeder()
	{ 
		type = SphereSeederType;
	};
	~SphereParticleSeeder() 
	{};

	virtual void createPoint(Point& p)
	{ p = RandomPointInSphere( center, size/2.0 ); };
	virtual void draw()
	{
		glColor3f(1.0, 0.6f, 0.8f);
		glPushMatrix();
		glTranslatef( center.x, center.y, center.z);
		glutWireSphere( size/2.0, 10, 10);
		glPopMatrix();
	};

};

class OnSphereParticleSeeder : public SphereParticleSeeder
{
public:
	OnSphereParticleSeeder()
	{
		type = OnSphereSeederType;
	};
	~OnSphereParticleSeeder()
	{};

	virtual void createPoint(Point& p)
	{ p = RandomPointOnSphere( center, size/2.0 ); };

	virtual void createPoints( float4* ps, int number )
	{
		assert( number > 0 );
		assert( ps );

		float rad = size/2.0f;
		Point p = center;

		float dl = M_PI*(3.0f-sqrt(5.0f));
		float l = 0;
		float dz = 2.0f/ (float)number;
		float z = 1 - dz/2.0f;

		for( int i=0; i<number; i++)
		{
			float r = sqrt(1-z*z);
			ps[i] = make_float4( rad*r*cos(l) + p.x, rad*r*sin(l) + p.y, rad*z + p.z, 1.0f);
		    z = z - dz;
		    l = l + dl;
		}
	};
};

class RingParticleSeeder : public SphereParticleSeeder
{
public:
	RingParticleSeeder()
	{ 
		type = RingSphereSeederType;
	};
	~RingParticleSeeder() 
	{};

	virtual void createPoint(Point& p)
	{ p = RandomPointOnRing( center, size/2.0 ); };

	void draw()
	{
		glColor3f(1.0, 0.6f, 0.8f);
		glPushMatrix();
		glTranslatef( center.x, center.y, center.z);
		glutWireTorus(0.01, size/2.0, 10, 15);
		glPopMatrix();
	};

	virtual void createPoints( float4* ps, int number )
	{
		assert( number > 0 );
		assert( ps );

		float dPhi = 2.0f*M_PI / (float) number;
		float phi = 0.0f;
		float r = size/2.0f;
		Point p = center;

		for( int i=0; i<number; i++)
		{
			ps[i] = make_float4( r*cosf(phi) + p.x , r*sinf(phi)+ p.y, p.z, 1.0f);
			phi += dPhi;
		}
	};
};

class CubeParticleSeeder : public ParticleSeeder
{
public:
	CubeParticleSeeder()
	{ 
		type = CubeSeederType;
	};
	~CubeParticleSeeder() {};

	virtual void createPoint(Point& p)
	{
		float lengthDiv2 = size/2.0f;
		float x = Random( center.x - lengthDiv2, center.x + lengthDiv2);
		float y = Random( center.y - lengthDiv2, center.y + lengthDiv2);
		float z = Random( center.z - lengthDiv2, center.z + lengthDiv2);
		p = Point( x,y,z );
	};
	void draw()
	{
		glColor3f(1.0, 0.6f, 0.8f);
		glPushMatrix();
		glTranslatef( center.x, center.y, center.z);
		glutWireCube( size);
		glPopMatrix();
	};

	virtual void createPoints( float4* ps, int number )
	{
		assert( number > 0 );
		assert( ps );

		// calculate the nearest value of the third root
		int row = (int) floor(pow(number, 1.0f/3.0f));
		int pos = 0;

		float dl = size / (float)(row-1);
		Point p = Point( center.x - size/2.0f, center.y - size/2.0f, center.z - size/2.0f);

		for( int i=0; i<row; i++)
			for( int j=0; j<row; j++)
				for( int k=0; k<row; k++)
					ps[pos++] = make_float4( p.x + i*dl, p.y + j*dl, p.z + k*dl, 1.0f );

		for( int i=0; i<(number-row*row*row); i++)
		{
			ps[pos++] = make_float4( center.x, center.y, center.z, 1.0f );
		}
	}
};

#endif
