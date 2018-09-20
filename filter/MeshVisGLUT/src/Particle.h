
#pragma once

#ifndef PARTICLE_H
#define PARTICLE_H

// includes, GL
#include <GL/glew.h>
#include <GL/glut.h>

#include "Typedefs.h"
#include "Vector.h"

#include "common.cuh"

class Particle {

public:

	Particle() : Cell(0), Stepsize(0.0)
	{
	};
	Particle( const Point& pos ) : Position(pos), Cell(0), Stepsize(0.0)
	{ 
		//Trace.push_back(pos);
	};
	Particle( const float4& pos ) : Position(pos), Cell(0), Stepsize(0.0)
	{
		//Trace.push_back(pos);
	};
	~Particle() 
	{ 
		//Trace.clear();
	};

	Point Position;
	int* Cell; // store Cell indices for the timesteps
	//std::list<Point> Trace;
	bool bOutOfField;
	float Stepsize;

	void move( Vector direction )
	{
		Position.x += direction.x;
		Position.y += direction.y;
		Position.z += direction.z;
		//Trace.push_back(Position);
	};
};

typedef std::list<Particle>				ParticleList;
typedef std::list<Particle>::iterator	ParticleListIterator;
typedef std::list<Point>				PointList;
typedef	std::list<Point>::iterator		PointListIterator;

#endif
