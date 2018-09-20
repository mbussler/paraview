#include "BoundingBox.h"

#include <vector>
using namespace std;

/// constructor
BoundingBox::BoundingBox()
{
	init();
};

BoundingBox::BoundingBox(Point min, Point max )
{
	// initialize
	point_max = max;
	point_min = min;
};

void BoundingBox::init()
{
	point_max = Point(-1e10,-1e10,-1e10);
	point_min = Point(1e10,1e10,1e10);
};

/// destructor
BoundingBox::~BoundingBox(){};

void BoundingBox::updateBoundingBox(const Point& point)
{
	if( point.x > point_max.x ) { point_max.x = point.x; }
	if( point.y > point_max.y ) { point_max.y = point.y; }
	if( point.z > point_max.z ) { point_max.z = point.z; }
	
	if( point.x < point_min.x ) { point_min.x = point.x; }
	if( point.y < point_min.y ) { point_min.y = point.y; }
	if( point.z < point_min.z ) { point_min.z = point.z; }
};

REAL BoundingBox::getLengthOfDim(char dim)
{
	return (point_max.c[dim] - point_min.c[dim]);
};

void BoundingBox::getVertices( vector<Point>& p )
{
	p.resize(8);
	p[0] = Point( point_min.x, point_min.y, point_min.z);
	p[1] = Point( point_max.x, point_min.y, point_min.z);
	p[2] = Point( point_max.x, point_max.y, point_min.z);
	p[3] = Point( point_min.x, point_max.y, point_min.z);
	p[4] = Point( point_min.x, point_min.y, point_max.z);
	p[5] = Point( point_max.x, point_min.y, point_max.z);
	p[6] = Point( point_max.x, point_max.y, point_max.z);
	p[7] = Point( point_min.x, point_max.y, point_max.z);
};

Point BoundingBox::RandomPointInBox()
{
	return Point(   Random( point_min.x, point_max.x ),
					Random( point_min.y, point_max.y ),
					Random( point_min.z, point_max.z ) );
};
