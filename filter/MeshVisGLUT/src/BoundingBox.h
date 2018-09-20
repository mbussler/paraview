/**
 *
 *	\brief BoundingBox class
 *  
 *  
 *
 *	\author Michael Bu�ler
 */

#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include "Typedefs.h"
#include "InlineFunctions.h"

#include <vector>
using namespace std;

class BoundingBox{

public:

	/// constructor
	BoundingBox();
	BoundingBox( Point min, Point max);

	/// destructor
	~BoundingBox();

	void init();

	void updateBoundingBox(const Point& point);
	REAL getLengthOfDim(char dim);

	void getVertices( vector<Point>& p );

	Point RandomPointInBox();

public:

	Point point_max;
	Point point_min;

private:
};

#endif
