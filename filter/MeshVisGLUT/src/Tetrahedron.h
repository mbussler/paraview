
#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

#include "Typedefs.h"
#include "Vector.h"

class Tetrahedron
{
public:
	Tetrahedron(int index);
	~Tetrahedron();

	void setPhysicalCoordinates(const Point& p1, const Point& p2, const Point& p3, const Point& p4);
	void calculateNaturalCoordinates( const Point& p, double& c1, double& c2, double& c3, int& next);

	void setNeighbors(int n1, int n2, int n3, int n4)	
		{ m_n1=n1;m_n2=n2;m_n3=n3;m_n4=n4;};

private:

	// reference point
	double x1,y1,z1;

	// interpolation matrix
	double a11,a12,a13, a21,a22,a23, a31,a32,a33;

	// Determinant
	double det;

	// Neighbor cell indices
	int m_n1,m_n2,m_n3,m_n4;

	// Tetrahedron cell index
	int m_index;

};



#endif
