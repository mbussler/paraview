#include "Tetrahedron.h"

Tetrahedron::Tetrahedron(int index)
{
	x1=y1=z1=0;
	m_n1=m_n2=m_n3=m_n4=0;
	m_index=index;
	det=-1;
};

Tetrahedron::~Tetrahedron() 
{
	x1=y1=z1=0;
	m_n1=m_n2=m_n3=m_n4=0;
	m_index=0;
	det=-1;
};

void Tetrahedron::setPhysicalCoordinates(const Point& p1, const Point& p2, const Point& p3, const Point& p4)
{
	double x2,y2,z2, x3,y3,z3, x4,y4,z4;
	
	x1=p1.x; y1=p1.y; z1=p1.z;
	x2=p2.x; y2=p2.y; z2=p2.z;
	x3=p3.x; y3=p3.y; z3=p3.z;
	x4=p4.x; y4=p4.y; z4=p4.z;

	det = \
		(x2-x1)*( (y3-y1)*(z4-z1)-(z3-z1)*(y4-y1) ) +
		(x3-x1)*( (y1-y2)*(z4-z1)-(z1-z2)*(y4-y1) ) +
		(x4-x1)*( (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1) );

	a11= (z4-z1)*(y3-y4) - (z3-z4)*(y4-y1);
	a21= (z4-z1)*(y1-y2) - (z1-z2)*(y4-y1);
	a31= (z2-z3)*(y1-y2) - (z1-z2)*(y2-y3);
	
	a12= (x4-x1)*(z3-z4) - (x3-x4)*(z4-z1);
	a22= (x4-x1)*(z1-z2) - (x1-x2)*(z4-z1);
	a32= (x2-x3)*(z1-z2) - (x1-x2)*(z2-z3);
	
	a13= (y4-y1)*(x3-x4) - (y3-y4)*(x4-x1);
	a23= (y4-y1)*(x1-x2) - (y1-y2)*(x4-x1);
	a33= (y2-y3)*(x1-x2) - (y1-y2)*(x2-x3);
};

void Tetrahedron::calculateNaturalCoordinates( const Point& p, double& c1, double& c2, double& c3, int& next)
{

	c1 = (a11*(p.x-x1) + a12*(p.y-y1) + a13*(p.z-z1)) / det;  // 1->2
	c2 = (a21*(p.x-x1) + a22*(p.y-y1) + a23*(p.z-z1)) / det;  // 1->3
	c3 = (a31*(p.x-x1) + a32*(p.y-y1) + a33*(p.z-z1)) / det;  // 1->4
	double c4 = 1-c1-c2-c3;

	//calculate worst violator
	if( (c1<0)||(c2<0)||(c3<0)||(c4<0) )
	{
		if ( c1<c2 && c1<c3 && c1<c4 ) {
			next=m_n2; // next cell shares points 0,2,3
		}
		else if( c2<c3 && c2<c4 ) {
			next=m_n3; // next cell shares points 0,1,3
		}
		else if( c3<c4 ) {
			next=m_n4; // next cell shares points 0,1,2
		}
		else {
			next=m_n1; // next cell shares points 1,2,3
		}
	}
	else
	{
		// point lies within tetrahedron -> terminate
		next=m_index;
	}

};


