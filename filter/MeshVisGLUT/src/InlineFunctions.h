#pragma once
#ifndef INLINEFUNCTIONS_H
#define INLINEFUNCTIONS_H

#include "Vector.h"
#include "Typedefs.h"
#include <list>
#include <vector>
#include <helper_math.h>

using namespace std;

//std::vector<bool> bytesToBitset( unsigned char *data, int numBytes )
//{
//    std::vector<bool> b;
//	b.resize( numBytes * CHAR_BIT );
//
//    for(int i = 0; i < numBytes; ++i)
//    {
//        unsigned char cur = data[i];
//        int offset = i * CHAR_BIT;
//
//        for(int bit = 0; bit < CHAR_BIT; ++bit)
//        {
//            b[offset] = cur & 1;
//            ++offset;   // Move to next bit in b
//            cur >>= 1;  // Move to next bit in array
//        }
//    }
//    return b;
//}


#define TOGGLE(x) x=!x
#define checkFlag(f,x) ((f&x)==x)
#define toggleFlag(f,x) f^=x

#define drawPoint(p) glVertex3f(p.x, p.y, p.z)
#define drawColor(c) glColor3f(c.r, c.g, c.b)

#define PI 3.1415926536

inline double Random() 
{ 
	return (rand() / (float) RAND_MAX); 
}

inline double Random(double lower_bound, double upper_bound) 
{
	return Random() * (upper_bound-lower_bound)+lower_bound;
}
inline int signof(double a) 
{ 
	return (a>0) ? 1 : (a<0 ? -1 : 0);
}

/// calculates the distance between two points
inline float distanceToPoint( Point p, Point r)
{
	return VectorLength(Vector(p.x -r.x ,p.y -r.y ,p.z -r.z ));
}

// Colas Schretter <cschrett@ulb.ac.be> 2006
inline Point RandomPointOnSphere(Point p, double r)
{
	float x,y,z,s;
	do {
		x = Random(-1,1);
		y = Random(-1,1);
		s = x*x+y*y;
	} while (s > 1);

	const float t = 2*sqrt(1-s);

	x *= t;
	y *= t;
	z = 2*s - 1;

	return Point(r*x+p.x,r*y+p.y,r*z+p.z);
};

inline Point RandomPointInSphere(Point p, double r)
{
	return RandomPointOnSphere( p, Random(0.0f, r));
};

inline Point RandomPointOnRing( Point p, double r)
{
	float phi=Random(0, 2*PI);
	return Point( r*cosf(phi)+p.x, r*sinf(phi)+p.y, p.z);
};

inline Color3f RandomColor()
{
	return Color3f( Random(), Random(), Random());
};

inline Vector Pt2Vec( Point p) { return Vector( p.x, p.y, p.z); };
inline Point Vec2Pt( const Vector v) { return Point(v.x, v.y, v.z); };


//template< typename T >
//inline void safeIncrement( list<T>& list, list<T>::iterator& iter )
//{
//	iter++;
//	if( iter == list.end() )
//		iter = list.begin();
//}

inline void safeIncrement( std::list<TimePoint>& list, std::list<TimePoint>::iterator& iter )
{
	iter++;
	if( iter == list.end() )
		iter = list.begin();
};


//struct first {
//	typedef T * pointer;
//};

//using namespace std;
//template<typedef T>
//
//static void safeIncrement( const list<typedef T>& ls, list<T>::iterator& iter )
//{
//	if( ++iter == ls.end() )
//		iter--;
//};

// create a color ramp
inline void colorRamp(float t, float *r)
{
    const int ncolors = 7;
    float c[ncolors][3] = {
        { 1.0, 0.0, 0.0, },
        { 1.0, 0.5, 0.0, },
	    { 1.0, 1.0, 0.0, },
	    { 0.0, 1.0, 0.0, },
	    { 0.0, 1.0, 1.0, },
	    { 0.0, 0.0, 1.0, },
	    { 1.0, 0.0, 1.0, },
    };
    t = t * (ncolors-1);
    int i = (int) t;
    float u = t - floor(t);
    r[0] = lerp(c[i][0], c[i+1][0], u);
    r[1] = lerp(c[i][1], c[i+1][1], u);
    r[2] = lerp(c[i][2], c[i+1][2], u);
}

#endif
