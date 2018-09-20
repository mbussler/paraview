//No Vector.cpp necessary
#ifndef VECTOR_H
#define VECTOR_H

#include "math.h"
#include "stdlib.h"	//needed for random vector
#include "Typedefs.h"

class Vector
{

public:
	union{
		struct{
			float x,y,z;
		};
		float c[3];
	};

	Vector()
	{
		x=0;
		y=0;
		z=0;
	};
	Vector(const float a, const float b, const float c):x(a),y(b),z(c)
	{
	};
	Vector(const Vector& v):x(v.x), y(v.y), z(v.z)
	{
	};
	Vector(Point& p):x(p.x), y(p.y), z(p.z)
	{
	};
	~Vector()
	{
	};

	const Point pt() const
	{
		return Point(x,y,z);
	};

	Vector& operator =(const Vector& v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
		return(*this);
	};

	Vector& operator +=(const Vector& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
		return(*this);

	};
	Vector& operator -=(const Vector& v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return(*this);

	};
	Vector& operator /=(const float f)
	{
		x /= f;
		y /= f;
		z /= f;
		return(*this);
	};
	Vector& operator *=(const float f)
	{
		x *= f;
		y *= f;
		z *= f;
		return(*this);
	};
};

	inline bool operator ==(Vector& a, Vector& b)
	{
		return ((a.x==b.x) && (a.y == b.y) && (a.z == b.z));
	};

	inline bool operator !=(Vector& a, Vector& b)
	{
		return ((a.x!=b.x) || (a.y != b.y) || (a.z != b.z));
	};

	inline Vector operator + (const Vector& a, const Vector& b)
	{
		return Vector(a.x+b.x, a.y+b.y, a.z+b.z);
	};

	inline Vector operator - (const Vector& a, const Vector& b)
	{
		return Vector(a.x-b.x, a.y-b.y, a.z-b.z);
	};

	inline Vector operator * (const Vector& a, const float f)
	{
		return Vector(a.x*f, a.y*f, a.z*f);
	};

	inline Vector operator * (const float f, const Vector& a)
	{
		return Vector(a.x*f, a.y*f, a.z*f);
	};

		inline Vector operator / (const Vector& a, float f)
	{
		return Vector(a.x/f, a.y/f, a.z/f);
	};

	inline float VectorLength(const Vector& v)
	{
		return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
	};

	inline float VectorLengthSquare(const Vector& v)
	{
		return (v.x*v.x + v.y*v.y + v.z*v.z);
	};
	inline Vector VectorNormalize(const Vector& v)
	{	
		return v / sqrtf(v.x*v.x + v.y*v.y + v.z*v.z);
	};
	inline Vector VectorNormalizeSafe(const Vector& v)
	{
		return v / (sqrtf(v.x*v.x + v.y*v.y + v.z*v.z)+0.0001f);
	};
	inline Vector VectorCross(const Vector& a, const Vector& b)
	{
		return Vector(a.y * b.z - a.z * b.y,
					  a.z * b.x - a.x * b.z,
					  a.x * b.y - a.y * b.x);
	};

	inline float  VectorDot(const Vector& a, const Vector& b)
	{
		return a.x*b.x + a.y*b.y + a.z*b.z;
	};
	inline float  VectorAngle(const Vector& a, const Vector& b)
	{
		return acosf((a.x*b.x + a.y*b.y + a.z*b.z) /
					 sqrtf((a.x*a.x + a.y*a.y + a.z*a.z) *
						   (b.x*b.x + b.y*b.y + b.z*b.z)));
	};

	inline Vector VectorInterpolateCoords(const Vector& a, const Vector& b, float f)
	{
		return a + f * (b - a);
	};

	inline Vector VectorRandom()
	{
		return Vector((float)(rand()),(float)(rand()),(float)(rand()));
	};
#endif
