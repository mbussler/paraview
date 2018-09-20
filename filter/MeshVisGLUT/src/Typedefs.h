
#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <vector_types.h>
#include <vector_functions.h>

typedef unsigned long DWORD;

#ifdef TETLIBRARY 
	#include "tetgen.h"
#else             
	#define REAL float
#endif

class Color3f {
public:
	union{
		struct {
			float r;
			float g;
			float b;
			float a;
		};
		float c[4];
	};

	Color3f() {};
	Color3f( float _r, float _g, float _b)
		{ r=_r; g=_g; b=_b;a=1.0;}
	Color3f( float _r, float _g, float _b, float _a)
		{ r=_r; g=_g; b=_b; a=_a;}
	~Color3f() {};
};

static Color3f C3fGreen = Color3f(0.5f, 1.0f, 0.5f);

const Color3f trace_colors[] = {
	Color3f( 0, 0, 0 ),
	Color3f( 1, 0, 0 ),
	Color3f( 0, 0.7, 0 ),
	Color3f( 0, 0, 1 ),
	Color3f( 0.5, 0.5, 0 ),
};
const int trace_colors_count = 5;

class Point 
{
public:
	union{
		struct {
			REAL x;
			REAL y;
			REAL z;
		};
		REAL c[3];
	};

	Point() : 
        x(0.0f), 
        y(0.0f), 
        z(0.0f) 
    {};
    Point( REAL _x, REAL _y, REAL _z) :
        x(_x),
        y(_y),
        z(_z)
	{};
    Point( const float4& v ) :
        x(v.x),
        y(v.y),
        z(v.z)
        {};
    ~Point() {};

	bool operator !=(Point const& other)
	{
		return (c[0] != other.c[0] || c[1] != other.c[1] || c[2] != other.c[2]);
	}
    float4 toFloat4()
    {
        return make_float4(x,y,z,1.0f);
    };
};

struct Node {
	int index;
	Point p;
};

/**
 * A TreeNode of the Kd-Tree
 */
struct TreeNode {
	/// Defines the split-plane for this treenode
	char splitDim;
	/// Defines the split-value for this treenode
	REAL splitValue;

	int index;
	/// The left successor
	TreeNode* leftNode;
	/// The right successor
	TreeNode* rightNode;

	int a,b; // lower and upper bound for the indices in this Treenode
};

class TimePoint
{
public:
	int index;
	double time;
};

static bool TimePointSortPredicate(const TimePoint& tp1, const TimePoint& tp2)
{
	return tp1.time < tp2.time;
};

#endif
