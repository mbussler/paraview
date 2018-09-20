
#pragma once

#ifndef TETMESHGL_H
#define TETMESHGL_H

// includes, GL
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <glh_linear.h>

#define BUFFER_OFFSET(bytes) ((GLubyte*) NULL + (bytes))

#include "TetMesh.h"
#include <list>
#include <boost/dynamic_bitset.hpp>

#include <boost/shared_ptr.hpp>

/**
 *	\brief TetMeshGL class
 *  
 *	Extension of TetMesh with additional OpenGL Visualization functions
 *
 *	\author Michael Bussler
*/

// Colors of the traversed cells
const int		tc_cellcolor_count = 4;
const Color3f	tc_cellcolor[] = {
	Color3f(.3, .3, 1 ), 
	Color3f(.4, .4, 1 ), 
	Color3f(.6, .6, 1 ), 
	Color3f(.8, .8, 1 ) 
};

const Color3f tc_startcell_color = Color3f(1.0, 0.3, 0.2);

// defines the order by which the triangles of a tetrahedron cell are drawn
const int triangles[][3] = { {0,2,1}, {1,2,3}, {2,0,3}, {3,0,1}};

enum {
	DRAW_BOX 				= 1,
	DRAW_VERTICES			= 2,
	DRAW_VELOCITIES 		= 4,
	DRAW_CELLS				= 8,
	DRAW_TRAVERSED_CELLS 	= 16,
	DRAW_OUTER_FACES		= 32
};

inline bool compare_depth ( std::pair<int,float> a, std::pair<int,float> b)
{
	return a.second < b.second;
};

/// Mesh type definition
class TetMeshGL : public TetMesh
{
public:

	TetMeshGL();
	~TetMeshGL();

	void init();
	void initVBOs();
	void deleteVBOs();

	void draw( DWORD flags);
	void createTraversedCellsVBO( const bool &use_modelview = false, const glh::matrix4 &m = glh::matrix4());


	boost::dynamic_bitset<>* getTraversedCellsList()
	{
		return &tc_indices;
	};

	void resetTraversedCells()
	{
		tc_indices.reset();
	};


private:

	void createPointsVBO();
	void createBoxVBO();
	void createCellsIBO();
	void createOuterFacesIBO();

	float vec_scale;

	// vbo variables
	GLuint	vbo_p,    // Holds point coordinates
			vbo_box,  // Box indices
			ibo_c,    // Hold cell indices
			ibo_tc,   // Hold traversed cell indices
			ibo_oc,   // Hold the indices of the outer faces
			vbo_tc;	  // VBO for traversed cells


	boost::dynamic_bitset<> tc_indices;

	int traversed_cells_count;
	int outer_faces_count;
	int draw_cells_count;
	int draw_nodes_count;
	int draw_vel_count;
};
typedef boost::shared_ptr<TetMeshGL> TetMeshGLPtr;


#endif
