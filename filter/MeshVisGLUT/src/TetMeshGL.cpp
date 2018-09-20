#include "TetMeshGL.h"

#include <vector>
using namespace std;

TetMeshGL::TetMeshGL()
{
	vbo_p   = 0;
	vbo_box = 0;
	vbo_tc  = 0;

	ibo_c   = 0;
	ibo_tc  = 0;
	ibo_oc  = 0;

	traversed_cells_count = 0;
	outer_faces_count = 0;
	draw_cells_count = 0;
	draw_nodes_count = 0;
	draw_vel_count = 0;
};

TetMeshGL::~TetMeshGL()
{
#ifndef NOGL
	deleteVBOs();
#endif
};

void TetMeshGL::init()
{
	vec_scale = 0.06;
	assert( cells_count > 0 );
	tc_indices.resize( cells_count );
#ifndef NOGL
	initVBOs();
#endif
};

void TetMeshGL::initVBOs()
{	
	/* Create VBO for vertices and velocities*/
	glGenBuffers(1, &vbo_p);

	/* Create VBO for the Bounding Box */
	glGenBuffers(1, &vbo_box);

	/* Create IBO for Cells */
	glGenBuffers(1, &ibo_c);

	// VBOs for traversed cells rendering
	glGenBuffers(1, &vbo_tc );
	glGenBuffers(1, &ibo_tc );

	glGenBuffers(1, &ibo_oc);

	createBoxVBO();
	createOuterFacesIBO();
	createPointsVBO();
	//createCellsIBO(); <- lazy evaluation
};

void TetMeshGL::deleteVBOs()
{
	// return if extensions were not loaded
	if( glIsBuffer == NULL)
		return;

	if( glIsBuffer(vbo_p))
	{
		glBindBuffer(1, vbo_p);
		glDeleteBuffers(1, &vbo_p);
	}
	if( glIsBuffer(vbo_box))
	{
		glBindBuffer(1, vbo_box);
		glDeleteBuffers(1, &vbo_box);
	}
	if( glIsBuffer(ibo_c))
	{
		glBindBuffer(1, ibo_c);
		glDeleteBuffers(1, &ibo_c);
	}
	if( glIsBuffer(vbo_tc))
	{
		glBindBuffer(1, vbo_tc);
		glDeleteBuffers(1, &vbo_tc);
	}
	if( glIsBuffer(ibo_tc))
	{
		glBindBuffer(1, ibo_tc);
		glDeleteBuffers(1, &ibo_tc);
	}
	if( glIsBuffer(ibo_oc))
	{
		glBindBuffer(1, ibo_oc);
		glDeleteBuffers(1, &ibo_oc);
	}
	vbo_p   = 0;
	ibo_c   = 0;
	vbo_box = 0;
	vbo_tc  = 0;
	ibo_tc  = 0;
	ibo_oc  = 0;
};

void TetMeshGL::draw( DWORD flags)
{
	glPushAttrib(GL_DEPTH_BUFFER_BIT
			| GL_ENABLE_BIT
			| GL_VIEWPORT_BIT
			| GL_POLYGON_BIT ) ;

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glDisable(GL_LIGHTING);

	if( checkFlag(flags, DRAW_BOX) )
	{
		glBindBuffer(GL_ARRAY_BUFFER, vbo_box);
		glVertexPointer(3, GL_FLOAT, 0, BUFFER_OFFSET(0));

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glColor3f(1.0, 1.0, 1.0);
		glDrawArrays( GL_QUADS, 0, 24);
	}

	// use the points-vbo from now on
	glBindBuffer(GL_ARRAY_BUFFER, vbo_p);

	if( checkFlag(flags, DRAW_VELOCITIES))
	{
		glVertexPointer(3, GL_FLOAT, 0, BUFFER_OFFSET(0));
	
		// draw velocity vectors as yellow lines
		glColor3f(1.0, 1.0, 0.0);
		glDrawArrays(GL_LINES, 0, draw_vel_count);
	}

	// set the vertex-pointer to node-coordinates
	glVertexPointer(3, GL_FLOAT, 6*sizeof(GLfloat), BUFFER_OFFSET(0));

	if( checkFlag(flags, DRAW_VERTICES))
	{
		// draw Mesh Nodes
		glColor3f(1.0, 0.0, 0.0);
		glPointSize(2.0f);
		glDrawArrays(GL_POINTS, 0, draw_nodes_count);
	}

	if( checkFlag(flags, DRAW_CELLS))
	{
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_c);
		// draw cells as dark blue lined triangles
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glColor3f(0.663f, 0.776f, 1.0f);
		glDrawElements(GL_TRIANGLES, 3*4*draw_cells_count, GL_UNSIGNED_INT, 0);
	}

	if( checkFlag(flags, DRAW_OUTER_FACES))
	{
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_oc);

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glColor3f(0.6, 0.6, .6);
		glDrawElements(GL_TRIANGLES, 3* outer_faces_count, GL_UNSIGNED_INT, 0);
	}

	if( checkFlag(flags, DRAW_TRAVERSED_CELLS))
	{
		glBindBuffer(GL_ARRAY_BUFFER, vbo_tc);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_tc);

		glVertexPointer(3, GL_FLOAT, 8*sizeof(GLfloat), BUFFER_OFFSET(0));
		glColorPointer( 4, GL_FLOAT, 8*sizeof(GLfloat), BUFFER_OFFSET(4 * sizeof(GLfloat)));

	    glDepthFunc(GL_LEQUAL);

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
//		glDepthRange(0.0, 1.0);
		glColor3f(0.1, 0.1, 1.0);
		glDrawElements(GL_TRIANGLES, 3*4* traversed_cells_count, GL_UNSIGNED_INT, 0);

		glEnableClientState(GL_COLOR_ARRAY);
		glPolygonMode(GL_BACK, GL_FILL);
//		glDepthRange(0.01, 1.0);
		glDrawElements(GL_TRIANGLES, 3*4* traversed_cells_count, GL_UNSIGNED_INT, 0);
		glDisableClientState(GL_COLOR_ARRAY);
	}

	// clean up
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glDisableClientState(GL_VERTEX_ARRAY);

	glPopAttrib();
}

void TetMeshGL::createPointsVBO()
{
	draw_nodes_count = node_count;
	draw_vel_count = 2*node_count;

	/* Buffer size:
		3 floats per vertex;
		3 floats per velocity vector;
		8*3 floats for the surrounding box; */

	unsigned int size = (6 * node_count );

	GLfloat* data = new GLfloat[size];

	for(int i=0; i<node_count; i++)
	{
		// vertex coordinates
		data[6*i+0] = nodes[i].x;
		data[6*i+1] = nodes[i].y;
		data[6*i+2] = nodes[i].z;

		// velocity vector
		data[6*i+3] = nodes[i].x + node_attributes[i].x*vec_scale;
		data[6*i+4] = nodes[i].y + node_attributes[i].y*vec_scale;
		data[6*i+5] = nodes[i].z + node_attributes[i].z*vec_scale;
	}

	glBindBuffer(GL_ARRAY_BUFFER, vbo_p);
	glBufferData(GL_ARRAY_BUFFER, size * sizeof(GLfloat), data, GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	delete data;
};

void TetMeshGL::createBoxVBO()
{
	// box vertices
	vector<Point> pts;	
	box.getVertices( pts );

	const int box_indices[24] = { 0,1,2,3, 4,7,6,5, 0,4,5,1, 3,2,6,7, 0,3,7,4, 1,5,6,2};

	// we have 24 box points with 3 coordinates each
	const int size = 3 * 24;
	GLfloat* data = new GLfloat[ size ];
	
	for( int i=0; i<24; i++) 
	{
		Point p = pts[ box_indices[i] ];
		data[3*i+0] = p.x;
		data[3*i+1] = p.y;
		data[3*i+2] = p.z;
	}

	glBindBuffer(GL_ARRAY_BUFFER, vbo_box);
	glBufferData(GL_ARRAY_BUFFER, size * sizeof(GLfloat), data, GL_STATIC_DRAW);
	
	/* cleanup */
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	delete data;

};

void TetMeshGL::createCellsIBO()
{
	draw_cells_count = cells_count;

	// a cell is drawn by four triangles
	int size = 12 * cells_count;
	
	unsigned int* data = new unsigned int[size];
	for(int i=0; i<cells_count; i++)
	{
		// first triangle
		data[12*i+0] = cells[i].indices[0];
		data[12*i+1] = cells[i].indices[2];
		data[12*i+2] = cells[i].indices[1];
		// second triangle
		data[12*i+3] = cells[i].indices[0];
		data[12*i+4] = cells[i].indices[1];
		data[12*i+5] = cells[i].indices[3];
		// third triangle
		data[12*i+6] = cells[i].indices[1];
		data[12*i+7] = cells[i].indices[2];
		data[12*i+8] = cells[i].indices[3];
		// fourth triangle
		data[12*i+9]  = cells[i].indices[0];
		data[12*i+10] = cells[i].indices[3];
		data[12*i+11] = cells[i].indices[2];
	};

	/* fill VBO for Cell rendering */
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_c);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, size * sizeof(int), data, GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	
	delete data;
};

void TetMeshGL::createOuterFacesIBO()
{
	// a face is drawn as triangle
	int size = 3 * outer_faces_count * sizeof(int);

	/* fill IBO for outer face rendering */
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_oc);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, size, (int*)outer_faces, GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
};

void TetMeshGL::createTraversedCellsVBO( const bool &use_modelview, const glh::matrix4 &m)
{
	glh::matrix4 mv;

	if ( !use_modelview )
	{
		glGetFloatv(GL_MODELVIEW_MATRIX, mv);
	}
	else
	{
		mv = m;
	}

	traversed_cells_count = tc_indices.count();

	// Back to Front Ordering
	std::list<std::pair<int, float> > toSort;

	for( int it = tc_indices.find_first(); it < tc_indices.npos; it = tc_indices.find_next(it))
	{
		std::pair<int,float> item;
		item.first = it;

		// calculate single point for current cell (e.g. center) and transform by MODELVIEW_MATRIX
		Vector vs[4];
		for( int i=0; i<4; i++) {
			vs[i] = Pt2Vec( nodes[ cells[it].indices[i] ] );
		}

		// Calculate Cell Center
		Vector c = vs[0] + 0.5 * (vs[1]-vs[0]) + 0.5 * (vs[2]-vs[0]) + 0.5 * (vs[3]-vs[0]);
	
		// Depth Value
		item.second = c.x*mv[2] + c.y*mv[6] + c.z*mv[10];

		toSort.push_back(item);
	}
	toSort.sort(compare_depth);

	// stores vertex indices, 12 per cell
	int num_cells = 12 * traversed_cells_count;
	unsigned int* cs = new unsigned int[num_cells];

	// store cell vertices and vertice color
	int num_vert = 2*4*4 * traversed_cells_count;
	GLfloat* vs = new GLfloat[num_vert];

	Color3f cellColor;
	float alpha = 0.3;

	int idx=0;

	std::list<std::pair<int,float> >::iterator c_it;
	for( c_it=toSort.begin(); c_it!=toSort.end(); c_it++)
	{
		int i_it = c_it->first;

		// check for starting cell property to set the correct cell color
		// use special color only for cells where the TetWalk ends
//		if (*i_it == tc_startCells.front()) {
//			tc_startCells.pop_front();
//			cellColor = tc_startcell_color;
//		} else {
//			cellColor = tc_cellcolor[*i_it%tc_cellcolor_count];
//		}

		cellColor = tc_cellcolor[i_it%tc_cellcolor_count];
		//cellColor = Color3f(idx/count, 0.0, 0.0);

		/* iterate over the four cell vertices */ 
		for( int i=0; i<4; i++) {

			/* v_idx is the index of the i'th Vertex of cell *c_it */ 
			int v_idx = cells[i_it].indices[i];

			/* store components of the current vertex with index v_idx */
            vs[32*idx +8*i +0] = nodes[v_idx].x;
			vs[32*idx +8*i +1] = nodes[v_idx].y;
			vs[32*idx +8*i +2] = nodes[v_idx].z;

			// store color components for the current vertex
    		vs[32*idx +8*i +4] = cellColor.r;
    		vs[32*idx +8*i +5] = cellColor.g;
    		vs[32*idx +8*i +6] = cellColor.b;

            // set the alpha color component
			vs[32*idx +8*i +7] = alpha;

			/* store indice j for the current triangle i */
			cs[12*idx +3*i +0] = 4*idx + triangles[i][0];
			cs[12*idx +3*i +1] = 4*idx + triangles[i][1];
			cs[12*idx +3*i +2] = 4*idx + triangles[i][2];

		}
		idx++;
	};

	/* fill VBO */
	glBindBuffer(GL_ARRAY_BUFFER, vbo_tc);
	glBufferData(GL_ARRAY_BUFFER, num_vert* sizeof(GLfloat), vs, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	/* fill IBO */
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_tc);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, num_cells* sizeof(int), cs, GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	delete[] cs;
	delete[] vs;
	
	//delete[] mv;
	return;
};

