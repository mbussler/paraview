#include "TetMeshCuda.h"

TetMeshCuda::TetMeshCuda()
{
    pointsVBO = 0;    
    boxVBO = 0;       
    traversedCellsIBO = 0;
    traversedCellsVBO = 0;
    outerFacesIBO = 0;
	outerFaceNormalsVBO = 0;

    traversed_cells_count = 0;
	outer_faces_count = 0;

	hNodes = 0;
	hNodeAttributes = 0;
	hCells = 0;
	hNeighbors = 0;

	nodesVBOMapped = false;
};

TetMeshCuda::~TetMeshCuda()
{
    //freeGPUMem();
};

void TetMeshCuda::freeGPUMem()
{
#ifndef NOGL
    deleteVBOs();
#else
	freeArray( m_meshGPU.nodes );
#endif

	freeArray( m_meshGPU.cells );
    freeArray( m_meshGPU.neighbors );
    freeArray( m_meshGPU.nodeAttributes );
	freeArray( m_meshGPU.traversedCells);

	delete[] hTraversedCells;

	// free page-locked host memory
	if(hNodes)          freePageLockedHostMemory(hNodes);
    if(hNodeAttributes) freePageLockedHostMemory(hNodeAttributes);
	if(hCells)          freePageLockedHostMemory(hCells);
    if(hNeighbors)      freePageLockedHostMemory(hNeighbors);

	hNodes = 0;
	hNodeAttributes = 0;
	hCells = 0;
	hNeighbors = 0;

};

void TetMeshCuda::initVBOs()
{
    /* Create VBO for vertices and velocities*/
	createPointsVBO();

	/* Create VBO for the Bounding Box */
	glGenBuffers(1, &boxVBO);

	// VBOs for traversed cells rendering
	glGenBuffers(1, &traversedCellsVBO );
	glGenBuffers(1, &traversedCellsIBO );
	glGenBuffers(1, &outerFacesIBO);
	glGenBuffers(1, &outerFaceNormalsVBO);

	createBoxVBO();
	createOuterFacesIBO();
	createOuterFaceNormalsVBO();
};

void TetMeshCuda::createPointsVBO()
{
	unsigned int size = node_count * sizeof( float4 );

	pointsVBO = createVBO(size);
	registerGLBufferObject( pointsVBO, &pointsVBO_CUDA );
};

uint TetMeshCuda::createVBO(uint size)
{
	GLuint vbo;
	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, size, 0, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	return vbo;
};


void TetMeshCuda::deleteVBOs()
{
	// return if extensions were not loaded
	if( glIsBuffer == NULL)
		return;

	if( nodesVBOMapped )
		unmapVBOtoCuda();

	unregisterGLBufferObject( pointsVBO_CUDA );
	if( glIsBuffer(pointsVBO))
	{
		glBindBuffer(1, pointsVBO);
		glDeleteBuffers(1, (const GLuint*)&pointsVBO);
	}

	if( glIsBuffer(boxVBO))
	{
		glBindBuffer(1, boxVBO);
		glDeleteBuffers(1, &boxVBO);
	}
	if( glIsBuffer(traversedCellsVBO))
	{
		glBindBuffer(1, traversedCellsVBO);
		glDeleteBuffers(1, &traversedCellsVBO);
	}
	if( glIsBuffer(traversedCellsIBO))
	{
		glBindBuffer(1, traversedCellsIBO);
		glDeleteBuffers(1, &traversedCellsIBO);
	}
	if( glIsBuffer(outerFacesIBO))
	{
		glBindBuffer(1, outerFacesIBO);
		glDeleteBuffers(1, &outerFacesIBO);
	}
	if( glIsBuffer(outerFaceNormalsVBO))
	{
		glBindBuffer(1, outerFaceNormalsVBO);
		glDeleteBuffers(1, &outerFaceNormalsVBO);
	}
	pointsVBO   = 0;
	boxVBO      = 0;
	traversedCellsVBO  = 0;
	traversedCellsIBO  = 0;
	outerFacesIBO      = 0;
	outerFaceNormalsVBO = 0;
};

MeshGPU& TetMeshCuda::mapVBOtoCuda()
{
	if( !nodesVBOMapped)
	{
#ifndef NOGL
		m_meshGPU.nodes = (float4*) mapGLBufferObject( &pointsVBO_CUDA );
#endif
		nodesVBOMapped = true;
	}

	return m_meshGPU;
};

void TetMeshCuda::unmapVBOtoCuda()
{
	if( nodesVBOMapped )
	{
#ifndef NOGL
		unmapGLBufferObject( pointsVBO_CUDA );
#endif
		nodesVBOMapped = false;
	}
};

void TetMeshCuda::copyDataToGPU()
{
	m_meshGPU.s_nodes = sizeof( float4 ) * node_count;
    m_meshGPU.s_cells = sizeof( int4 ) * cells_count;
    m_meshGPU.num_cells = cells_count;

#ifndef NOGL
	initVBOs();
	mapVBOtoCuda();
#else
    allocateArray( (void**)&m_meshGPU.nodes, m_meshGPU.s_nodes );
#endif

	// alloc mem: nodes and node attributes
	copyArrayToDevice( m_meshGPU.nodes, nodes, m_meshGPU.s_nodes );

    allocateArray( (void**)&m_meshGPU.nodeAttributes, m_meshGPU.s_nodes );
    copyArrayToDevice( m_meshGPU.nodeAttributes, node_attributes, m_meshGPU.s_nodes );

    // alloc mem: cells and neighbors
    allocateArray( (void**)&m_meshGPU.cells, m_meshGPU.s_cells );
    copyArrayToDevice( m_meshGPU.cells, (int4*) cells, m_meshGPU.s_cells );

    allocateArray( (void**)&m_meshGPU.neighbors, m_meshGPU.s_cells );
    copyArrayToDevice( m_meshGPU.neighbors, (int4*) neighbors, m_meshGPU.s_cells );

	// alloc mem: traversed cells
	allocateArray( (void**)&m_meshGPU.traversedCells, cells_count * sizeof( char ) );
	hTraversedCells = new char[ cells_count ];
    resetTraversedCells();

#ifndef NOGL
	unmapVBOtoCuda();
#endif

};

void TetMeshCuda::copyDataToGPUAsync(cudaStream_t stream)
{
	// calculate memory usage
    int nSize = sizeof( float4 ) * node_count;
    int cSize = sizeof( int4 ) * cells_count;

    m_meshGPU.s_nodes = nSize;
    m_meshGPU.s_cells = cSize;
    m_meshGPU.num_cells = cells_count;

#ifndef NOGL
    initVBOs();
    mapVBOtoCuda();
#else
	allocateArray( (void**)&m_meshGPU.nodes, nSize );
#endif

    allocateArray( (void**)&m_meshGPU.nodeAttributes, nSize );
    allocateArray( (void**)&m_meshGPU.cells, cSize );
    allocateArray( (void**)&m_meshGPU.neighbors, cSize );

	bool writeCombined = true;
	allocatePageLockedArrayPortable( (void**)&hNodes, nSize, writeCombined );
    allocatePageLockedArrayPortable( (void**)&hNodeAttributes, nSize, writeCombined );
    allocatePageLockedArrayPortable( (void**)&hCells, cSize, writeCombined );
    allocatePageLockedArrayPortable( (void**)&hNeighbors, cSize, writeCombined );

	// allocate page-locked host memory and copy asynchronously to device
    copyArrayToPageLockedHostMemory( hNodes, this->nodes, nSize);
	copyArrayToPageLockedHostMemory( hNodeAttributes, this->node_attributes, nSize);
	copyArrayToPageLockedHostMemory( hCells, this->cells, cSize);
	copyArrayToPageLockedHostMemory( hNeighbors, this->neighbors, cSize);

	printf("Copy grid: "); startTimer(stream);

	copyArrayToDeviceAsync( m_meshGPU.nodes, hNodes, nSize, stream);
	copyArrayToDeviceAsync( m_meshGPU.nodeAttributes, hNodeAttributes, nSize, stream);
    copyArrayToDeviceAsync( m_meshGPU.cells, hCells, cSize, stream);
    copyArrayToDeviceAsync( m_meshGPU.neighbors, hNeighbors, cSize, stream);

    stopTimer(stream); printTimer(); printf("\n"); destroyTimer(stream);

    // alloc mem: traversed cells
	hTraversedCells = new char[ cells_count ];
	allocateArray( (void**)&m_meshGPU.traversedCells, cells_count * sizeof( char ) );
    resetTraversedCells();

#ifndef NOGL
	//unmapVBOtoCuda();
#endif
};


void TetMeshCuda::resetTraversedCells()
{
	// reset mem for traversed cells
	setArray( m_meshGPU.traversedCells, 0, cells_count * sizeof( char));
};

void TetMeshCuda::PerformTetWalk( float4* points, int numPoints, int* cells, int* dBlocks, cudaStream_t stream)
{
	mapVBOtoCuda();
	PerformTetWalkCuda( points, numPoints, cells, dBlocks, m_meshGPU, stream );
	unmapVBOtoCuda();
};


void TetMeshCuda::draw( DWORD flags)
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
        glBindBuffer(GL_ARRAY_BUFFER, boxVBO);
		glVertexPointer(3, GL_FLOAT, 0, BUFFER_OFFSET(0));

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		//glColor3f(1.0, 1.0, 1.0);
		glColor3f(.1, .1, .1);
		glDrawArrays( GL_QUADS, 0, 24);
	}

	// use the points-vbo from now on
    glBindBuffer(GL_ARRAY_BUFFER, pointsVBO);

	// set the vertex-pointer to node-coordinates
	glVertexPointer(3, GL_FLOAT, sizeof(float4), BUFFER_OFFSET(0));

	if( checkFlag(flags, DRAW_VERTICES))
	{
		// draw Mesh Nodes
		glColor3f(1.0, 0.0, 0.0);
		glPointSize(2.0f);
        glDrawArrays(GL_POINTS, 0, node_count);
	}

	if( checkFlag(flags, DRAW_OUTER_FACES))
	{
		
// Wireframe rendering
/*
		glBindBuffer(GL_ARRAY_BUFFER, outerFacesIBO);
		glVertexPointer(3, GL_FLOAT, 0, BUFFER_OFFSET(0));

		glEnable(GL_LINE_SMOOTH);
		glLineWidth(.5);

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		//glColor3f(0.4, 0.4, .4);
		glColor3f(0.6, 0.6, .6);
		glDrawArrays( GL_TRIANGLES, 0, 3 * outer_faces_count);
*/
// Nice rendering with lighting

		glEnable(GL_LIGHTING);
		glEnable(GL_NORMALIZE);
		glEnable(GL_CULL_FACE);
		glEnable(GL_LINE_SMOOTH);
		glEnableClientState(GL_NORMAL_ARRAY);

		glBindBuffer(GL_ARRAY_BUFFER, outerFacesIBO);
		glVertexPointer(3, GL_FLOAT, 0, BUFFER_OFFSET(0));

		glBindBuffer(GL_ARRAY_BUFFER, outerFaceNormalsVBO);
		glNormalPointer( GL_FLOAT, 0, BUFFER_OFFSET(0));

		glShadeModel( GL_SMOOTH );
		glLineWidth(0.8f);

		// Draw Front Faces as solid, smooth shaded polygons
		glColor3f(1.0f, 0.676f, 0.272f);

		glPolygonMode(GL_FRONT, GL_FILL);
		glCullFace(GL_BACK);

		glDepthRange(0.01, 1.0);
		glDrawArrays( GL_TRIANGLES, 0, 3 * outer_faces_count);

		glDisable(GL_LIGHTING);
		
		// Draw Back Faces as wireframe
		glColor3f(0.2, 0.2, .2);

		glPolygonMode(GL_FRONT, GL_LINE);
		glCullFace(GL_BACK);

		glDepthRange(0.0, 1.0);
		glDrawArrays( GL_TRIANGLES, 0, 3 * outer_faces_count);

		glDisable(GL_CULL_FACE);
		glDisable(GL_NORMALIZE);
		glDisableClientState(GL_NORMAL_ARRAY);

	}

	if( checkFlag(flags, DRAW_TRAVERSED_CELLS))
	{
        glBindBuffer(GL_ARRAY_BUFFER, traversedCellsVBO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, traversedCellsIBO);

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
};

void TetMeshCuda::createTraversedCellsVBO( const bool &use_modelview, const glh::matrix4 &m)
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

	copyArrayFromDevice( hTraversedCells, m_meshGPU.traversedCells, cells_count * sizeof( char));

	int tc_count = 0;

	// Back to Front Ordering
	list<pair<int, float> > toSort;

	for( int i=0; i < cells_count; i++)
	{
		if( hTraversedCells[i] == 1 )
		{
			tc_count++;

			pair<int,float> item;
			item.first = i;

			// calculate single point for current cell (e.g. center) and transform by MODELVIEW_MATRIX
			Vector vs[4];
			for( int j=0; j<4; j++) {
				vs[j] = Pt2Vec( nodes[ cells[i].indices[j] ] );
			}

			// Calculate Cell Center
			Vector c = vs[0] + 0.5 * (vs[1]-vs[0]) + 0.5 * (vs[2]-vs[0]) + 0.5 * (vs[3]-vs[0]);
	
			// Depth Value
			item.second = c.x*mv[2] + c.y*mv[6] + c.z*mv[10];

			toSort.push_back(item);
		}
	}
	toSort.sort(compare_depth);

	traversed_cells_count = tc_count;

	// stores vertex indices, 12 per cell
	int num_cells = 12 * tc_count;
	unsigned int* cs = new unsigned int[num_cells];

	// store vertices and vertex colors: 2*4*4 floats per cell
	int num_vert = 2*4*4 * traversed_cells_count;
	GLfloat* vs = new GLfloat[num_vert];

	Color3f cellColor;
	float alpha = 0.3;

	list< pair< int,float > >::iterator c_it;
	int idx=0;

	for( c_it=toSort.begin(); c_it!=toSort.end(); c_it++)
	{
		int cell_idx = c_it->first;

		cellColor = tc_cellcolor[ cell_idx % tc_cellcolor_count];

		/* iterate over the four cell vertices */ 
		for( int i=0; i<4; i++) 
		{
			/* v_idx is the index of the i'th Vertex of cell *c_it */ 
			int v_idx = cells[cell_idx].indices[i];

			/* store components of the current vertex with index v_idx */
            vs[32*idx +8*i +0] = nodes[v_idx].x;
			vs[32*idx +8*i +1] = nodes[v_idx].y;
			vs[32*idx +8*i +2] = nodes[v_idx].z;

			// store color components for the current vertex
    		vs[32*idx +8*i +4] = cellColor.r;
    		vs[32*idx +8*i +5] = cellColor.g;
    		vs[32*idx +8*i +6] = cellColor.b;
			vs[32*idx +8*i +7] = alpha;

			/* store indice j for the current triangle i */
			cs[12*idx +3*i +0] = 4*idx + triangles[i][0];
			cs[12*idx +3*i +1] = 4*idx + triangles[i][1];
			cs[12*idx +3*i +2] = 4*idx + triangles[i][2];

		}
		idx++;
	};

	/* fill VBO */
	glBindBuffer(GL_ARRAY_BUFFER, traversedCellsVBO);
	glBufferData(GL_ARRAY_BUFFER, num_vert* sizeof(GLfloat), vs, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	/* fill IBO */
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, traversedCellsIBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, num_cells* sizeof(int), cs, GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	delete[] cs;
	delete[] vs;
	
	//delete[] mv;
	return;
};

void TetMeshCuda::createBoxVBO()
{
    // box vertices
	vector<Point> pts;	
	box.getVertices( pts );

	const int box_indices[24] = { 0,1,2,3, 4,7,6,5, 0,4,5,1, 3,2,6,7, 0,3,7,4, 1,5,6,2 };

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

	glBindBuffer(GL_ARRAY_BUFFER, boxVBO);
	glBufferData(GL_ARRAY_BUFFER, size * sizeof(GLfloat), data, GL_STATIC_DRAW);
	
	/* cleanup */
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	delete data;
};

void TetMeshCuda::createOuterFacesIBO()
{
	//if( !outer_faces ) return;

	//// a face is drawn as triangle
	//int size = 3 * outer_faces_count * sizeof(int);

	///* fill IBO for outer face rendering */
	//glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, outerFacesIBO);
	//glBufferData(GL_ELEMENT_ARRAY_BUFFER, size, (int*)outer_faces, GL_STATIC_DRAW);
	//glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	if( !face_vertices ) return;

	// a face is drawn as triangle
	int size = 3 * outer_faces_count * sizeof(Vector);

	/* fill IBO for outer face rendering */
    glBindBuffer(GL_ARRAY_BUFFER, outerFacesIBO);
	glBufferData(GL_ARRAY_BUFFER, size, (GLfloat*)face_vertices , GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
};

void TetMeshCuda::createOuterFaceNormalsVBO()
{
	if( !outer_face_normals ) return;

	// the three vertices of each face have the same normal
	int size = 3 * outer_faces_count * sizeof(Vector);

	/* fill IBO for outer face rendering */
    glBindBuffer(GL_ARRAY_BUFFER, outerFaceNormalsVBO);
	glBufferData(GL_ARRAY_BUFFER, size, (GLfloat*)outer_face_normals, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
};

