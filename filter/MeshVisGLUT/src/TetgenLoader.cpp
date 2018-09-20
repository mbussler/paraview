#include "TetgenLoader.h"
#include "InlineFunctions.h"
#include "Matrix.h"
#include "vector.h"

#include "tetgen.h"

#ifdef _MEM_LEAK_DETECT
#define new new(_NORMAL_BLOCK, __FILE__, __LINE__)
#endif

/// constructor
TetgenLoader::TetgenLoader()
{
	a=0.0f;
	b=1.0f;
	c=0.1f;
	m_domain = BoundingBox( Point( -1, -1, -1), Point( 1, 1, 1));
};

/// destructor
TetgenLoader::~TetgenLoader()
{
};

/// load nodes from file
void TetgenLoader::loadNodesFromFile(char* filename, TetMeshPtr mesh)
{
	printf("[TetgenLoader] Loading nodes from %s.node.\n", filename);

	tetgenio tetgen_mesh;

	// load mesh with tetgen and copy nodes, cells and attributes to tetmesh
	if (tetgen_mesh.load_node(filename))
	{
		copyNodesFromTetgen( mesh, tetgen_mesh );
		printf("[TetgenLoader] %d nodes loaded.\n", mesh->node_count);
		printf("[TetgenLoader] %d attributes per node.\n", mesh->attribute_count);
	}
};

/// load simulation data from file
void TetgenLoader::loadMeshFromFiles(char* base_filename, TetMeshPtr mesh)
{
	sprintf( mesh->filename, "%s", base_filename);

	printf("[TetgenLoader] Loading with tetgen from %s.*.\n", base_filename);

	tetgenio tetMesh;

	if (tetMesh.load_tetmesh(base_filename))
	{
		copyNodesFromTetgen(mesh, tetMesh);
		copyCellsFromTetgen(mesh, tetMesh);

		if( !loadNeighborsFromFile( base_filename, mesh ) )
		{
			tetrahedralize("rn", &tetMesh, &tetMesh);
			copyNeighborsFromTetgen(mesh, tetMesh);
			saveNeighborsToFile( base_filename, mesh );
		}

		printf("[TetgenLoader] %d nodes loaded.\n", mesh->node_count);
		printf("[TetgenLoader] %d attributes per node.\n", mesh->attribute_count);
		printf("[TetgenLoader] %d cells loaded.\n", mesh->cells_count);
	}
};

void TetgenLoader::saveAsTetgen(char* filename, TetMeshPtr mesh)
{
	printf("[TetgenLoader] Saving mesh with tetgen to %s.*.\n", filename);

	tetgenio tetgen_in;
	tetgenio tetgen_out;

	tetgen_in.numberofpoints = mesh->node_count;
	tetgen_in.numberoftetrahedra = mesh->cells_count;

	tetgen_in.pointlist = new REAL[3*mesh->node_count];
	tetgen_in.tetrahedronlist = new int[4*mesh->cells_count];

	Matrix m = MatrixRotationY(PI/4.0f);
	//m = m * MatrixRotationZ(PI/4.0f);
	Vector v;

	for(int i=0; i<mesh->node_count; i++)
	{
		float4 node = mesh->nodes[i];
		v.x = node.x;
		v.y = node.y;
		v.z = node.z;
		v = m*v;
		tetgen_in.pointlist[3*i+0] = (REAL) v.x;
		tetgen_in.pointlist[3*i+1] = (REAL) v.y;
		tetgen_in.pointlist[3*i+2] = (REAL) v.z;
	}

	for(int i=0; i<mesh->cells_count; i++)
	{
		Cell cell = mesh->cells[i];
		tetgen_in.tetrahedronlist[4*i+0] = cell.indices[0];
		tetgen_in.tetrahedronlist[4*i+1] = cell.indices[1];
		tetgen_in.tetrahedronlist[4*i+2] = cell.indices[2];
		tetgen_in.tetrahedronlist[4*i+3] = cell.indices[3];
	}

	printf("[TetgenLoader] Reconstructing Mesh...\n");
	tetrahedralize("r", &tetgen_in, &tetgen_out);

	printf("[TetgenLoader] Writing to files...\n");

	// Output mesh to files
	tetgen_out.save_nodes(filename);
	tetgen_out.save_faces(filename);
	tetgen_out.save_elements(filename);
};


void TetgenLoader::createRandomPointSet(TetMeshPtr mesh, int num_points, bool save)
{
	delete[] mesh->nodes;
	delete[] mesh->node_attributes;

	/// the Tetgen Mesh used to tetraheralize the domain
	tetgenio tetgen_in;

	tetgen_in.numberofpoints = num_points;
	tetgen_in.pointlist = new REAL[tetgen_in.numberofpoints * 3]; // 3 coords per point

	tetgen_in.numberofpointattributes = 3;
	tetgen_in.pointattributelist = new REAL[tetgen_in.numberofpoints * tetgen_in.numberofpointattributes];

	// create random point set and calculate flow field
	for( int i=0; i<num_points; i++)
	{
		Point rnd = m_domain.RandomPointInBox();

		tetgen_in.pointlist[3*i+0] = rnd.x;
		tetgen_in.pointlist[3*i+1] = rnd.y;
		tetgen_in.pointlist[3*i+2] = rnd.z;
		
		Point fl = flow( rnd.x, rnd.y, rnd.z );

		tetgen_in.pointattributelist[3*i+0] = fl.x;
		tetgen_in.pointattributelist[3*i+1] = fl.y;
		tetgen_in.pointattributelist[3*i+2] = fl.z;
	}

	char filename[50];
	sprintf(filename, "RandomPointSet_%d_a%4.2f_b%4.2f_c%4.2f.in", num_points, a,b,c);

	tetgenio tetgen_out;

	printf("[TetgenLoader] Calling TetGen to create Tetrahedral Mesh.\n");

#ifdef _WIN32
	tetrahedralize("n", &tetgen_in, &tetgen_out);
#else
	tetrahedralize("", &tetgen_in, &tetgen_in);
	tetrahedralize("rn", &tetgen_in, &tetgen_out);
#endif

	if ( save )
	{
		// save mesh to files
		sprintf(filename, "RandomPointSet_%d_a%4.2f_b%4.2f_c%4.2f",num_points, a,b,c);
		tetgen_out.save_nodes(filename);
		tetgen_out.save_faces(filename);
		tetgen_out.save_elements(filename);
	}

	mesh->updateBoundingBox( false );
	mesh->box = m_domain;
	mesh->grid = Point( 1, 1, 1 );

	// copy values
	copyNodesFromTetgen( mesh, tetgen_out );
	copyCellsFromTetgen( mesh, tetgen_out );
	copyNeighborsFromTetgen( mesh, tetgen_out );

	if( save )
	{
		saveNeighborsToFile( filename, mesh );
	}

	printf("[TetgenLoader] %d random nodes created.\n",num_points);

};


void TetgenLoader::createTetrahedralGrid( TetMeshPtr mesh, int gridX, int gridY, int gridZ, float variance, bool save_to_file)
{
	printf("[TetgenLoader] Create Tetrahedral Grid of %dx%dx%d.\n",gridX, gridY, gridZ);

	/// the Tetgen Mesh used to tetraheralize the domain
	tetgenio tetgen_in;

	tetgen_in.numberofpoints = gridX*gridY*gridZ;
	tetgen_in.pointlist = new REAL[tetgen_in.numberofpoints * 3]; // 3 coords per point

	tetgen_in.numberofpointattributes = 3;
	tetgen_in.pointattributelist = new REAL[tetgen_in.numberofpoints * tetgen_in.numberofpointattributes];

	int index = 0;
	int i,j,k;
	
	for( k = 0; k < gridZ; k++)
	{
		for( j = 0; j < gridY; j++)
		{
			for( i = 0; i < gridX; i++)
			{
				//Point p = Point(i,j,k);
				Point p;
				p.x = ((REAL)i / (REAL)(gridX-1)) * m_domain.getLengthOfDim(0) + m_domain.point_min.x;
				p.y = ((REAL)j / (REAL)(gridY-1)) * m_domain.getLengthOfDim(1) + m_domain.point_min.y;
				p.z = ((REAL)k / (REAL)(gridZ-1)) * m_domain.getLengthOfDim(2) + m_domain.point_min.z;

				// create random offset
				p.x += ((i==0) || (i==gridX-1)) ? 0 : Random(-variance, variance);
				p.y += ((j==0) || (j==gridY-1)) ? 0 : Random(-variance, variance);
				p.z += ((k==0) || (k==gridZ-1)) ? 0 : Random(-variance, variance);

				// calculate synthetic flow
				Point f = flow( p.x, p.y, p.z);
				
				for( int r=0; r<3;r++) {
					tetgen_in.pointlist[3*index+r]		  = p.c[r];
					tetgen_in.pointattributelist[3*index+r] = f.c[r];
				}
				index++;
			}
		}
	}
	
	// Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
	//   do quality mesh generation (q) with a specified quality bound
	//   (1.414), and apply a maximum volume constraint (a0.1).
	//tetrahedralize("q1.414a0.1", &in, &out);

	tetgenio tetgen_out;

	printf("[TetgenLoader] Calling TetGen to create Tetrahedral Mesh.\n");

#ifdef _WIN32
	tetrahedralize("n", &tetgen_in, &tetgen_out);
#else
	tetrahedralize("", &tetgen_in, &tetgen_in);
	tetrahedralize("rn", &tetgen_in, &tetgen_out);
#endif

	char filename[50];
	sprintf(filename, "TetMesh_%dx%dx%d_a%4.2f_b%4.2f_c%4.2f",gridX, gridY, gridZ, a,b,c);

	if( save_to_file )
	{
		// Output mesh to files
		tetgen_out.save_nodes(filename);
		tetgen_out.save_faces(filename);
		tetgen_out.save_elements(filename);
	}

	mesh->updateBoundingBox(false);
	mesh->box = m_domain;
	mesh->grid = Point( gridX, gridY, gridZ );

	// copy values
	copyNodesFromTetgen( mesh, tetgen_out);
	copyCellsFromTetgen( mesh, tetgen_out);
	copyNeighborsFromTetgen( mesh, tetgen_out);

	if( save_to_file )
	{
		saveNeighborsToFile( filename, mesh);
	}
};

// copy vertices and vertex attributes from tetgen
void TetgenLoader::copyNodesFromTetgen( TetMeshPtr mesh, const tetgenio& tetgenMesh )
{
	mesh->node_count		= tetgenMesh.numberofpoints;
	mesh->attribute_count	= 3; //tetgenMesh.numberofpointattributes;

	if( mesh->nodes) 
		delete[] mesh->nodes;
	mesh->nodes = new Vertex[mesh->node_count];

	if( mesh->node_attributes) 
		delete[] mesh->node_attributes;
	mesh->node_attributes	= new VertexAttribute[mesh->node_count];

	for( int i=0; i<mesh->node_count; i++)
	{
		Vertex v;
		v.x = tetgenMesh.pointlist[3*i+0];
		v.y = tetgenMesh.pointlist[3*i+1];
		v.z = tetgenMesh.pointlist[3*i+2];
		v.w = 1.0;

		mesh->nodes[i] = v;

		if( mesh->getUpdateBoundingBox() )
			mesh->box.updateBoundingBox( v );

		VertexAttribute va;
		va.x = tetgenMesh.pointattributelist[3*i+0];
		va.y = tetgenMesh.pointattributelist[3*i+1];
		va.z = tetgenMesh.pointattributelist[3*i+2];
		va.w = 1.0;

		mesh->node_attributes[i] = va;
	}
};
void TetgenLoader::copyCellsFromTetgen( TetMeshPtr mesh, const tetgenio& tetgenMesh )
{
	mesh->cells_count	= tetgenMesh.numberoftetrahedra;

	if( mesh->cells )
		delete[] mesh->cells;
	mesh->cells = new Cell[mesh->cells_count];

	// copy cells from tetgen
	for( int i=0; i<mesh->cells_count; i++ )
	{
		Cell c; 
		for( int j=0; j<4; j++) 
			c.indices[j] = tetgenMesh.tetrahedronlist[4*i+j];

		mesh->cells[i] = c;
	}
};

void TetgenLoader::copyNeighborsFromTetgen( TetMeshPtr mesh, const tetgenio& tetgenMesh )
{
	if( mesh->neighbors)
		delete[] mesh->neighbors;
	mesh->neighbors = new Cell[mesh->cells_count];

	// copy neighbors from tetgen
	for( int i=0; i<mesh->cells_count; i++ )
	{
		Cell n; 
		for( int j=0; j<4; j++)
		{
			n.indices[j] = tetgenMesh.neighborlist[4*i+j];
		}
		mesh->neighbors[i] = n;
	}
};

bool TetgenLoader::copyMeshToTetgen( TetMeshPtr mesh, tetgenio& tetgenMesh )
{
	delete tetgenMesh.pointlist;
	delete tetgenMesh.pointattributelist;
	delete tetgenMesh.tetrahedronlist;

	tetgenMesh.numberofpoints			= mesh->node_count;
	tetgenMesh.numberofpointattributes	= mesh->attribute_count;
	tetgenMesh.numberoftetrahedra		= mesh->cells_count;

	// alloc mem
	tetgenMesh.pointlist			= new REAL[mesh->node_count * 3];
	tetgenMesh.pointattributelist	= new REAL[mesh->node_count * mesh->attribute_count];;
	tetgenMesh.tetrahedronlist		= new int[ mesh->cells_count * 4 ];

	// copy vertices and vertex attributes to tetgen
	for( int i= 0; i< mesh->node_count; i++ )
	{
		tetgenMesh.pointlist[3*i+0] = mesh->nodes[i].x;
		tetgenMesh.pointlist[3*i+1] = mesh->nodes[i].y;
		tetgenMesh.pointlist[3*i+2] = mesh->nodes[i].z;

		tetgenMesh.pointattributelist[3*i+0] = mesh->node_attributes[i].x;
		tetgenMesh.pointattributelist[3*i+1] = mesh->node_attributes[i].y;
		tetgenMesh.pointattributelist[3*i+2] = mesh->node_attributes[i].z;

		if( mesh->getUpdateBoundingBox() )
			mesh->box.updateBoundingBox( mesh->nodes[i] );
	}

	// copy cells to tetgen
	for( int i= 0; i< mesh->cells_count; i++ )
	{
		tetgenMesh.tetrahedronlist[4*i+0] = mesh->cells[i].indices[0];
		tetgenMesh.tetrahedronlist[4*i+1] = mesh->cells[i].indices[1];
		tetgenMesh.tetrahedronlist[4*i+2] = mesh->cells[i].indices[2];
		tetgenMesh.tetrahedronlist[4*i+3] = mesh->cells[i].indices[3];
	}

	return true;
};

bool TetgenLoader::copyNeighborsToTetgen( TetMeshPtr mesh, tetgenio& tetgenMesh )
{
	assert(	tetgenMesh.tetrahedronlist && 
		tetgenMesh.numberoftetrahedra	>0 );

	delete[] tetgenMesh.neighborlist;

	// four neighbours per cell, non existing cells are marked as -1
	tetgenMesh.neighborlist = new int[ 4* mesh->cells_count ];

	// copy neighbors to tetgen
	for( int i= 0; i< mesh->cells_count; i++ )
	{
		tetgenMesh.neighborlist[4*i+0] = mesh->neighbors[i].indices[0];
		tetgenMesh.neighborlist[4*i+1] = mesh->neighbors[i].indices[1];
		tetgenMesh.neighborlist[4*i+2] = mesh->neighbors[i].indices[2];
		tetgenMesh.neighborlist[4*i+3] = mesh->neighbors[i].indices[3];
	};

	return true;
};


REAL TetgenLoader::flow_u(REAL x, REAL y, REAL z)
{
	return a*x-b*y;
};

REAL TetgenLoader::flow_v(REAL x, REAL y, REAL z)
{
	return b*x+a*y;
};

REAL TetgenLoader::flow_w(REAL x, REAL y, REAL z)
{
	return -2*a*z+c;
};

Point TetgenLoader::flow(REAL x, REAL y, REAL z)
{
	Point p;
	p.x = flow_u(x,y,z);
	p.y = flow_v(x,y,z);
	p.z = flow_w(x,y,z);
	return p;
};
