#include "TetMesh.h"

TetMesh::TetMesh() 
{
	nodes=0;
	node_attributes=0;
	cells=0;
	neighbors=0;
	node_count=0;
	cells_count=0;
	attribute_count=0;
	outer_faces=0;
	outer_face_normals=0;
	face_vertices=0;
	m_updateBoundingBox=true;
    
    filename[0] = 0;
	
	grid = Point( 1, 1, 1 ); // unsteady
};

TetMesh::~TetMesh() 
{
	delete[] nodes;
	delete[] node_attributes;
	delete[] cells;
	delete[] neighbors;

	delete[] outer_faces;
	delete[] outer_face_normals;
};

void TetMesh::copyOuterFaces(const list<Face>& fs)
{
	// copy list to array
	if( outer_faces)
		delete[] outer_faces;

	outer_faces_count = fs.size();
	outer_faces = new Face[outer_faces_count];

	int idx = 0;
	list<Face>::const_iterator f = fs.begin();
	while( f!=fs.end() )
	{
		outer_faces[idx].pts[0] = f->pts[0];
		outer_faces[idx].pts[1] = f->pts[1];
		outer_faces[idx].pts[2] = f->pts[2];
		idx++;
		f++;
	}
};

void TetMesh::copyNormals(const list<Vector>& ns)
{
	// copy list to array
	if( outer_face_normals )
		delete[] outer_face_normals;

	outer_face_normals = new Vector[ 3*outer_faces_count ];

	int idx = 0;
	list<Vector>::const_iterator n = ns.begin();
	
	while( n!=ns.end() )
	{
		outer_face_normals[idx++] = *n;
		outer_face_normals[idx++] = *n;
		outer_face_normals[idx++] = *n;
		n++;
	}
};

void TetMesh::copyFaceVertices(const list<Vector>& vs)
{
	// copy list to array
	if( face_vertices)
		delete[] face_vertices;

	face_vertices = new Vector[ 3*outer_faces_count ];

	int idx = 0;
	list<Vector>::const_iterator v = vs.begin();
	
	while( v != vs.end() )
	{
		face_vertices[idx++] = *v; v++;
		face_vertices[idx++] = *v; v++;
		face_vertices[idx++] = *v; v++;
	}
};
