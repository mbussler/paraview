
#pragma once

#ifndef TETMESH_H
#define TETMESH_H

#include "Typedefs.h"
#include "BoundingBox.h"

#include <boost/shared_ptr.hpp>

// CUDA Vector types
#include <vector_types.h>


/**
 *	\brief TetMesh class
 *  This class defines a Mesh consisting of tetrahedron cells.
 *  The Mesh is described by an array of nodes and an array of cell indices.
 *
 *	\author Michael Bussler
*/

typedef float4	Vertex;
typedef float4	VertexAttribute;
struct	Cell	{ int indices[4]; };
struct  Face	{ int pts[3]; }; // a face consists of three point indices

// Reihenfolge, in der die faces eines tetrahedron iteriert werden fuer die nachbarschaftssuche
const int faces[4][3] = { {1,2,3}, {0,2,3}, {0,1,3}, {0,1,2}};

/// Mesh type definition
class TetMesh {

public:

	TetMesh();
	~TetMesh();

	/// the nodes of the mesh
	Vertex* nodes;
	
	/// the list of attributes per vertex
	VertexAttribute* node_attributes;

	/// the number of nodes of the loaded mesh
	int node_count;

	/// number of attributes per node
	int attribute_count;

	/// the cells of the loaded mesh
	Cell* cells;

	/// the number of cells of the loaded mesh
	int cells_count;

	Cell* neighbors;

	void copyOuterFaces(const list<Face>& fs);
	void copyNormals(const list<Vector>& ns);
	void copyFaceVertices(const list<Vector>& ns);

	BoundingBox box;
	Point grid;
	char filename[255];

	void updateBoundingBox( bool update )
	{ m_updateBoundingBox = update; }

	bool getUpdateBoundingBox()
	{ return m_updateBoundingBox; }

protected:
	/// The boundary mesh defined as faces
	Face* outer_faces;
	Vector* outer_face_normals;
	Vector* face_vertices;
	int outer_faces_count;

private:
	bool m_updateBoundingBox;

};

typedef boost::shared_ptr<TetMesh> TetMeshPtr;

#endif
