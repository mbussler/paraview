#pragma once

/**
 *
 *	\brief VTK Loader
 *  
 *  Load VTK Meshes
 *
 *	\author Michael Bu�ler
 */

#ifndef VTK_LOADER_H
#define VTK_LOADER_H

// definition of points, cells, etc
#include "Typedefs.h"

// needed to update the meshes' bounding box
#include "BoundingBox.h"

// the definition of a Tetrahedral Mesh
#include "TetMesh.h"

// the definition of an unstructured grid of the VTK
#include <vtkUnstructuredGrid.h>

#include "MeshLoader.h"

class VTKLoader : public MeshLoader
{

public:

	/// constructor
	VTKLoader();

	/// destructor
	~VTKLoader();

	/// load nodes from file
	void loadFromFile( const char* filename, TetMeshPtr mesh );

	// *** Synthetic dataset creation ***
	void setDomain( BoundingBox domain) { m_domain = domain;}
	void setSynthFlowParameter( float* _f) { a=_f[0]; b=_f[1]; c=_f[2]; };

	void createRandomPointSet(TetMeshPtr mesh, int num_points, bool save = true);
	void createTetrahedralGrid(TetMeshPtr mesh, int gridX, int gridY, int gridZ, float variance, bool save_to_file=false);

private:

	void calculateNeighbors( vtkUnstructuredGrid *grid, TetMeshPtr mesh );
	void VtkGridToMesh( vtkUnstructuredGrid* vtkGrid, TetMeshPtr mesh );

	BoundingBox m_domain;

	/// Synthetic flow functions
	REAL flow_u(REAL x, REAL y, REAL z);
	REAL flow_v(REAL x, REAL y, REAL z);
	REAL flow_w(REAL x, REAL y, REAL z);
	Point flow(REAL x, REAL y, REAL z);

	// parameter for the synthetic flow field
	float a,b,c;

};

#endif
