#pragma once

/**
 *
 *	\brief Mesh Loader
 *  
 *  This MeshLoader is the base class for other loaders
 *
 *	\author Michael Bussler
 */

#ifndef MESH_LOADER_H
#define MESH_LOADER_H

#include "Typedefs.h"
#include "TetMesh.h"
#include "BoundingBox.h"


class MeshLoader{

public:

	/// constructor
	MeshLoader();

	/// destructor
	~MeshLoader();

protected:

	/// load neighbors from file
	bool loadNeighborsFromFile( const char* filename, TetMeshPtr mesh);
	void saveNeighborsToFile( const char* filename, TetMeshPtr mesh );

};

#endif
