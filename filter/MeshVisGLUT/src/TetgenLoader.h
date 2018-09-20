#pragma once

/**
 *
 *	\brief Mesh Loader
 *  
 *  This TetgenLoader uses the Tetgen Library to load nodes and meshes.
 *
 *	\author Michael Bu�ler
 */

#ifndef TETGEN_LOADER_H
#define TETGEN_LOADER_H

#include "tetgen.h"
#include "Typedefs.h"
#include "BoundingBox.h"
#include "TetMesh.h"
#include "MeshLoader.h"

#include <boost/shared_ptr.hpp>

class TetgenLoader : public MeshLoader
{

public:

	/// constructor
	TetgenLoader();

	/// destructor
	~TetgenLoader();

	/// load nodes from file
	void loadNodesFromFile(char* filename, TetMeshPtr mesh);

	/// load nodes and cells from files
	void loadMeshFromFiles(char* base_filename, TetMeshPtr mesh);

	void saveAsTetgen(char* filename, TetMeshPtr mesh);

	// *** Synthetic dataset creation ***
	void createRandomPointSet(TetMeshPtr mesh, int num_points, bool save = true);
	void createTetrahedralGrid(TetMeshPtr mesh, int gridX, int gridY, int gridZ, float variance, bool save_to_file=false);

	void copyNodesFromTetgen( TetMeshPtr mesh, const tetgenio& tetgenMesh );
	void copyCellsFromTetgen( TetMeshPtr mesh, const tetgenio& tetgenMesh );
	void copyNeighborsFromTetgen( TetMeshPtr mesh, const tetgenio& tetgenMesh );

	bool copyMeshToTetgen( TetMeshPtr mesh, tetgenio& tetgenMesh );
	bool copyNeighborsToTetgen( TetMeshPtr mesh, tetgenio& tetgenMesh );

	void setDomain( BoundingBox domain) { m_domain = domain;}

	// parameter for the synthetic flow field
	float a,b,c;

	void setSynthFlowParameter( float* _f)
	{ 
		//assert( _f[0] && _f[1] && _f[2] ); 
		
		a=_f[0]; 
		b=_f[1]; 
		c=_f[2]; 
	};

private:

	BoundingBox m_domain;

	/// Synthetic flow funtions
	REAL flow_u(REAL x, REAL y, REAL z);
	REAL flow_v(REAL x, REAL y, REAL z);
	REAL flow_w(REAL x, REAL y, REAL z);
	Point flow(REAL x, REAL y, REAL z);
};

typedef boost::shared_ptr<TetgenLoader> TetgenLoaderPtr;

#endif
