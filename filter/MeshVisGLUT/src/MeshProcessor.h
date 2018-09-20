/**
 *
 *	\brief Mesh Processor
 *
 *	The Mesh Processor class, used to calculate adjacent cells and the lookup table.
 *  Adjacent cell calculation is done by using the tetgen library.
 *
 *	\author Michael Buﬂler
 */

#include "TetMesh.h"

#include <boost/shared_ptr.hpp>

#ifndef MESH_PROCESSOR_H
#define MESH_PROCESSOR_H

class MeshProcessor{

public:

	/// constructor
	MeshProcessor();

	/// destructor
	~MeshProcessor();

	/// calculate lookup table
	int* calculateLookupTable( TetMeshPtr mesh );

private:
	
};
typedef boost::shared_ptr<MeshProcessor> MeshProcessorPtr;


#endif

