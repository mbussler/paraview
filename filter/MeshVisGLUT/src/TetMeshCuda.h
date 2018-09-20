
#pragma once

#ifndef TETMESH_CUDA_H
#define TETMESH_CUDA_H

#include "TetMeshGL.h"
#include <boost/shared_ptr.hpp>

// includes, CUDA
#include "common.cuh"

// includes for Tetwalk
#include "Tetwalk.cuh"

using namespace std;

/**
 *	\brief TetMeshCuda class
 *  
 *	Extension of TetMesh with CUDA support
 *
 *	\author Michael Bussler
*/


/// TetMeshCu class definition
class TetMeshCuda : public TetMesh
{
public:

	TetMeshCuda();
	~TetMeshCuda();

   	void draw( DWORD flags);

    // Asynchronous host to device copy using page-locked memory
    void copyDataToGPUAsync(cudaStream_t stream);

    void copyDataToGPU();
    void freeGPUMem();


	void PerformTetWalk( float4* points, int numPoints, int* cells, int* dBlocks, cudaStream_t stream);

	MeshGPU& mapVBOtoCuda();
	void unmapVBOtoCuda();

	void createTraversedCellsVBO( const bool &use_modelview = false, const glh::matrix4 &m = glh::matrix4());
	void resetTraversedCells();

	void initVBOs();

private:

    MeshGPU m_meshGPU;

    // pointers to page-locked host memory
	float4* hNodes;
	float4* hNodeAttributes;
	int4* hCells;
	int4* hNeighbors;

    // vbo variables
	GLuint	pointsVBO,    // Holds point coordinates
			boxVBO,       // Box indices
			traversedCellsIBO,   // Hold traversed cell indices
			traversedCellsVBO,   // VBO for traversed cells
			outerFacesIBO,       // Hold the indices of the outer faces
			outerFaceNormalsVBO; // Hold the normals of the outer faces

    struct cudaGraphicsResource* pointsVBO_CUDA;

    bool nodesVBOMapped;

	// VBO functions
    void deleteVBOs();
    void createBoxVBO();
	void createOuterFacesIBO();
	void createOuterFaceNormalsVBO();
	void createPointsVBO();

	uint createVBO(uint size);

	// counter
	int traversed_cells_count;

    // list of traversed cells on host
    char* hTraversedCells;

};
typedef boost::shared_ptr<TetMeshCuda> TetMeshCudaPtr;


#endif
