
#include "TetrahedronCreator.h"
#include "InlineFunctions.h"

TetrahedronCreator::TetrahedronCreator()
{
	tetList=0;
}; 

TetrahedronCreator::~TetrahedronCreator()
{
	if( tetList) {
		for( int i=0; i<m_mesh->cells_count; i++) {
			delete tetList[i];  } }

	delete[] tetList;
};

void TetrahedronCreator::setMesh( TetMeshPtr input_mesh ) 
{ 
	m_mesh = input_mesh;
};

void TetrahedronCreator::init() 
{
	tetList = new Tetrahedron* [m_mesh->cells_count];
	for( int i=0; i < m_mesh->cells_count; i++)
	{
		tetList[i] = NULL;
	}
}


Tetrahedron* TetrahedronCreator::getTetrahedron( int index )
{
	// test bounds
	if( index < 0 || index > m_mesh->cells_count) 
		return NULL;
	
	// create new tetrahedron if not in buffer
	if ( tetList[index] == NULL )
	{
		tetList[index] = new Tetrahedron(index);

		Vertex ps[4];
		int ns[4];

		// get physical coordinates and neighbors
		for( int i=0; i<4; i++)
		{
			int point_index = m_mesh->cells[index].indices[i];
			ps[i] = m_mesh->nodes[point_index];
			ns[i] = m_mesh->neighbors[index].indices[i];
		}
		tetList[index]->setPhysicalCoordinates(ps[0], ps[1], ps[2], ps[3]);
		tetList[index]->setNeighbors(ns[0], ns[1], ns[2], ns[3]);
	}

	return tetList[index];
};

void TetrahedronCreator::interpolateFieldValue(const int& cell, const double& c1, const double& c2, const double& c3, Vector& velocity)
{
	Vector vs[4];
	for( int i=0; i<4; i++)
	{
		int point_index = m_mesh->cells[cell].indices[i];
		vs[i] = Pt2Vec( m_mesh->node_attributes[point_index] );
	}

	velocity = vs[0] + (vs[1]-vs[0])*c1 + (vs[2]-vs[0])*c2 + (vs[3]-vs[0])*c3;
};
