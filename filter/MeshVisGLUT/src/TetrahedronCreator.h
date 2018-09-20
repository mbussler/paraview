
#ifndef TETRAHEDRON_CREATOR_H
#define TETRAHEDRON_CREATOR_H

#include "Tetrahedron.h"
#include "TetMesh.h"

#include <boost/shared_ptr.hpp>

class TetrahedronCreator
{
public:
	TetrahedronCreator();
	~TetrahedronCreator();

	Tetrahedron* getTetrahedron( int index );
	void interpolateFieldValue( const int& cell, 
								const double& c1, 
								const double& c2, 
								const double& c3, 
								Vector& velocity);

	void setMesh( TetMeshPtr input_mesh );
	void init();

private:
	TetMeshPtr m_mesh;
	Tetrahedron** tetList;

};
typedef boost::shared_ptr<TetrahedronCreator> TetrahedronCreatorPtr;


#endif
