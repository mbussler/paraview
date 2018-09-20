#include "TetWalker.h"

TetWalker::TetWalker() 
{
	_traversed_cells = 0;
};

TetWalker::~TetWalker() 
{
};

int TetWalker::locateCell(const Point& p)
{
	if( !_kdt || !_tc) {
		return -1;
	}

	/* find a grid cell next to the point */
	int current_cell = _kdt->findPointNextToEstimate(p);
	
	assert( current_cell != -1 );

	/* natural coordinates of p */
	double c1,c2,c3;
	performTetWalk( p, current_cell, c1, c2, c3 );

	return current_cell;
};

// Tetrahedral Walk implementation
bool TetWalker::performTetWalk( const Point& p, int& startCell, double& c1, double& c2, double& c3)
{
	if ( !_tc ) return false;

	int nextCell = startCell;
	
	int count = 0;
	const int maxHops = 20;

	do  
	{ 
		// store the current Cell
		startCell = nextCell;

		// get the current Tetrahedron
		Tetrahedron* tetrahedron = _tc->getTetrahedron( startCell );

		// No neighboring cell exists -> return false
		if( tetrahedron == NULL ) {
			return false;
		}

		// store the cell index
		assert( startCell != -1);
		if( _traversed_cells ) _traversed_cells->set( startCell );

		// Calculating the natural coordinates of p gives the next cell to test
		tetrahedron->calculateNaturalCoordinates(p, c1, c2, c3, nextCell);
		
		if( count > maxHops )
		{
			startCell = -1; // HACK: kill this particle		
			return false;
		} 

		count++;
	} 
	while (  startCell != nextCell);
	
	return true;
};


void TetWalker::interpolateVelocityAt( const Point& position, int& startCell, Vector& velocity)
{
	if( startCell == -1 ) {
		velocity = Vector(0.0f, 0.0f, 0.0f);
		return;
	}

	double c1,c2,c3;
	
	bool success = performTetWalk( position, startCell, c1, c2, c3 );

	if( success ) {
		// calculate interpolated Velocity
		_tc->interpolateFieldValue( startCell, c1, c2, c3, velocity );
	} else {
		velocity = Vector(0.0f, 0.0f, 0.0f);
	}
};
