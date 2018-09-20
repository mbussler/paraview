
#pragma once

#ifndef TETWALKER_H
#define TETWALKER_H

#include "KdTreeTest.h"
#include "MeshProcessor.h"
#include "TetrahedronCreator.h"

#include <list>
#include <boost/dynamic_bitset.hpp>
#include <boost/shared_ptr.hpp>

class TetWalker
{
public:

	/// constructor
	TetWalker();
	/// destructor
	~TetWalker();

	/// The TetWalker needs a reference to a KD-Tree implementation for global searching
	void setKDTree(KdTreeTestPtr kdt) {_kdt = kdt;};

	/// The TetWalker needs a reference to a TetrahedronCreator implementation
	void setTetrahedronCreator( TetrahedronCreatorPtr tc) {_tc=tc;};
	
	/// Perform a global cellsearch for the given point
	int locateCell(const Point& p);

	/// calculate the interpolated velocity for a given position with a given start-cell for the TetWalk
	void interpolateVelocityAt( const Point& position, int& startCell, Vector& velocity);

	/// Give the TetWalker a reference to an int-list to store traversed cells at
	void setupVis( boost::dynamic_bitset<>* vis_traversed_cells) { _traversed_cells = vis_traversed_cells; };

private:

	// Tetrahedral Walk implementation
	bool performTetWalk( const Point& p, int& startCell, double& c1, double& c2, double& c3);

	/// reference to the KD-Tree
	KdTreeTestPtr _kdt;

	/// reference to the TetrahedronCreator
	TetrahedronCreatorPtr _tc;

	/// the list of traversed cells
	boost::dynamic_bitset<>* _traversed_cells;
};

typedef boost::shared_ptr<TetWalker> TetWalkerPtr;

#endif
