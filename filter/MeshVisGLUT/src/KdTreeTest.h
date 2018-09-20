#ifndef KD_TREE_TEST_H
#define KD_TREE_TEST_H

#include "KdTree.h"
#include <boost/shared_ptr.hpp>

/*
 * Class to test the Kd-Tree implementation
 *
 */

#define NUM_TESTS 5

class KdTreeTest : public KdTree
{
public:
	enum Output {
		Informal,
		Statistical,
		Debug
	};

	KdTreeTest() : KdTree() {};
	~KdTreeTest() {};

	/* Point location Methods */
	int findPointNextToEstimate(Point p);
	int findPointNextToRandomized(Point p);
	int findPointNextToBacktracking(Point p);

	void runTest(Output o);

private:

	int findNearestNeighborEstimate(Point p, Point& nn, float& dist, int& cycles);
	int findNearestNeighborSecondPass(Point p, Point& nn, float& dist, int& cycles);
	int findNearestNeighborBacktracking(Point p, Point& nn, float& dist, int& cycles);

	/* Iterate over all nodes and find the nearest neighbor with "brute force" */
	void findNearestNeighborBF(Point p, Point& nn, float& dist, int& cycles); 

	int findNearestNeighborRandomized(Point p, Point& nn, float& dist, int& cycles, float c=2.0f, int iterations=10);
};

typedef boost::shared_ptr<KdTreeTest> KdTreeTestPtr;

#endif

