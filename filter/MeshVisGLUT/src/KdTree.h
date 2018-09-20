/**
 *
 *	\brief The kd-tree class
 *  
 *  
 *
 *	\author Michael Bussler
 */

#ifndef KD_TREE_H
#define KD_TREE_H

#include "Typedefs.h"
#include "BoundingBox.h"
#include "quickselect.h"
#include "Vector.h"
#include "TetMesh.h"
#include <queue>
#include <stack>

#include <boost/shared_ptr.hpp>

// calculate Point-to-Cell lookup table
#include "MeshProcessor.h"

class KdTree
{

public:

	/// constructor
	KdTree();

	/// destructor
	~KdTree();

	void createKdTree(TetMeshPtr mesh);

public:

	/// the loaded points to sort into the kd tree
	Node* points;

	/// The split values
	float* S;

	/// The dimension values
	char* D;

	/// The cell by pointindex
	int* L;

    /// The Point indices
	int* I;

	/// The number of levels the final Kd-Tree has
	int levels;

	/// the number of loaded points
	int point_count;

protected:

	/// load nodes and create an intern copy
	void loadNodes(TetMeshPtr mesh);

	bool ReadKdTreeFromFile( char* filename );
	void WriteKdTreeToFile( char* filename );

	/// create a KdTree from the given nodes and store result in S,D,L
	void buildKdTree();

	/// create a KdTree with a maximum of max_nodes nodes per leave
	void buildKdTree( int max_nodes);

	/// get the point at index
	Point PointFromIndex( int index );

	/// the input nodes
	const Vertex* in_nodes;
	
	/// the number of input nodes
	int nodes_count;

	/// The box to start the calculation with
	BoundingBox start_box;

private:
	
	/// Build the Kd Tree from points[a..b] with depth levels and associated bounding box 
	TreeNode* buildKdTree( int a, int b, int levels, BoundingBox box);

	BoundingBox BoundingBoxOfRandomPoints(int a, int b, int count=1000);
	BoundingBox calculateBoundingBox(int a, int b);
	char calculateSplitDimension(BoundingBox box, BoundingBox boundingBox);

	void storeAtEndOfSplitDimArray(char splitDim);
	void storeAtEndOfSplitValueArray(float splitValue);
	void storeAtEndOfLeaves(int index);

	BoundingBox computeLeftBox(BoundingBox box, char splitDim, REAL splitValue);
	BoundingBox computeRightBox(BoundingBox box, char splitDim, REAL splitValue);

	REAL SplitPointsAtMedian(int a, int b, char splitDim);

	int ptrS;
	int ptrD;
	int ptrL;

	TreeNode* root;

    int* lookupTable;

};

typedef boost::shared_ptr<KdTree> KdTreePtr;

#endif
