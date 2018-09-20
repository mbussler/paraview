#include "KdTree.h"

// file I/O
#include <iostream>
#include <fstream>
using namespace std;

/// constructor
KdTree::KdTree() 
{
	points = (Node*) NULL;
	in_nodes = (Vertex*) NULL;
	S = (REAL*) NULL;
	D = (char*) NULL;
    L = (int*) NULL;
    I = (int*) NULL;
	ptrD = 0;
	ptrS = 0;
	ptrL = 0;
	point_count = 0;
	nodes_count = 0;
	levels = 0;

    lookupTable=0;
};

/// destructor
KdTree::~KdTree()
{
	delete[] S; S=0;
	delete[] D; D=0;
	delete[] L; L=0;
	delete[] I; I=0;
	delete[] points; points=0;
};

void KdTree::createKdTree(TetMeshPtr mesh)
{
	// store references to nodes and length for traversal
	in_nodes = mesh->nodes;
	nodes_count = mesh->node_count;

    if( mesh->filename[0] == 0 )
    {
		loadNodes( mesh );
		buildKdTree();
    } 
    else if( !ReadKdTreeFromFile( mesh->filename) )
	{
		loadNodes( mesh );
		buildKdTree();
		WriteKdTreeToFile( mesh->filename );
	}
};

bool KdTree::ReadKdTreeFromFile( char* filename )
{
	char fn[255];
	sprintf(fn,"%s.kdt",filename);

	printf("[KdTree] Reading Kd-Tree from %s..", fn);

	ifstream kdFile ( fn, ios::binary );

	if( kdFile.is_open() )
	{
		kdFile.seekg(0);
		kdFile.read((char*)&point_count, sizeof(int));

		levels = log((double)point_count) / log(2.0);

		/* create arrays */
		S = new float[ point_count ];
		D = new char[ point_count ];
		L = new int[ point_count ];
        I = new int[ point_count ];
			
		/* read binary */
		kdFile.read((char*)S, point_count * sizeof(float));
		kdFile.read((char*)D, point_count * sizeof(char));
		kdFile.read((char*)L, point_count * sizeof(int));
		kdFile.read((char*)I, point_count * sizeof(int));

		kdFile.close();
		printf("done.\n");
		return( true );
	}
	else
	{
		printf("failed!.\n", fn);
		return false;
	}
};

void KdTree::WriteKdTreeToFile( char* filename )
{
	char fn[255];
	sprintf(fn,"%s.kdt",filename);

	printf("[KdTree] Writing Kd-Tree to %s..", fn);

	ofstream kdFile ( fn, ios::binary | ios::trunc );

	if( kdFile.is_open() )
	{
		kdFile.seekp(0);
		kdFile.write((char*)&point_count, sizeof(int));

		/* write binary */
		kdFile.write((char*)S, point_count * sizeof(float));
		kdFile.write((char*)D, point_count * sizeof(char));
		kdFile.write((char*)L, point_count * sizeof(int));
		kdFile.write((char*)I, point_count * sizeof(int));

		kdFile.close();
		printf("done.\n");
	}
	else
	{
		printf("error: could not open %s for write.\n", fn);
	}
};

void KdTree::loadNodes( TetMeshPtr mesh)
{
    // calculate and store point-to-cell lookup table
    MeshProcessorPtr mp( new MeshProcessor() );
    lookupTable = mp->calculateLookupTable( mesh );
    
	printf("[KdTree] \n");
	printf("[KdTree] *** loading nodes ***\n");
	printf("[KdTree] \n");
	printf("[KdTree] %d nodes to be loaded.\n", mesh->node_count);
	
	//calculate levels of final Kd-Tree
	levels = (int)ceil(log((double)nodes_count)/log(2.0f));

	point_count = (int)pow(2.0f, levels);
	points = new Node[point_count];

	// copy node values to point array
	for( int k=0; k<nodes_count; k++)
	{
		points[k].index = k;
		
		// copy values
		points[k].p	= mesh->nodes[k];

		start_box.updateBoundingBox(points[k].p);
	}

	// insert additional points at the end
	int add_pc = point_count - nodes_count;
	if( add_pc > 0)
	{
		double off = (nodes_count / (double) add_pc);
		for( int k=0; k<add_pc; k++)
		{
			points[nodes_count+k] = points[(int)floor(k*off)];
		}
	}


};


void KdTree::buildKdTree()
{
	// build a kd-tree of maximum size
	buildKdTree(1);
}

void KdTree::buildKdTree( int max_nodes )
{
	printf("[KdTree] \n");
	printf("[KdTree] *** Build up Kd-Tree ***\n");
	printf("[KdTree] \n");
		
	printf("[KdTree] Final Kd-Tree has %d levels.\n",levels);
	
	printf("[KdTree] Allocating memory...");

		S = new float[ point_count ];
		D = new char[ point_count ];
		L = new int[ point_count ];
		I = new int[ point_count ];
	
	printf("done!\n");

	printf("\n[KdTree] *** Memory Consumption: ***\n");
	printf("[KdTree] S: %d bytes\n", point_count * sizeof(float));
	printf("[KdTree] D: %d bytes\n", point_count * sizeof(char));
	printf("[KdTree] L: %d bytes\n\n", point_count * sizeof(int));

	printf("[KdTree] Perform recursive build-up...");

		root = buildKdTree(0, point_count-1, levels, start_box);
	
	printf("done!\n");
	
	printf("[KdTree] Creating Arrays...");
	
		// Traverse tree breadth-first order to store in arrays
		std::queue<TreeNode*> traverse_list;
		traverse_list.push(root);

		while( !traverse_list.empty() )
		{
			TreeNode* curr = traverse_list.front();
			if( curr->splitDim == -1) // leave node
			{
				storeAtEndOfLeaves(curr->index);
			}
			else
			{
				storeAtEndOfSplitDimArray(curr->splitDim);
				storeAtEndOfSplitValueArray(curr->splitValue);
				traverse_list.push(curr->leftNode);
				traverse_list.push(curr->rightNode);
			}
			/* free mem */
			delete curr;
			traverse_list.pop();
		}

	printf("done!\n");

};

TreeNode* KdTree::buildKdTree( int a, int b, int level, BoundingBox box)
{
	//printf("[KdTree] buildKdTree called with a=%d, b=%d, level=%d.\n",a,b,level);

	BoundingBox boundingBox, leftBox, rightBox;
	TreeNode* result = new TreeNode();
	
	/* store boundaries */
	result->a = a;
	result->b = b;

	if ((b-a) > 1000)
		boundingBox = BoundingBoxOfRandomPoints(a,b);
	else
		boundingBox = calculateBoundingBox(a,b);

	char splitDim = calculateSplitDimension(box, boundingBox);
	REAL splitValue = SplitPointsAtMedian(a,b,splitDim);

		//printf("[KdTree] Split dimension is dim=%d.\n",splitDim);
		//printf("[KdTree] Split value is %.2f.\n",splitValue);
	
	result->splitDim = splitDim;
	result->splitValue = splitValue;
	
	level--;
	if (level == 0)
	{
		// Store index in splitvalue of left successor node
		result->leftNode=new TreeNode();
		result->leftNode->index = points[a].index;
		result->leftNode->splitDim = -1;

		// Store index in splitvalue of right successor node
		result->rightNode=new TreeNode();
		result->rightNode->index = points[b].index;
		result->rightNode->splitDim = -1;

		//printf("[KdTree] Number of points in this leave: %d\n", b-a);
		//printf("[KdTree] Bounding Box:\n");
		//printf("[KdTree] min: (%f,%f,%f)\n",box.point_min.x,box.point_min.y,box.point_min.z);
		//printf("[KdTree] max: (%f,%f,%f)\n",box.point_max.x,box.point_max.y,box.point_max.z);
	} 
	else
	{
		leftBox = computeLeftBox( box, splitDim, splitValue);
		result->leftNode = 
			buildKdTree( a, (a+b)/2, level, leftBox);
		
		rightBox = computeRightBox( box, splitDim, splitValue);
		result->rightNode =
			buildKdTree( ((a+b)/2)+1, b, level, rightBox);
	}

	return result;
};

BoundingBox KdTree::BoundingBoxOfRandomPoints(int a, int b, int count)
{
	int span = b-a;
	BoundingBox result;

	for( int i=0; i<count; i++)
	{
		int index = (rand()%span + a);
		result.updateBoundingBox(points[index].p);
	};
	
	return result;
};

BoundingBox KdTree::calculateBoundingBox(int a, int b)
{
	BoundingBox result;
	for( int i=a; i<=b; i++)
	{
		result.updateBoundingBox(points[i].p);
	}
	return result;
};

char KdTree::calculateSplitDimension(BoundingBox box, BoundingBox boundingBox)
{
	char resultDim = 0;
	REAL max_length = 0.0f;

	for( char dim = 0; dim < 3; dim++)
	{
		REAL length = box.getLengthOfDim(dim) * boundingBox.getLengthOfDim(dim);
		
		if ( length > max_length)
		{
			max_length = length;
			resultDim = dim;
		}
	}
	return resultDim;
};

void KdTree::storeAtEndOfSplitDimArray(char splitDim)
{
	assert( ptrD < point_count );
	D[ptrD] = splitDim;
	ptrD++;
};

void KdTree::storeAtEndOfSplitValueArray(REAL splitValue)
{
	assert( ptrS < point_count );
	S[ptrS] = splitValue;
	ptrS++;
};

void KdTree::storeAtEndOfLeaves(int index)
{
	assert( ptrL < point_count );
    I[ptrL] = index;
    L[ptrL] = lookupTable[index];
	ptrL++;
};

REAL KdTree::SplitPointsAtMedian(int a, int b, char splitDim)
{
	/* printf("[KdTree] vor Sortierung f�r dim=%d:\n",splitDim); */
		// for( int i=a; i<=b;i++) printf("%.2f,",points[i].p.c[splitDim]);

	int median = quick_select(points, a, b, splitDim);

	REAL result=0;
	if ((a+b)%2 == 1)
	{
		result = (points[median].p.c[splitDim] + points[median+1].p.c[splitDim]) / 2.0f;
	}
	else
	{
		result = points[median].p.c[splitDim];
	}

	/* printf("\n[KdTree] nach Sortierung:\n"); */
		//for( int i=a; i<=b;i++) printf("%.2f,",points[i].p.c[splitDim]);
		//printf("\n[KdTree] Median ist %.2f",result);

	return result;
};

BoundingBox KdTree::computeLeftBox(BoundingBox box, char splitDim, REAL splitValue)
{
	BoundingBox result;
	
	// copy box coords
	result.point_min = box.point_min;
	result.point_max = box.point_max;

	// adjust coords of point_max
	assert(splitDim >= 0 && splitDim < 3);
	result.point_max.c[splitDim] = splitValue;

	return result;
};
BoundingBox KdTree::computeRightBox(BoundingBox box, char splitDim, REAL splitValue)
{
	BoundingBox result;
	
	// copy box coords
	result.point_min = box.point_min;
	result.point_max = box.point_max;

	// adjust coords of point_min
	assert(splitDim >= 0 && splitDim < 3);
	result.point_min.c[splitDim] = splitValue;

	return result;
};


Point KdTree::PointFromIndex( int index )
{
	if( index > point_count) {
		return Point();
	}

    return in_nodes[ index ];
};
