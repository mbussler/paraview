#include "MeshProcessor.h"

/// constructor
MeshProcessor::MeshProcessor()
{
};

/// destructor
MeshProcessor::~MeshProcessor()
{
};

/// calculate lookup table
int* MeshProcessor::calculateLookupTable( TetMeshPtr mesh )
{
	// create references
	Vertex* nodes	= mesh->nodes;
	Cell* cells		= mesh->cells;
	int num_nodes	= mesh->node_count;
	int num_cells	= mesh->cells_count;

	// create arrays
	int* N		= new int[num_nodes];
	int* count	= new int[num_nodes];

	// initialize
	for( int i=0; i<num_nodes; i++) {
		N[i]=0; count[i]=0;
	}

	/* Fill the counter array: for each Vertex count the number of cells */
	/* iterate over all cells */
	for( int ci=0; ci<num_cells; ci++) 
	{
		/* iterate over all point indices of cell c_i */
		for( int pi=0; pi<4; pi++)
		{
			/* read index of pi's point of cell ci and increment counter for this point */
			int point = cells[ci].indices[pi]; 
			assert( point < num_nodes && point >= 0 );
			count[ point ]++;
		}
	}

	/* calculate the offset for each point and store it in array N */
	int current_offset=0;
	for( int i=0; i<num_nodes; i++) 
	{
		N[i] = current_offset;
        
        // HACK: remove points that don't belong to a cell
        if( count[i] == 0){
          nodes[i].x = nodes[i].y = nodes[i].z = nodes[i].w = -1000000000.0f;          
        } else {
          current_offset += count[i];
          /* reset the counter */
          count[i]=0;
        }
	}

	/* current_offset is now the sum of the counters which is the length of C */ 
	int length_of_C = current_offset;

	/* initialize C */
	int* C = new int[length_of_C];


	/* Iterate over all cells again but this time cell indices are stored in array C */
	for( int ci=0; ci<num_cells; ci++) 
	{
		for( int pi=0; pi<4; pi++) 
		{
			int point = cells[ci].indices[pi]; 
			assert( point < num_nodes && point >= 0 );
			
			/* store cell index ci in C using offset N[point] and cell counter count[point] */
			C[ N[point] + count[point]  ] = ci;

			/* increment counter */
			count[point]++;
		}
	}

    int* res = new int[ num_nodes ];
    for( int i=0; i<num_nodes; i++ )
    {
        res[i] = C[N[i]];
    }

	delete[] N;
    delete[] C;
    delete[] count;

    return res;
};

