#include "KdTreeTest.h"
#include "InlineFunctions.h"

#include <stdio.h>

int KdTreeTest::findPointNextToEstimate(Point p)
{
	Point nn;
	float dist;
	int cycles;

	//printf("\n[KdTree] *** Tree-traversal for (%.2f,%.2f,%.2f). ***\n",p.x,p.y,p.z);

	//int point_index = findNearestNeighborRandomized(p, nn, dist, cycles);
	//int point_index = findNearestNeighborBacktracking(p, nn, dist, cycles);
	int point_index = findNearestNeighborEstimate(p, nn, dist, cycles);
	
	//printf("[KdTree] Nearest Point found at (%.2f,%.2f,%.2f) with a distance of %.2f @ %d cycles.\n",
	//	nn.x,nn.y,nn.z, dist, cycles);

	return point_index;
};

int KdTreeTest::findPointNextToRandomized(Point p)
{
	Point nn;
	float dist;
	int cycles;

	printf("\n[KdTree] *** Tree-traversal for (%.2f,%.2f,%.2f). ***\n",p.x,p.y,p.z);

	int point_index = findNearestNeighborRandomized(p, nn, dist, cycles);
	
	printf("[KdTree] Nearest Point found at (%.2f,%.2f,%.2f) with a distance of %.2f @ %d cycles.\n",
		nn.x,nn.y,nn.z, dist, cycles);

	return point_index;
};

int KdTreeTest::findPointNextToBacktracking(Point p)
{
	Point nn;
	float dist;
	int cycles;

	//printf("\n[KdTree] *** Tree-traversal for (%.2f,%.2f,%.2f). ***\n",p.x,p.y,p.z);

	int point_index = findNearestNeighborBacktracking(p, nn, dist, cycles);
	
	printf("[KdTree] Nearest Point found at (%.2f,%.2f,%.2f) with a distance of %.2f @ %d cycles.\n",
		nn.x,nn.y,nn.z, dist, cycles);

	return point_index;
};




void KdTreeTest::runTest(Output output)
{
	const int   test_run_count = 1000;
	const float rnd_param_c    = 2.0f;
	const int   rnd_iterations = 10;

	switch( output ) {
		case Debug:
		case Informal:
			printf("\n[KdTree] Running Kd-Tree Test\n");
		break;
		case Statistical:
			printf("SP:dist\tSP:cycles\t2P:dist\t2P:cycles\tBT:dist\tBT:cycles\tRND:dist\tRND:cycles\n");
		break;
	}

	/* statistics */
	int* min_c = new int[NUM_TESTS];
	int* max_c = new int[NUM_TESTS];
	int* sum_c = new int[NUM_TESTS];

	double* min_d = new double[NUM_TESTS];
	double* max_d = new double[NUM_TESTS];
	double* sum_d = new double[NUM_TESTS];

	int* fail = new int[NUM_TESTS];
	int random_within_bounds=0;

	/* init */
	for (int i=0; i<NUM_TESTS; i++) {
		min_c[i]=1000000; max_c[i]=0; sum_c[i]=0;
		min_d[i]=100000000; max_d[i]=0; sum_d[i]=0;
		fail[i]=0;
	}

	Point* nn = new Point[NUM_TESTS];
	float* dist = new float[NUM_TESTS];
	int* cycles = new int[NUM_TESTS];

	for( int i=0; i<test_run_count; i++)
	{
		REAL x = Random() * start_box.getLengthOfDim(0) + start_box.point_min.c[0];
		REAL y = Random() * start_box.getLengthOfDim(1) + start_box.point_min.c[1];
		REAL z = Random() * start_box.getLengthOfDim(2) + start_box.point_min.c[2];

		Point p = Point(x,y,z);
		
		findNearestNeighborEstimate(p, nn[0], dist[0], cycles[0]);
		findNearestNeighborSecondPass(p, nn[1], dist[1], cycles[1]);
		findNearestNeighborBacktracking (p, nn[2], dist[2], cycles[2]);
		findNearestNeighborRandomized (p, nn[3], dist[3], cycles[3], rnd_param_c, rnd_iterations);
		findNearestNeighborBF(p, nn[4], dist[4], cycles[4]);

		/* update statistics */
		for( int i=0; i<NUM_TESTS; i++)
		{
			if( dist[i] < min_d[i]) min_d[i] = dist[i];	
			if( dist[i] > max_d[i]) max_d[i] = dist[i];	
			if( cycles[i] < min_c[i]) min_c[i] = cycles[i]; 
			if( cycles[i] > max_c[i]) max_c[i] = cycles[i]; 
			sum_d[i] += dist[i];
			sum_c[i] += cycles[i];
		}
		if( nn[0] != nn[4] ) { 
			fail[0]++;
		}
		if( nn[1] != nn[4] ) { 
			fail[1]++;
			if( output == Debug ) {
				printf("*2P found pt @dist=%.4f,\n",dist[2]);
				printf(" but best is @dist=%.4f!\n",dist[4]);
			}
		}
		if( nn[2] != nn[4] ) { 
			fail[2]++;
			if( output == Debug ) {
				printf("*BT found pt @dist=%.4f,\n",dist[3]);
				printf(" but best is @dist=%.4f!\n",dist[4]);
			}
		}
		if( nn[3] != nn[4] ) {
			fail[3]++;
			if( dist[3] > rnd_param_c * dist[4])
				random_within_bounds++;
		}

		switch ( output ) {
			case Debug:
				printf("\n[KdTree] *** Tree-traversal for (%.2f,%.2f,%.2f). ***\n",p.x,p.y,p.z);
				printf("[KdTree] Single Pass:  (%.2f,%.2f,%.2f) with a distance of %.2f @ %d cycles.\n",nn[0].x,nn[0].y,nn[0].z, dist[0], cycles[0]);
				printf("[KdTree] Second Pass:  (%.2f,%.2f,%.2f) with a distance of %.2f @ %d cycles.\n",nn[1].x,nn[1].y,nn[1].z, dist[1], cycles[1]);
				printf("[KdTree] Backtracking: (%.2f,%.2f,%.2f) with a distance of %.2f @ %d cycles.\n",nn[2].x,nn[2].y,nn[2].z, dist[2], cycles[2]);
				printf("[KdTree] Randomized:   (%.2f,%.2f,%.2f) with a distance of %.2f @ %d cycles.\n",nn[3].x,nn[3].y,nn[3].z, dist[3], cycles[3]);
				printf("[KdTree] Brute Force : (%.2f,%.2f,%.2f) with a distance of %.2f @ %d cycles.\n",nn[4].x,nn[4].y,nn[4].z, dist[4], cycles[4]);
			break;
			case Statistical:
				printf("%.4f\t%d\t",dist[0], cycles[0]);
				printf("%.4f\t%d\t",dist[1], cycles[1]);
				printf("%.4f\t%d\t",dist[2], cycles[2]);
				printf("%.4f\t%d\n",dist[3], cycles[3]);
			break;
		}
	}
	
	if( output == Debug || output == Informal )
	{
		printf("[KdTree] *** Statistics ***\n");
		printf("\n");
		printf("[KdTree] * Single Pass * \n");
		printf("[KdTree]   Distance: min=%.2f, max=%.2f, avg=%.2f \n",min_d[0], max_d[0], (sum_d[0] / (double)test_run_count));
		printf("[KdTree]   Cycles: min=%d, max=%d, avg=%d \n",min_c[0], max_c[0], (sum_c[0] / test_run_count));
		printf("\n");
		printf("[KdTree] * Second Pass * \n");
		printf("[KdTree]   Distance: min=%.4f, max=%.4f, avg=%.4f \n",min_d[1], max_d[1], (sum_d[1] / (double)test_run_count));
		printf("[KdTree]   Cycles: min=%d, max=%d, avg=%d \n",min_c[1], max_c[1], (sum_c[1] / test_run_count));
		printf("\n");
		printf("[KdTree] * Backtracking * \n");
		printf("[KdTree]   Distance: min=%.4f, max=%.4f, avg=%.4f \n",min_d[2], max_d[2], (sum_d[2] / (double)test_run_count));
		printf("[KdTree]   Cycles: min=%d, max=%d, avg=%d \n",min_c[2], max_c[2], (sum_c[2] / test_run_count));
		printf("\n");
		printf("[KdTree] * Randomized * \n");
		printf("[KdTree]   Distance: min=%.4f, max=%.4f, avg=%.4f \n",min_d[3], max_d[3], (sum_d[3] / (double)test_run_count));
		printf("[KdTree]   Cycles: min=%d, max=%d, avg=%d \n",min_c[3], max_c[3], (sum_c[3] / test_run_count));
		printf("\n");
		printf("[KdTree] Single Pass failed %d times.\n",fail[0]);
		printf("[KdTree] Second Pass failed %d times.\n",fail[1]);
		printf("[KdTree] Backtracking failed %d times.\n",fail[2]);
		printf("[KdTree] Randomized failed %d times but failed for the %.1f-best point only %d times.\n",fail[3], rnd_param_c, random_within_bounds);
	}

};





void KdTreeTest::findNearestNeighborBF(Point p, Point& nn, float& dist, int& cycles)
{
	dist = 1e10;

	for( int i=0; i<nodes_count; i++)
	{
		Point r = in_nodes[i];
		REAL dx = r.x - p.x;
		REAL dy = r.y - p.y;
		REAL dz = r.z - p.z;

		REAL dist_sq = (dx*dx)+(dy*dy)+(dz*dz);

		if ( dist_sq < dist)
		{
			dist = dist_sq;
			nn = r;
		}
	}

	cycles = nodes_count;
	dist = sqrt( dist);
	//printf("[KdTree]     Best Point for (%.2f,%.2f,%.2f): ",p.x,p.y,p.z);
	//printf("(%.2f,%.2f,%.2f) dist=%.2f\n",best.x,best.y,best.z,sqrt(dist));
};


int KdTreeTest::findNearestNeighborEstimate(Point p, Point& nn, float& dist, int& cycles)
{
	int l=0;
	int index=0;
	
	// Traverse tree to find the right index to access index-array L
	while (l < levels)
	{
		char splitDim = D[index];
		REAL splitValue = S[index];
		
		int succ = (p.c[splitDim] < splitValue) ? 1 : 2;
		index = 2*index + succ;

		l++;
	}
	
	// calculate offset to access L
	int leave = index-point_count+1;
	
	nn = PointFromIndex( I[leave] );
	dist = distanceToPoint(p,nn);
	cycles = l;

	return L[leave];
};

/* this method uses a two-pass aproach and a queue to find the nearest neighbor */
int KdTreeTest::findNearestNeighborSecondPass(Point p, Point& _nn, float& _dist, int& _cycles)
{
	/* first pass */

	int l=0;
	int index=0;
	int cycles=0;

	int nn_idx=0;

	/* Traverse tree to find the current best distance to 
	   the nearest neighbor in the current bucket */
	while (l < levels)
	{
		char splitDim = D[index];
		REAL splitValue = S[index];
		
		int succ = (p.c[splitDim] < splitValue) ? 1 : 2;
		index = 2*index + succ;

		l++;
		cycles++;
	}
	
	// calculate offset to access L
	int idx = I[(index-point_count+1)];

	Point r = PointFromIndex(idx);

	float current_best_distance = distanceToPoint(p,r);
	nn_idx = idx;

	/* second pass with boundary check*/

	Point best=r;

	// a queue is used to emulate the recursion
	std::queue<int> indices;
	indices.push(0);

	while( !indices.empty() )
	{
		index = indices.front();
		while (index < (point_count-1))
		{
			cycles++;
			char splitDim = D[index];
			REAL splitValue = S[index];
			
			if (p.c[splitDim] < splitValue)
			{
				index = 2*index + 1;
				if( p.c[splitDim] + current_best_distance > splitValue)
				{
					/* also test the other successor */
					indices.push(index+1);
				}
			}
			else
			{
				index = 2*index + 2;
				if( p.c[splitDim] - current_best_distance < splitValue)
				{
					/* also test the other successor */
					indices.push(index-1);
				}
			}
		}
		
		/* update current-best value */
		int idx = I[(index-point_count+1)];
		
		Point r = PointFromIndex(idx);
		float new_dist = distanceToPoint(p,r);

		if( new_dist < current_best_distance)
		{
			current_best_distance = new_dist;
			best = r;
			nn_idx = idx;
			//printf("[KdTree] Better Point at distance %.2f found after %d cycles.\n",new_dist,cycles);
		}
		indices.pop();
	}

	_cycles = cycles;
	_nn = best;
	_dist = current_best_distance;

	return nn_idx;
};



/*
 * 
 *
 */
int KdTreeTest::findNearestNeighborBacktracking(Point p, Point& _nn, float& _dist, int& _cycles)
{
	/* reference to the current processed node */
	int index=0;
	
	/* number of computing cycles, e.g. node processing steps */
	int cycles=0;

	/* the distance to the current best neighbor of p */
	float current_best_dist=1e10;

	/* the nearest neighbor of p */
	Point NearestNeighbor;
	int nn_idx=0;

	/* recursion is emulated by two stacks:	
	 * the indices stack collects the traversed tree nodes for back tracking */
	std::stack<int> indices;

	/* the alternatives stack collects the (sub)-root nodes of the alternative paths to traverse */
	std::queue<int> alt;

	/* the first traversal path  starts at the root node at index 0 */
	alt.push(0);

	/* process all alternatives */
	while( !alt.empty() )
	{
		index=alt.front(); alt.pop();
		//printf("root: %d\n",index);

		/* Traverse the tree starting at the current root node to 
		 * find the current best distance to the nearest neighbor in the current bucket. 
		 * The visited nodes are collected in the indices stack */

		//printf("[stack] ");
		while( index < (point_count-1))
		{
			indices.push(index); 
			//printf("%d ",index);
			
			char splitDim = D[index];
			REAL splitValue = S[index];
			
			/* choose successor node and add offset to index */
			int succ = (p.c[splitDim] < splitValue) ? 1 : 2;
			index = 2*index + succ;

			cycles++;
		}
		//printf(": %d\n",indices.size());	

		/* calculate the current best distance */
		int point_idx = I[(index-point_count+1)];
		Point r = PointFromIndex( point_idx );
		float dist = distanceToPoint(p, r);

		if( dist < current_best_dist)
		{
			current_best_dist = dist;
			NearestNeighbor = r;
			nn_idx = point_idx;
			//printf("new best distance: %.2f\n",dist);
		}

		/* the traversed nodes are visited again in reversed order to 
		 * find buckets which may include a nearer neighbor */

		while( !indices.empty() )
		{
			index = indices.top(); indices.pop();

			char splitDim = D[index];
			REAL splitValue = S[index];
			
			if( p.c[splitDim] < splitValue)
			{
				if( p.c[splitDim] + current_best_dist > splitValue) 
				{
					//alt.push( 2*index+2 );
					int idx = 2*index+2;
					alt.push( idx );
					//printf(" %d added to alternatives\n",idx);
				}
			}
			else
			{
				if( p.c[splitDim] - current_best_dist < splitValue) 
				{
					//alt.push( 2*index+1 );
					int idx = 2*index+1;
					alt.push( idx );
					//printf(" %d added to alternatives\n",idx);
				}
			}
			cycles++;
		}
	}

	_nn = NearestNeighbor;
	_dist = current_best_dist;
	_cycles = cycles;

	return nn_idx;
};


int KdTreeTest::findNearestNeighborRandomized(Point p, Point& _nn, float& _dist, int& _cycles, float c, int iterations)
{
	Point nn;
	float dist;
	int cycles;
	int nn_i=0;
	
	// first kd tree traversal to estimate distance
	nn_i = findNearestNeighborEstimate(p, nn, dist, cycles);
	Point rnd_nn;
	float rnd_dist;
	int rnd_cycles;

	// traverse with random points
	for( int i=0; i<iterations; i++) 
	{
		Point rnd = RandomPointOnSphere(p, dist/2.0f);
		int rnd_i = findNearestNeighborEstimate( rnd, rnd_nn, rnd_dist, rnd_cycles);
		
		rnd_dist = distanceToPoint(p, rnd_nn);

		if( rnd_dist < dist) {
			nn=rnd_nn;
			nn_i=rnd_i;
			dist=rnd_dist;
		}
		cycles += rnd_cycles;
	}

	_nn = nn;
	_dist = dist;
	_cycles = cycles;

	return nn_i;
};
