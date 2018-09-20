

#ifndef MESH_FUNCTIONS_CU
#define MESH_FUNCTIONS_CU

// USE_TEX
#include "common.cuh"

// textures for cells etc.
texture<float4, 1, cudaReadModeElementType> nodes1Tex;
texture<float4, 1, cudaReadModeElementType> nodeAttributes1Tex;
texture<int4,   1, cudaReadModeElementType> cells1Tex;
texture<int4,   1, cudaReadModeElementType> neighbors1Tex;

texture<float4, 1, cudaReadModeElementType> nodes2Tex;
texture<float4, 1, cudaReadModeElementType> nodeAttributes2Tex;
texture<int4,   1, cudaReadModeElementType> cells2Tex;
texture<int4,   1, cudaReadModeElementType> neighbors2Tex;

__device__ __constant__ int num_cells_m1;
__device__ __constant__ int num_cells_m2;

// number of hops per particle while the flow field evaluation
__device__ int num_hops1[MAX_PARTICLES];
__device__ int num_hops2[MAX_PARTICLES];

// maximum number of traversed cells in tetrahedral walk
__device__ __constant__ int tw_max_hops = 40;

static __device__ float
calculateDeterminant(float4 p1, float4 p2, float4 p3, float4 p4)
{
	return (p2.x-p1.x)*((p3.y-p1.y)*(p4.z-p1.z)-(p3.z-p1.z)*(p4.y-p1.y)) +
		(p3.x-p1.x)*((p1.y-p2.y)*(p4.z-p1.z)-(p1.z-p2.z)*(p4.y-p1.y)) +
		(p4.x-p1.x)*((p2.y-p1.y)*(p3.z-p1.z)-(p2.z-p1.z)*(p3.y-p1.y));
}

static __device__ float4
calculateNaturalCoordinates( float4 p, float4 p1, float4 p2, float4 p3, float4 p4 )
{
	float4 res;

	float det = calculateDeterminant(p1, p2, p3, p4);

	if( det > 0 + EPS)
	{
		float a11 = (p4.z-p1.z)*(p3.y-p4.y) - (p3.z-p4.z)*(p4.y-p1.y);
		float a21 = (p4.z-p1.z)*(p1.y-p2.y) - (p1.z-p2.z)*(p4.y-p1.y);
		float a31 = (p2.z-p3.z)*(p1.y-p2.y) - (p1.z-p2.z)*(p2.y-p3.y);

		float a12 = (p4.x-p1.x)*(p3.z-p4.z) - (p3.x-p4.x)*(p4.z-p1.z);
		float a22 = (p4.x-p1.x)*(p1.z-p2.z) - (p1.x-p2.x)*(p4.z-p1.z);
		float a32 = (p2.x-p3.x)*(p1.z-p2.z) - (p1.x-p2.x)*(p2.z-p3.z);

		float a13 = (p4.y-p1.y)*(p3.x-p4.x) - (p3.y-p4.y)*(p4.x-p1.x);
		float a23 = (p4.y-p1.y)*(p1.x-p2.x) - (p1.y-p2.y)*(p4.x-p1.x);
		float a33 = (p2.y-p3.y)*(p1.x-p2.x) - (p1.y-p2.y)*(p2.x-p3.x);

		res.x = (a11*(p.x-p1.x) + a12*(p.y-p1.y) + a13*(p.z-p1.z)) / det;  // 1->2
		res.y = (a21*(p.x-p1.x) + a22*(p.y-p1.y) + a23*(p.z-p1.z)) / det;  // 1->3
		res.z = (a31*(p.x-p1.x) + a32*(p.y-p1.y) + a33*(p.z-p1.z)) / det;  // 1->4
		res.w = 1-res.x-res.y-res.z;
	}
	return res;
};

static __device__ bool
naturalCoordinatesValid( float4 crds )
{
	bool res = true;
	res &= ( crds.x+EPS > 0 );
	res &= ( crds.y+EPS > 0 );
	res &= ( crds.z+EPS > 0 );
	res &= ( crds.w+EPS > 0 );
	return res;
};

static __device__ int
calculateNextCell( float4 nc, int4 ns, int c )
{
	int res = c;

  // Calculate worst violator
	if((nc.x + EPS < 0) || (nc.y + EPS < 0) || 
     (nc.z + EPS < 0) || (nc.w + EPS < 0) )
  {
    if ( nc.x<nc.y && nc.x<nc.z && nc.x<nc.w )  
	{
		res = ns.y; // next cell shares points 1,3,4
	}
    else if( nc.y<nc.z && nc.y<nc.w ) 
	{
		res = ns.z; // next cell shares points 1,2,4
	}
  	else if( nc.z<nc.w ) 
	{
		res = ns.w; // next cell shares points 1,2,3
	}
	  else
    {
	  	res = ns.x; // next cell shares points 2,3,4
  	}
  }
  
	return res;
};

static __device__ float4
interpolateVelocity1( int4 indices, float4 crds)
{
	float4 v = make_float4( 0.0, 0.0, 0.0, 0.0);

	float4 v1 = tex1Dfetch(nodeAttributes1Tex, indices.x );
	v = v1;

	float4 v2 = tex1Dfetch(nodeAttributes1Tex, indices.y );
	v += (v2-v1) * clamp(crds.x,0.0f,1.0f);

	float4 v3 = tex1Dfetch(nodeAttributes1Tex, indices.z );
	v += (v3-v1) * clamp(crds.y,0.0f,1.0f);

	float4 v4 = tex1Dfetch(nodeAttributes1Tex, indices.w );
	v += (v4-v1) * clamp(crds.z,0.0f,1.0f);

	return v;
}

static __device__ float4
InterpolateVelocityAt1
(
 float4	pos, 
 int*	startCells, 
 char*	tc,
 int	tid,
 bool*	outOfField = 0
 )
{
	// cell index of current cell
	int	current_cell = startCells[ tid];

	// iterations counter
	int	stepCount = 0;

	// flags
	bool oof = ( current_cell < 0 );
	if( outOfField ) oof |= *outOfField;

	bool crdsValid = false;

	// natural coodinates
	float4 crds;

	// points indices of current cell
	int4 indices;

	while( !( crdsValid || oof ))
	{
		// write index of current cell to traversed cells field
		if( current_cell > 0 && current_cell < num_cells_m1) tc[ current_cell ] = 1;

		// read cell node indices from mem
		indices = tex1Dfetch(cells1Tex, current_cell);

		// read physical coordinates from mem
		float4 p1 = tex1Dfetch( nodes1Tex, indices.x);
		float4 p2 = tex1Dfetch( nodes1Tex, indices.y);
		float4 p3 = tex1Dfetch( nodes1Tex, indices.z);
		float4 p4 = tex1Dfetch( nodes1Tex, indices.w);

		// calculate natural coordinates of search point in current cell
		crds = calculateNaturalCoordinates( pos, p1, p2, p3, p4);

		int4 ns = tex1Dfetch( neighbors1Tex, current_cell);
		int next_cell = calculateNextCell( crds, ns, current_cell);

		crdsValid = ( next_cell == current_cell || stepCount > tw_max_hops);

		if( !crdsValid )
		{
			current_cell = next_cell;
			oof = ( current_cell < 0 );
			stepCount++;
		}
	}

	// save step count
	num_hops1[tid] += stepCount;

	// update startcell index
	startCells[ tid] = current_cell;

	// write oof state
	if( outOfField ) *outOfField = oof;

	// write interpolated velocity to global mem
	float4 res = make_float4( 0.0, 0.0, 0.0, 0.0 );

	if( crdsValid )
	{
		res = interpolateVelocity1( indices, crds ); 
	}
	return res;
};

 static __device__ float4
interpolateVelocity2( int4 indices, float4 crds)
{
	float4 v = make_float4( 0.0, 0.0, 0.0, 0.0);

	float4 v1 = tex1Dfetch(nodeAttributes2Tex, indices.x );
	v = v1;

	float4 v2 = tex1Dfetch(nodeAttributes2Tex, indices.y );
	v += (v2-v1) * clamp(crds.x,0.0f,1.0f);

	float4 v3 = tex1Dfetch(nodeAttributes2Tex, indices.z );
	v += (v3-v1) * clamp(crds.y,0.0f,1.0f);

	float4 v4 = tex1Dfetch(nodeAttributes2Tex, indices.w );
	v += (v4-v1) * clamp(crds.z,0.0f,1.0f);

	return v;

}

 static __device__ float4
InterpolateVelocityAt2
(
 float4	pos, 
 int*	startCells, 
 char*	tc,
 int	tid,
 bool*	outOfField = 0
 )
{
	// cell index of current cell
	int	current_cell = startCells[ tid];

	// iterations counter
	int	stepCount = 0;

	// flags
	bool oof = ( current_cell < 0 );
	if( outOfField ) oof |= *outOfField;

	bool crdsValid = false;

	// natural coodinates
	float4 crds;

	// points indices of current cell
	int4 indices;

	while( !( crdsValid || oof ))
	{
		// write index of current cell to traversed cells field
		if( current_cell > 0 && current_cell < num_cells_m2 ) tc[ current_cell ] = 1;

		// read cell node indices from mem
		indices = tex1Dfetch(cells2Tex, current_cell);

		// read physical coordinates from mem
		float4 p1 = tex1Dfetch( nodes2Tex, indices.x);
		float4 p2 = tex1Dfetch( nodes2Tex, indices.y);
		float4 p3 = tex1Dfetch( nodes2Tex, indices.z);
		float4 p4 = tex1Dfetch( nodes2Tex, indices.w);

		// calculate natural coordinates of search point in current cell
		crds = calculateNaturalCoordinates( pos, p1, p2, p3, p4);

		int4 ns = tex1Dfetch( neighbors2Tex, current_cell);
		int next_cell = calculateNextCell( crds, ns, current_cell);

		crdsValid = ( next_cell == current_cell || stepCount > tw_max_hops );

		if( !crdsValid )
		{
			current_cell = next_cell;
			oof = ( current_cell < 0 );
			stepCount++;
		}
	}

	// save step count
	num_hops2[tid] += stepCount;

	// update startcell index
	startCells[ tid] = current_cell;

	// write oof state
	if( outOfField ) *outOfField = oof;

	// write interpolated velocity to global mem
	float4 res = make_float4( 0.0, 0.0, 0.0, 0.0 );

	if( crdsValid )
	{
		res = interpolateVelocity2( indices, crds ); 
	}
	return res;
};


#endif
