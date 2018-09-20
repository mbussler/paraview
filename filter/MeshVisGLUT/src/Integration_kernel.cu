

#ifndef INTEGRATION_KERNEL_CU
#define INTEGRATION_KERNEL_CU

// USE_TEX
#include "common.cuh"

// Interpolate velocity in tetrahedral grid
#include "MeshFunctions.cu"

// visualization time parameters
__device__ __constant__ vistime_t vt;

__device__ int iterations[MAX_PARTICLES];

// Parameters for the adaptive time stepping
__device__ __constant__ ATSParams params;

__global__ void
IntegrateVelocitiesEulerKernel
(
	float4* positions,
	float4* velocities,
	int numParticles,
	char* tc1, char* tc2,
	int* sc1, int* sc2
)
{
	const unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if( tid < numParticles )
	{
		float4 pos = positions[tid];

		if( pos.w > 0 ) // test, if position is inside field
		{
			iterations[tid]++;

			bool outOfField = false;

			// flow field data evaluation
			float4 v1 = InterpolateVelocityAt1( pos, sc1, tc1, tid, &outOfField );
			float4 v2 = InterpolateVelocityAt2( pos, sc2, tc2, tid, &outOfField );

			// temporal interpolation
			velocities[tid] = ((1-vt.a)*v1) + (vt.a*v2);

			// test whether particle has exited the field
			if( outOfField ) 
			{
				positions[tid].w = -1.0f;
			}
		}
	}
};

__global__ void
IntegrateVelocitiesRK3Kernel
(
	float4* positions,
	float4* velocities,
	int numParticles,
	char* tc1, char* tc2,
	int* sc1, int* sc2
)
{

	const unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if( tid < numParticles )
	{
		float4 pos = positions[tid];

		if( pos.w > 0 ) // test, if position is oof
		{
			bool outOfField = false;

			iterations[tid]++;
			float4 vel = make_float4(0.0f,0.0f,0.0f,0.0f);

			float4 k1, kx, kx1, yx;

			kx  = InterpolateVelocityAt1( pos, sc1, tc1, tid, &outOfField );
			kx1 = InterpolateVelocityAt2( pos, sc2, tc2, tid, &outOfField );

			// interpolate between timesteps
			k1 = (1-vt.a)*kx + vt.a*kx1;
			vel += k1/6.0f;

			// calculate next postition to evaluate
			yx = pos + vt.h*(1.0f/2.0f*k1);

			kx  = InterpolateVelocityAt1( yx, sc1, tc1, tid, &outOfField );
			kx1 = InterpolateVelocityAt2( yx, sc2, tc2, tid, &outOfField );

			// for the temporal interpolation, the timestep is increased relative
			kx = (1-(vt.a+vt.ah/2.0f))*kx + (vt.a+vt.ah/2.0f)*kx1;
			vel += 4.0f/6.0f*kx;

			yx = pos + vt.h*(-1.0f*k1 + 2.0f*kx);

			kx  = InterpolateVelocityAt1( yx, sc1, tc1, tid, &outOfField );
			kx1 = InterpolateVelocityAt2( yx, sc2, tc2, tid, &outOfField );

			kx = (1-(vt.a+vt.ah))*kx + (vt.a+vt.ah)*kx1;
			vel += kx/3.0f;

			velocities[tid] = vel;

			if( outOfField ) 
			{
				positions[tid].w = -1.0f;
			}
		}
	}
};

__global__ void
IntegrateVelocitiesRK4Kernel
(
	float4* positions,
	float4* velocities,
	int numParticles,
	char* tc1, char* tc2,
	int* sc1, int* sc2
)
{

	const unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if( tid < numParticles )
	{
		float4 pos = positions[tid];

		if( pos.w > 0 ) // test, if position is oof
		{
			bool outOfField = false;
			iterations[tid]++;

			float4 vel = make_float4(0.0f,0.0f,0.0f,0.0f);
			float4 yx, kx, kx1;
			
			// velocity field interpolation for mesh 1 and mesh 2
			kx  = InterpolateVelocityAt1( pos, sc1, tc1, tid, &outOfField );
			kx1 = InterpolateVelocityAt2( pos, sc2, tc2, tid, &outOfField );

			// temporal interpolation
			kx = (1-vt.a)*kx + vt.a*kx1;
			vel += kx/6.0f;

			yx = pos + (vt.h/2.0f*kx);
			kx  = InterpolateVelocityAt1( yx, sc1, tc1, tid, &outOfField );
			kx1 = InterpolateVelocityAt2( yx, sc2, tc2, tid, &outOfField );

			kx = (1-(vt.a+vt.ah/2.0f))*kx + (vt.a+vt.ah/2.0f)*kx1;
			vel += kx/3.0f;

			yx = pos + (vt.h/2.0f*kx);
			kx  = InterpolateVelocityAt1( yx, sc1, tc1, tid, &outOfField );
			kx1 = InterpolateVelocityAt2( yx, sc2, tc2, tid, &outOfField );

			kx = (1-(vt.a+vt.ah/2.0f))*kx + (vt.a+vt.ah/2.0f)*kx1;
			vel += kx/3.0f;

			yx = pos + vt.h*kx;
			kx  = InterpolateVelocityAt1( yx, sc1, tc1, tid, &outOfField );
			kx1 = InterpolateVelocityAt2( yx, sc2, tc2, tid, &outOfField );
			
			kx = (1-(vt.a+vt.ah))*kx + (vt.a+vt.ah)*kx1;
			vel += kx/6.0f;

			velocities[tid] = vel;

			if( outOfField ) 
			{
				positions[tid].w = -1.0f;
			}
		}
	}
};

__device__ void IntegrateVelocitiesDopri5 ( float4 pos, int tid, 
											float h_curr, float a_curr, float ah_curr,
											char* tc1, char* tc2, int* sc1, int* sc2, 
											float4* velocity_rk4, float4* velocity_rk5, bool* outOfField );

__global__ void
IntegrateVelocitiesDopri5Kernel
(
	float4* positions,
	float4* velocities,
	int numParticles,
	char* tc1, char* tc2,
	int* sc1, int* sc2
)
{
	const unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;

	if( tid < numParticles )
	{
		float4 pos = positions[tid];

		// test, if position is inside flow field
		if( pos.w > 0 )
		{
			bool outOfField = false;
			float4 velocity_rk4, velocity_rk5;

			iterations[tid]++;

			IntegrateVelocitiesDopri5(	pos, tid,
										vt.h, vt.a, vt.ah,
										tc1, tc2, sc1, sc2,
										&velocity_rk4, &velocity_rk5, &outOfField );
			

			velocities[tid] = velocity_rk5;

			if( outOfField ) 
			{
				positions[tid].w = -1.0f;
			}
		}
	}
};

__global__ void
IntegrateVelocitiesDopri5_ATS_Kernel
(
	float4* positions,
	float4* velocities,
	int numParticles,
	char* tc1, char* tc2,
	int* sc1, int* sc2,
	float* stepSizes
)
{

	const unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;

	float T = vt.t + vt.h;  // maximum time for the current advection step
	float t_curr = vt.t; // current time during the adaptive timestep estimation
	float h_curr = vt.h; // the current stepsize, initialized with the global step size

	// the conversion factor between relative and real timestep
	__shared__ float f;
	f = vt.ah / vt.h;

	float a_curr = vt.a; // the current relative time

	if( tid < numParticles )
	{
		// set current stepsize to last one, if inside bounds
		float h_last = stepSizes[tid];

		if(( h_last>0 ) && ( h_last<h_curr )) { 
			h_curr = h_last; 
		} else {
			h_last = h_curr;
		}

		// read starting position from 
		float4 pos = positions[tid];

		float4 vel = make_float4( 0.0, 0.0, 0.0, 0.0);

		// test, if position is inside flow field
		bool outOfField = (pos.w < 0);

		while( !outOfField && h_curr>0 )
		{
	    iterations[tid]++;

			float4 velocity_rk4, velocity_rk5;

			// Calculate 4-th and 5-th order integration with Dopri-5 integration scheme
			IntegrateVelocitiesDopri5(	pos, tid,
										h_curr, a_curr, (f*h_curr),
										tc1, tc2, sc1, sc2,
										&velocity_rk4, &velocity_rk5, &outOfField );

			// calculate of error
			float error = length( velocity_rk5 - velocity_rk4 );

			// Accept current stepsize if error is 
			// below threshold or minimum stepsize reached.
			if( ((error / h_curr) <= params.tol ) || ( h_curr <= params.h_min ))
			{
				// set new position to evaluate
				pos += h_curr * velocity_rk5;
				vel += h_curr * velocity_rk5;

				// increment time and relative time
				t_curr += h_curr;
				a_curr += (f*h_curr);

				// estimate new stepsize
				h_curr = max( params.h_min, min( params.eta * h_curr, params.rho* pow((params.tol/error * pow(h_curr,5)),0.25f)));

				if(( t_curr+h_curr ) > T ) { 
					h_curr = T-t_curr; // The remaining timestep, typically much smaller than h_curr
				} else { 
					h_last = h_curr;
				}
			}
			// Otherwise the stepsize is halved
			else 
			{ 
				h_curr = h_curr/2.0f; 
			}
		} 

		// store last stepsize
		stepSizes[tid] = h_last;
		velocities[tid] = vel / vt.h;

		if( outOfField ) 
		{
			positions[tid].w = -1.0f;
		}
	}
};


__device__ void
IntegrateVelocitiesDopri5
(
	float4 pos,
	int tid,
	float h_curr,
	float a_curr,
	float ah_curr,
	char* tc1, char* tc2,
	int* sc1, int* sc2,
	float4* velocity_rk4,
	float4* velocity_rk5,
	bool* outOfField
)
{

	float ax;
	float4 y0, y1, y2, y3, y4, y5, y6;
	float4 k0, k1, k2, k3, k4, k5, k6;
	float4 kx1, kx2;

	// flow field evaluation in mesh 1 and mesh 2
	y0 = pos;
	kx1 = InterpolateVelocityAt1( y0, sc1, tc1, tid, outOfField );
	kx2 = InterpolateVelocityAt2( y0, sc2, tc2, tid, outOfField );

	// temporal interpolation
	ax = a_curr;

	k0 = (1-ax)*kx1 + ax*kx2;
	y1 = y0 + h_curr * (k0/5.0f);

	kx1 = InterpolateVelocityAt1( y1, sc1, tc1, tid, outOfField );
	kx2 = InterpolateVelocityAt2( y1, sc2, tc2, tid, outOfField );
	ax = a_curr + (ah_curr/5.0f);

	k1 = (1-ax)*kx1 + ax*kx2;
	y2 = y0 + h_curr*(	3.0f/40.0f*k0 + 
						9.0f/40.0f*k1 );

	kx1 = InterpolateVelocityAt1( y2, sc1, tc1, tid, outOfField );
	kx2 = InterpolateVelocityAt2( y2, sc2, tc2, tid, outOfField );
	ax = a_curr + ah_curr*(3.0f/10.0f );

	k2 = (1-ax)*kx1 + ax*kx2;
	y3 = y0 + h_curr*(	44.0f/45.0f * k0 - 
						56.0f/15.0f * k1 + 
						32.0f/9.0f  * k2 );

	kx1 = InterpolateVelocityAt1( y3, sc1, tc1, tid, outOfField );
	kx2 = InterpolateVelocityAt2( y3, sc2, tc2, tid, outOfField );
	ax = a_curr + ah_curr*(4.0f/5.0f);

	k3 = (1-ax)*kx1 + ax*kx2;
	y4 = y0 + h_curr*(	19372.0f/6561.0f  * k0 -
						25360.0f/2187.0f  * k1 + 
						64448.0f/6561.0f  * k2 -
						212.0f/729.0f     * k3 );
	kx1 = InterpolateVelocityAt1( y4, sc1, tc1, tid, outOfField );
	kx2 = InterpolateVelocityAt2( y4, sc2, tc2, tid, outOfField );
	ax = a_curr + ah_curr*(8.0f/ 9.0f);

	k4 = (1-ax)*kx1 + ax*kx2;
	y5 = y0 + h_curr*(	9017.0f/3168.0f   * k0 -
						355.0f/33.0f      * k1 + 
						46732.0f/5247.0f  * k2 +
						49.0f/176.0f      * k3 -
						5103.0f/18656.0f  * k4);

	kx1 = InterpolateVelocityAt1( y5, sc1, tc1, tid, outOfField );
	kx2 = InterpolateVelocityAt2( y5, sc2, tc2, tid, outOfField );
	ax = a_curr + ah_curr;

	k5 = (1-ax)*kx1 + ax*kx2;
	y6 = y0 + h_curr*(	35.0f/384.0f     * k0 +
						//0.0f           * k1 + 
						500.0f/1113.0f   * k2 +
						125.0f/192.0f    * k3 -
						2187.0f/6784.0f  * k4 +
						11.0f/84.0f      * k5 );

	kx1 = InterpolateVelocityAt1( y6, sc1, tc1, tid, outOfField );
	kx2 = InterpolateVelocityAt2( y6, sc2, tc2, tid, outOfField );
	ax = a_curr + ah_curr;

	k6 = (1-ax)*kx1 + ax*kx2;
	*velocity_rk4 = \
		5179.0f / 57600.0f 	* k0 +
		// 0.0f 			* k1 +
		7571.0f / 16695.0f 	* k2 +
		393.0f / 640.0f 	* k3 -
		92097.0f / 339200.0f * k4 +
		187.0f / 2100.0f 	* k5 +
		1.0f / 40.0f 		* k6 ;

	*velocity_rk5 = \
		35.0f / 384.0f 		* k0 + 
		//0.0f 				* k1 +
		500.0f / 1113.0f 	* k2 +
		125.0f / 192.0f 	* k3 -
		2187.0f / 6784.0f 	* k4 +
		11.0f / 84.0f 		* k5 ;

};


#endif
