
#ifndef PARTICLES_CU
#define PARTICLES_CU

// includes, project
#include <helper_cuda.h>

#include "Particles_kernel.cu"
#include "Integration_kernel.cu"


extern "C"
{
	void setATSParameters( ATSParams *hostParams)
	{
		// copy parameters to constant memory
		checkCudaErrors( cudaMemcpyToSymbol( params, hostParams, sizeof(ATSParams)) );
	}
	
	void setMaxParticleSteps( int num_steps )
	{
		checkCudaErrors( cudaMemcpyToSymbol( max_particle_steps, &num_steps, sizeof(int),0));
	}


	void AdvectParticlesCuda( float4* dPositions, float4* dVelocities, int numParticles, float timestep, uchar* occupiedBlocks, cudaStream_t stream)
    {
		if( numParticles < 1 ) return;

		// setup execution parameters
        int threadsPerBlock = 256; // 256 threads per block
        int numBlocks = (numParticles + threadsPerBlock-1) / threadsPerBlock;

        CREATE_TIMER
        START_TIMER(stream)

        //// execute the kernel
        AdvectParticlesKernel<<< numBlocks, threadsPerBlock, 0, stream >>> ( dPositions, dVelocities, timestep, numParticles, occupiedBlocks );
        getLastCudaError("Kernel execution failed");

        STOP_TIMER(stream)
        PRINT_TIMER("Particle Advection: ", "\n")
        STORE_TIMER(advectionTimer)
        DESTROY_TIMER
    }; 

	void IntegrateVelocitiesCuda
	(
		IntegrationScheme integration,
		float4* dPositions, float4* dVelocities, int numParticles,
		vistime_t* vistime,
		MeshGPU mesh1, MeshGPU mesh2,
		int* dStartCells1, int* dStartCells2, 
		float* dStepSizes, cudaStream_t stream
	)
    {

		if( numParticles < 1 ) return;

		checkCudaErrors(cudaBindTexture(0, nodes1Tex,		mesh1.nodes,		mesh1.s_nodes));
		checkCudaErrors(cudaBindTexture(0, cells1Tex,		mesh1.cells,		mesh1.s_cells));
		checkCudaErrors(cudaBindTexture(0, neighbors1Tex,	mesh1.neighbors,	mesh1.s_cells));
		checkCudaErrors(cudaBindTexture(0, nodeAttributes1Tex, mesh1.nodeAttributes, mesh1.s_nodes));

		checkCudaErrors(cudaBindTexture(0, nodes2Tex,		mesh2.nodes,		mesh2.s_nodes));
		checkCudaErrors(cudaBindTexture(0, cells2Tex,		mesh2.cells,		mesh2.s_cells));
		checkCudaErrors(cudaBindTexture(0, neighbors2Tex,	mesh2.neighbors,	mesh2.s_cells));
		checkCudaErrors(cudaBindTexture(0, nodeAttributes2Tex, mesh2.nodeAttributes, mesh2.s_nodes));

		// copy vistime to constant device memory
		checkCudaErrors( cudaMemcpyToSymbolAsync( vt,	vistime, sizeof(vistime_t), 0, cudaMemcpyHostToDevice, stream ));

		int* nc1;
		allocatePageLockedArray((void**)&nc1, sizeof(int), false);
		*(nc1) = mesh1.num_cells;
		checkCudaErrors( cudaMemcpyToSymbolAsync( num_cells_m1, nc1, sizeof(int), 0, cudaMemcpyHostToDevice, stream));

		int* nc2;
		allocatePageLockedArray((void**)&nc2, sizeof(int), false);
		*(nc2) = mesh2.num_cells;
		checkCudaErrors( cudaMemcpyToSymbolAsync( num_cells_m2, nc2, sizeof(int), 0, cudaMemcpyHostToDevice, stream));

		// reset number of hops for current iteration
		int* dSymbol;
		checkCudaErrors( cudaGetSymbolAddress((void**)&dSymbol, num_hops1));
		checkCudaErrors( cudaMemset( dSymbol, 0, numParticles * sizeof(int)));
		checkCudaErrors( cudaGetSymbolAddress((void**)&dSymbol, num_hops2));
		checkCudaErrors( cudaMemset( dSymbol, 0, numParticles * sizeof(int)));

		// reset iterations counter
		checkCudaErrors( cudaGetSymbolAddress((void**)&dSymbol, iterations));
		checkCudaErrors( cudaMemset( dSymbol, 0, numParticles * sizeof(int)));

		// setup execution parameters
		int threadsPerBlock = 128; // number of Particles per block
		int numBlocks = ( numParticles + threadsPerBlock-1) / threadsPerBlock;

        CREATE_TIMER
        START_TIMER(stream)

		// execute the kernel
		switch( integration )
		{
		default:
		case Euler:
			IntegrateVelocitiesEulerKernel<<< numBlocks, threadsPerBlock, 0, stream >>>
			(
				dPositions, dVelocities, numParticles,
				mesh1.traversedCells, mesh2.traversedCells,
				dStartCells1, dStartCells2
			);
			break;
		case RK3:
			IntegrateVelocitiesRK3Kernel<<< numBlocks, threadsPerBlock, 0, stream >>>
			(
				dPositions, dVelocities, numParticles,
				mesh1.traversedCells, mesh2.traversedCells,
				dStartCells1, dStartCells2
			);
			break;
		case RK4:
			IntegrateVelocitiesRK4Kernel<<< numBlocks, threadsPerBlock, 0, stream >>>
			(
				dPositions, dVelocities, numParticles,
				mesh1.traversedCells, mesh2.traversedCells,
				dStartCells1, dStartCells2
			);
			break;
		case Dopri5:
			IntegrateVelocitiesDopri5Kernel<<< numBlocks, threadsPerBlock, 0, stream >>>
			(
				dPositions, dVelocities, numParticles,
				mesh1.traversedCells, mesh2.traversedCells, 
				dStartCells1, dStartCells2
			);
			break;
		case Dopri5_ATS:
			IntegrateVelocitiesDopri5_ATS_Kernel<<< numBlocks, threadsPerBlock, 0, stream >>>
			(
				dPositions, dVelocities, numParticles,
				mesh1.traversedCells, mesh2.traversedCells,
				dStartCells1, dStartCells2,
				dStepSizes
			);
			break;
		}
   
		getLastCudaError("Kernel execution failed");

        STOP_TIMER(stream)
		PRINT_TIMER("Integration: ", "\n")
		STORE_TIMER(integrationTimer)
        DESTROY_TIMER


		checkCudaErrors(cudaUnbindTexture(nodes1Tex));
		checkCudaErrors(cudaUnbindTexture(nodeAttributes1Tex));
		checkCudaErrors(cudaUnbindTexture(cells1Tex));
		checkCudaErrors(cudaUnbindTexture(neighbors1Tex));

		checkCudaErrors(cudaUnbindTexture(nodes2Tex));
		checkCudaErrors(cudaUnbindTexture(nodeAttributes2Tex));
		checkCudaErrors(cudaUnbindTexture(cells2Tex));
		checkCudaErrors(cudaUnbindTexture(neighbors2Tex));

		freePageLockedHostMemory(nc1);
		freePageLockedHostMemory(nc2);
	};

	void getNumberOfHops( int count, int* hNumHops1, int* hNumHops2, int* hIterations )
	{
		cudaMemcpyFromSymbol(hNumHops1, num_hops1, count * sizeof(int), 0, cudaMemcpyDeviceToHost );
		cudaMemcpyFromSymbol(hNumHops2, num_hops2, count * sizeof(int), 0, cudaMemcpyDeviceToHost );
		cudaMemcpyFromSymbol(hIterations, iterations, count * sizeof(int), 0, cudaMemcpyDeviceToHost );
	}
}

#endif
