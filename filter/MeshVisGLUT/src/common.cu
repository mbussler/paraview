
#ifndef COMMON_CU
#define COMMON_CU

// includes
#include "common.cuh"

extern "C"
{
    void cudaInit(int argc, char **argv)
    {   
        // use command-line specified CUDA device, otherwise use device with highest Gflops/s
        if( checkCmdLineFlag(argc, (const char**)argv, "device") ) {
            cudaDeviceInit(argc, argv);
        } else {
            cudaSetDevice( getMaxGflopsDeviceId() );
        }
    };

    void cudaGLInit(int argc, char **argv)
    {   
        // use command-line specified CUDA device, otherwise use device with highest Gflops/s
        if( checkCmdLineFlag(argc, (const char**)argv, "device") ) {
            cudaDeviceInit(argc, argv);
        } else {
            cudaGLSetGLDevice( getMaxGflopsDeviceId() );
        }
    };

    void allocateArray(void **devPtr, size_t size)
    {
		cudaError err = cudaMalloc(devPtr, size);
	    if( cudaSuccess != err) {
			printf((stderr, "%s(%i) : cudaSafeCall() Runtime API error : %s.\n",
					__FILE__, __LINE__, cudaGetErrorString( err) ));
			exit(-1);
	    }
    };


    void freeArray(void *devPtr)
    {
		cudaError err = cudaFree(devPtr);
	    if( cudaSuccess != err) {
			printf((stderr, "%s(%i) : cudaSafeCall() Runtime API error : %s.\n",
					__FILE__, __LINE__, cudaGetErrorString( err) ));
			exit(-1);
	    }
    };

    void allocatePageLockedArray(void **hostPtr, size_t size, bool wc)
    {
#if CUDART_VERSION >= 2020
		checkCudaErrors( cudaHostAlloc( hostPtr, size, (wc) ? cudaHostAllocWriteCombined : 0));
#else
		cutilSafeCall( cudaHostAlloc( hostPtr, size );
#endif
    };

    void allocatePageLockedArrayPortable(void **hostPtr, size_t size, bool wc)
    {
		//cutilSafeCall( cudaHostAlloc( hostPtr, size, cudaHostAllocPortable | (wc) ? cudaHostAllocWriteCombined : 0));
		checkCudaErrors( cudaHostAlloc( hostPtr, size, (wc) ? cudaHostAllocWriteCombined : 0));
    };

    void freePageLockedHostMemory(void *hostPtr)
    {
		checkCudaErrors( cudaFreeHost( hostPtr));
    };

    void createStreams( int numStreams, cudaStream_t* streams )
    {
    	for (int i = 0; i < numStreams; ++i)
    	{
    	    cudaStreamCreate(&streams[i]);
    	}
    }


    void copyArrayToPageLockedHostMemory(void *hostPtr, void* src, size_t size)
    {
		checkCudaErrors( cudaMemcpy( hostPtr, src, size, cudaMemcpyHostToHost ));
    };

    void setArray(void* devPtr, int value, size_t count)
    {
		checkCudaErrors( cudaMemset( devPtr, value, count));
    };

	void copyArrayToDevice(void* device, const void* host, int size)
    {
        checkCudaErrors(cudaMemcpy(device, host, size, cudaMemcpyHostToDevice) );
    };

	void copyArrayToDeviceAsync(void* device, const void* host, int size, cudaStream_t stream)
    {
        checkCudaErrors(cudaMemcpyAsync(device, host, size, cudaMemcpyHostToDevice, stream));
    };

    void copyArrayFromDevice(void* host, const void* device, int size, struct cudaGraphicsResource **cuda_vbo_resource)
    {   

	#ifndef NOGL
		if (cuda_vbo_resource)
			device = mapGLBufferObject(cuda_vbo_resource);
	#endif

        checkCudaErrors(cudaMemcpy(host, device, size, cudaMemcpyDeviceToHost));
        
	#ifndef NOGL
		if (cuda_vbo_resource)
			unmapGLBufferObject(*cuda_vbo_resource);
	#endif
    };

    void copyArrayFromDeviceAsync(void* host, const void* device, int size, cudaStream_t stream)
    {
        checkCudaErrors(cudaMemcpyAsync(host, device, size, cudaMemcpyDeviceToHost, stream));
    };


    void copyToConstantMem( const char* symbol, const void* src, size_t count )
    {   
        checkCudaErrors(cudaMemcpyToSymbol( symbol, src, count, 0, cudaMemcpyHostToDevice ));
    };

    void copyToSymbolAsync( const char *symbol, const void *src, size_t count, cudaStream_t stream)
    {
    	checkCudaErrors( cudaMemcpyToSymbolAsync( symbol, src, count, 0, cudaMemcpyHostToDevice, stream ));
    }

    void resetSymbol( const char *symbol, size_t count )
    {
		int* dSymbol;
    	checkCudaErrors( cudaGetSymbolAddress((void**)&dSymbol, symbol ));
		checkCudaErrors( cudaMemset( dSymbol, 0, count ));
    }

    void threadSync()
	{
        checkCudaErrors(cudaThreadSynchronize());
    };

    void registerGLBufferObject(uint vbo, struct cudaGraphicsResource **cuda_vbo_resource)
    {
	#ifndef NOGL
		checkCudaErrors(cudaGraphicsGLRegisterBuffer(cuda_vbo_resource, vbo, cudaGraphicsMapFlagsNone));
	#endif
    };

    void unregisterGLBufferObject(struct cudaGraphicsResource *cuda_vbo_resource)
    {
	#ifndef NOGL
        checkCudaErrors(cudaGraphicsUnregisterResource(cuda_vbo_resource));
	#endif
    };

    void *mapGLBufferObject(struct cudaGraphicsResource **cuda_vbo_resource)
    {
        void *ptr = 0;
        
	#ifndef NOGL
		checkCudaErrors(cudaGraphicsMapResources(1, cuda_vbo_resource, 0));
        size_t num_bytes; 
        checkCudaErrors(cudaGraphicsResourceGetMappedPointer((void **)&ptr, &num_bytes, *cuda_vbo_resource));
	#endif
    
		return ptr;
    };

    void unmapGLBufferObject(struct cudaGraphicsResource *cuda_vbo_resource)
    {
	#ifndef NOGL
       checkCudaErrors(cudaGraphicsUnmapResources(1, &cuda_vbo_resource, 0));
	#endif
    };

	void getGPUMemoryUsage(size_t* free, size_t* total, int divBy)
	{
		checkCudaErrors( cudaMemGetInfo( free, total ));
		*free /= divBy;
		*total /= divBy;
	};
 
	//void getTimerMedian(vector<float>* timer, float* med_time )
	//{
	//	if( !timer || !med_time || timer->empty()) return;

	//	vector<float>::iterator midpoint;
	//	midpoint = timer->begin() + (timer->end() - timer->begin())/2;
	//	nth_element(timer->begin(), midpoint, timer->end());
	//	*med_time = *midpoint;
	//};

	//void getTimerAveraged(vector<float>* timer, float* avg_time )
	//{
	//	if( !timer || !avg_time || timer->empty()) return;

	//    float total = 0.0;
	//    for (int i=0; i< timer->size(); i++) {
	//        total += timer->at(i);
	//    }
	//    total /= timer->size();
	//    *avg_time = total;
	//};

	void startTimer(cudaStream_t s)
	{
		cudaEventCreate(&gen_start); cudaEventCreate(&gen_stop);
		cudaEventRecord(gen_start, s);
	};
	void stopTimer(cudaStream_t s)
	{
		cudaEventRecord(gen_stop, s); cudaEventSynchronize(gen_stop);
	};
	void printTimer()
	{
		cudaEventElapsedTime(&gen_elapsed, gen_start, gen_stop);
		printf("%.4f ms", gen_elapsed);
	};
	void destroyTimer(cudaStream_t s)
	{
		cudaEventDestroy(gen_start);
		cudaEventDestroy(gen_stop);
	};


	//#define STORE_TIMER(t)   float elapsedTime2; cudaEventElapsedTime(&elapsedTime2, start, stop); t.push_back(elapsedTime2);

	/////////// functions added for CUDA 5.0 compability ///////////////

	inline int cudaDeviceInit(int ARGC, char **ARGV)
	{
		int cuDevice = 0;
		int deviceCount = 0;
		CUresult err = cuInit(0);

		if (CUDA_SUCCESS == err)
		{
			checkCudaErrors(cuDeviceGetCount(&deviceCount));
		}

		if (deviceCount == 0)
		{
			fprintf(stderr, "cudaDeviceInit error: no devices supporting CUDA\n");
			exit(EXIT_FAILURE);
		}

		int dev = 0;
		dev = getCmdLineArgumentInt(ARGC, (const char **) ARGV, "device=");

		if (dev < 0)
		{
			dev = 0;
		}

		if (dev > deviceCount-1)
		{
			fprintf(stderr, "\n");
			fprintf(stderr, ">> %d CUDA capable GPU device(s) detected. <<\n", deviceCount);
			fprintf(stderr, ">> cudaDeviceInit (-device=%d) is not a valid GPU device. <<\n", dev);
			fprintf(stderr, "\n");
			return -dev;
		}

		checkCudaErrors(cuDeviceGet(&cuDevice, dev));
		char name[100];
		cuDeviceGetName(name, 100, cuDevice);

		if (checkCmdLineFlag(ARGC, (const char **) ARGV, "quiet") == false)
		{
			printf("> Using CUDA Device [%d]: %s\n", dev, name);
		}

		return dev;
	}

		// This function returns the best GPU based on performance
	inline int getMaxGflopsDeviceId()
	{
		CUdevice current_device = 0, max_perf_device = 0;
		int device_count     = 0, sm_per_multiproc = 0;
		int max_compute_perf = 0, best_SM_arch     = 0;
		int major = 0, minor = 0, multiProcessorCount, clockRate;

		cuInit(0);
		checkCudaErrors(cuDeviceGetCount(&device_count));

		// Find the best major SM Architecture GPU device
		while (current_device < device_count)
		{
			checkCudaErrors(cuDeviceComputeCapability(&major, &minor, current_device));

			if (major > 0 && major < 9999)
			{
				best_SM_arch = MAX(best_SM_arch, major);
			}

			current_device++;
		}

		// Find the best CUDA capable GPU device
		current_device = 0;

		while (current_device < device_count)
		{
			checkCudaErrors(cuDeviceGetAttribute(&multiProcessorCount,
												 CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT,
												 current_device));
			checkCudaErrors(cuDeviceGetAttribute(&clockRate,
												 CU_DEVICE_ATTRIBUTE_CLOCK_RATE,
												 current_device));
			checkCudaErrors(cuDeviceComputeCapability(&major, &minor, current_device));

			if (major == 9999 && minor == 9999)
			{
				sm_per_multiproc = 1;
			}
			else
			{
				sm_per_multiproc = _ConvertSMVer2Cores(major, minor);
			}

			int compute_perf  = multiProcessorCount * sm_per_multiproc * clockRate;

			if (compute_perf  > max_compute_perf)
			{
				// If we find GPU with SM major > 2, search only these
				if (best_SM_arch > 2)
				{
					// If our device==dest_SM_arch, choose this, or else pass
					if (major == best_SM_arch)
					{
						max_compute_perf  = compute_perf;
						max_perf_device   = current_device;
					}
				}
				else
				{
					max_compute_perf  = compute_perf;
					max_perf_device   = current_device;
				}
			}

			++current_device;
		}

		return max_perf_device;
	}


} // extern "C"

#endif
