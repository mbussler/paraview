
#ifndef COMMON_CUH
#define COMMON_CUH

// includes
#include <GL/glew.h>
#include <GL/freeglut.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_math.h>
#include <helper_functions.h>

#include <cuda_gl_interop.h>
#include <vector_types.h>

#define MAX_PARTICLES 2000000

#define USE_TEX 1

#ifdef KERNEL_VERBOSE
	#define CREATE_TIMER cudaEvent_t start, stop; cudaEventCreate(&start); cudaEventCreate(&stop);
	#define START_TIMER(s) cudaEventRecord(start, s);
	#define STOP_TIMER(s) cudaEventRecord(stop, s); cudaEventSynchronize(stop);
	#define PRINT_TIMER(s,t) float elapsedTime; cudaEventElapsedTime(&elapsedTime, start, stop); printf(s); printf("%.4f ms", elapsedTime); printf(t);
	#define STORE_TIMER(t) // float elapsedTime2; cudaEventElapsedTime(&elapsedTime2, start, stop); t.push_back(elapsedTime2);
	#define DESTROY_TIMER cudaEventDestroy(start); cudaEventDestroy(stop);
#else
	#define CREATE_TIMER
	#define START_TIMER(s)
	#define STOP_TIMER(s)
	#define PRINT_TIMER(s,t)
	#define STORE_TIMER(t)
	#define DESTROY_TIMER
#endif

// storage for timer values
//static vector<float> advectionTimer;
//static vector<float> integrationTimer;
//static vector<float> tetwalkTimer;
//static vector<float> kdTreeTimer;
//static vector<float> dataCopyTimer;

// global timer event
static cudaEvent_t gen_start, gen_stop;
static float gen_elapsed;

#if USE_TEX
    #define FETCH(t, i) tex1Dfetch(t##Tex, i)
#else
    #define FETCH(t, i) t[i]
#endif

#define EPS 0.00001f
#define clamp(x,lo,hi) (x<lo?lo:(x>hi?hi:x))
#define clamp_lo(x,lo) (x<lo?lo:x)
#define clamp_hi(x,hi) (x>hi?hi:x)


// Integration Schemes
typedef enum {
	// Euler integration Scheme
	Euler = 0,
	// Runge-Kutta 3rd-order integration Scheme
	RK3,
	// Runge-Kutta 4rd-order integration Scheme
	RK4,
	// Dormand-Prince 4th- and 5th-order integration Scheme
	Dopri5,
	// DoPri-5 with adaptive timestepping
	Dopri5_ATS
} IntegrationScheme;

typedef unsigned int uint;
typedef unsigned char uchar;

typedef struct {
    float4* nodes;
    float4* nodeAttributes;
    int4*   cells;
    int4*   neighbors;
	char*	traversedCells;
    // size of mesh data on device
    size_t  s_nodes;
    size_t  s_cells;

    int		num_cells;
} MeshGPU;

typedef struct {

	// device pointers
    float*  dS;	// Split values
    char*   dD; // Split dimensions
	int*	dI; // Point indices
    int*    dL; // Cell indices

	// mem sizes
    size_t sS;
    size_t sD;
	size_t sI;
    size_t sL;

	// number of levels of kd tree
	int levels;

} KdTreeGPU;

 typedef struct {
	float t;
	float h;
	float a;
	float ah;
} vistime_t;

// Parameters for adaptive time stepping with dopri-5
struct ATSParams {
	// maximum tolerance for error
	float tol; // = 10e-3;
	// minimum step width 
	float h_min; // = 0.005f;
	// Safetyfactor \rho \in (0, 1];
	float rho; // = 0.9f;
	// magnification barrier \nu >= 1
	float eta; // = 2;
};

extern "C"
{
    void cudaInit(int argc, char **argv);
    void cudaGLInit(int argc, char **argv);

    // Synchronous Memory Allocation and Transfer functions
    void allocateArray(void **devPtr, size_t size);
    void freeArray(void *devPtr);
    void setArray(void* devPtr, int value, size_t count);

    void copyArrayToDevice(void* device, const void* host, int size);
	void copyArrayFromDevice(void* host, const void* device, int size, struct cudaGraphicsResource **cuda_vbo_resource = 0);
    void copyToConstantMem( const char* symbol, const void* src, size_t count );

    // Asynchronous Memory Allocation and Transfer functions
    void allocatePageLockedArray(void **hostPtr, size_t size, bool wc = false);
    void allocatePageLockedArrayPortable(void **hostPtr, size_t size, bool wc);
    void freePageLockedHostMemory(void *hostPtr);
    void createStreams( int numStreams, cudaStream_t* streams );

    void copyArrayToPageLockedHostMemory(void *hostPtr, void* src, size_t size);
	void copyArrayToDeviceAsync(void* device, const void* host, int size, cudaStream_t stream);
	void copyArrayFromDeviceAsync(void* host, const void* device, int size, cudaStream_t stream);

    void copyToSymbolAsync( const char *symbol, const void *src, size_t count, cudaStream_t stream = 0);
    void resetSymbol( const char *symbol, size_t count );

	// Synchronize all current threads
	void threadSync();

    void registerGLBufferObject(uint vbo, struct cudaGraphicsResource **cuda_vbo_resource);
    void unregisterGLBufferObject(struct cudaGraphicsResource *cuda_vbo_resource);
    void *mapGLBufferObject(struct cudaGraphicsResource **cuda_vbo_resource);
    void unmapGLBufferObject(struct cudaGraphicsResource *cuda_vbo_resource);

	void getGPUMemoryUsage(size_t* free, size_t* total, int divBy = 1);
	//void getTimerAveraged(vector<float>* timer, float* avg_time );
	//void getTimerMedian(vector<float>* timer, float* med_time );

	// General purpose timer functions
	void startTimer(cudaStream_t s);
	void stopTimer(cudaStream_t s);
	void printTimer();
	void destroyTimer(cudaStream_t s);

	int cudaDeviceInit(int ARGC, char **ARGV);
	int getMaxGflopsDeviceId();
}
#endif
