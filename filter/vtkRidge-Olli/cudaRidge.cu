// System includes
#include "cudaRidge.h"

#include "mctable.h"
#include "kernels.cu"

//#define CUDA_OUTPUT_STEPS

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

extern "C" bool cudaGradient(const double* d, uint dim[3], double** gradients) {
    // for (uint i=0; i<dim[0]*dim[1]*dim[2]; ++i)
    // 	cout << d[i] << endl;
#ifdef CUDA_OUTPUT_STEPS
    std::cout << "    CUDA: Starting cuda method..." << std::endl;
#endif

    if (cudaSetDevice(0) != cudaSuccess) {
	std::cout << "    CUDA: Could not set cuda device." << std::endl;
    }
    else {
#ifdef CUDA_OUTPUT_STEPS
	std::cout << "    CUDA: Cuda device set." << std::endl;
#endif
    }

    //timespec start, end;

    uint numPoints=dim[0]*dim[1]*dim[2];
    uint sg = 3*numPoints*sizeof(double);
    uint sd = numPoints*sizeof(double);
    uint sdim = 3*sizeof(uint);
    uint *dim_d;
    double *g_h = new double[3*numPoints];
    double *g_d, *d_d;

    //clock_gettime(CLOCK_MONOTONIC, &start);

    gpuErrchk(cudaMalloc((void**)&d_d, sd));
    gpuErrchk(cudaMalloc((void**)&g_d, sg));
    gpuErrchk(cudaMalloc((void**)&dim_d, sdim));
#ifdef CUDA_OUTPUT_STEPS
    std::cout << "    CUDA: Memory allocated." << std::endl;
#endif

    gpuErrchk(cudaMemcpy(d_d, d, sd, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(dim_d, dim, sdim, cudaMemcpyHostToDevice));
#ifdef CUDA_OUTPUT_STEPS
    std::cout << "    CUDA: Memory copied." << std::endl;
#endif
    
    uint blockSize = 128;
    int nBlocks = numPoints/blockSize + (numPoints%blockSize == 0?0:1);
    processPoints <<< nBlocks, blockSize >>> (numPoints, d_d, g_d, dim_d);

    gpuErrchk(cudaMemcpy(g_h, g_d, sg, cudaMemcpyDeviceToHost));

    //clock_gettime(CLOCK_MONOTONIC, &end);
    
    gpuErrchk(cudaFree(dim_d));
    gpuErrchk(cudaFree(d_d));
    gpuErrchk(cudaFree(g_d));

    *gradients = g_h;
    //cout << "Time used for GPU calculation in milliseconds:" << double(end.tv_nsec-start.tv_nsec)/1000000.0 << endl;
    return true;
}


extern "C" bool cudaIsosurface(const double c, const double* in_data, double* in_grads, int in_size[3], uint* out_numVerts, double** out_verts, double** out_grads)
{
#ifdef CUDA_OUTPUT_STEPS
    std::cout << "    CUDA: Starting cuda method..." << std::endl;
#endif

    if (cudaSetDevice(0) != cudaSuccess) {
	std::cout << "    CUDA: Could not set cuda device." << std::endl;
    }
    else {
#ifdef CUDA_OUTPUT_STEPS
	std::cout << "    CUDA: Cuda device set." << std::endl;
#endif
    }
    //---Init---------------------------------------------------------------------    
    //timespec start, end;
    uint sdataB, sgradsB, sdimB, sbitmasksB, eTB, tTB, nVTB;
    double *data, *g;
    int *dim;
    uint *bitmasks, bitmaskCnt, *eT, *tT, *nVT;
    bitmaskCnt = in_size[0]*in_size[1]*in_size[2];
    uint *bm = new uint[2*bitmaskCnt];
    sdataB = in_size[0]*in_size[1]*in_size[2]*sizeof(double);
    sgradsB = 3*in_size[0]*in_size[1]*in_size[2]*sizeof(double);
    sdimB = 3*sizeof(int);
    sbitmasksB = 2*bitmaskCnt*sizeof(uint);
    eTB = 256*sizeof(uint);
    tTB = 256*16*sizeof(uint);
    nVTB = 256*sizeof(uint);
    
    gpuErrchk(cudaMalloc(     (void**)&data,      sdataB));
    gpuErrchk(cudaMalloc(      (void**)&dim,       sdimB));
    gpuErrchk(cudaMalloc( (void**)&bitmasks,  sbitmasksB));
    gpuErrchk(cudaMalloc(       (void**)&eT,         eTB));
    gpuErrchk(cudaMalloc(       (void**)&tT,         tTB));
    gpuErrchk(cudaMalloc(      (void**)&nVT,        nVTB));
    gpuErrchk(cudaMalloc(        (void**)&g,     sgradsB));
#ifdef CUDA_OUTPUT_STEPS
    std::cout << "    CUDA: Memory allocated." << std::endl;
#endif

    gpuErrchk(cudaMemcpy(     data,       in_data,       sdataB, cudaMemcpyHostToDevice));    
    gpuErrchk(cudaMemcpy(      dim,       in_size,        sdimB, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy( bitmasks,            bm,   sbitmasksB, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(       eT,     edgeTable,          eTB, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(       tT,      triTable,          tTB, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(      nVT, numVertsTable,         nVTB, cudaMemcpyHostToDevice));
    gpuErrchk(cudaMemcpy(        g,      in_grads,      sgradsB, cudaMemcpyHostToDevice));
#ifdef CUDA_OUTPUT_STEPS
    std::cout << "    CUDA: Memory copied." << std::endl;
#endif

    //---Calculation--------------------------------------------------------------
    //clock_gettime(CLOCK_MONOTONIC, &start);

    int blockSize = 128;
    int nBlocks = bitmaskCnt/blockSize + (bitmaskCnt%blockSize == 0?0:1);
    processDataArray <<< nBlocks, blockSize >>> (bitmaskCnt, c, dim, data, bitmasks, nVT);

#ifdef CUDA_OUTPUT_STEPS
    std::cout << "    CUDA: Calculation done." << std::endl;
#endif

    gpuErrchk(cudaMemcpy(bm, bitmasks, sbitmasksB, cudaMemcpyDeviceToHost));
#ifdef CUDA_OUTPUT_STEPS
    std::cout << "    CUDA: Result retrieved." << std::endl;
#endif

    int addOne = 0;
    int addTris = 0;
    if (bm[2*bitmaskCnt-1]) {
	addOne = 1;                    //doing exclusive scan, if last element is set, we need one more space
	addTris = bm[2*bitmaskCnt-1];  //and the tris of course
    }

    uint *pfArrays = new uint[2*bitmaskCnt];

    thrust::exclusive_scan(&bm[bitmaskCnt], bm+2*bitmaskCnt, &pfArrays[bitmaskCnt]);  //index of next cubes
    thrust::exclusive_scan(&bm[0], bm+bitmaskCnt, &pfArrays[0]);                      

    uint numIsoSurfaceCubes=pfArrays[2*bitmaskCnt-1]+addOne;
    uint numVertices=pfArrays[bitmaskCnt-1]+addTris;

#ifdef CUDA_OUTPUT_STEPS
    cout << "    CUDA: numIsoSurfaceCubes: " << numIsoSurfaceCubes << endl;
    cout << "    CUDA: numVertices: " << numVertices << endl;
#endif

    size_t iAB, vAB, gAB, pfAB;
    uint *pfA, *iA;
    double *vA, *vertexArray = new double[numVertices*3];
    double *gA, *gradientArray = new double[numVertices*3];

    iAB = numIsoSurfaceCubes*sizeof(uint);
    vAB = numVertices*3*sizeof(double);
    gAB = numVertices*3*sizeof(double);
    pfAB = sbitmasksB;

    gpuErrchk(cudaMalloc( (void**)&iA,  iAB));
    gpuErrchk(cudaMalloc( (void**)&vA,  vAB));
    gpuErrchk(cudaMalloc( (void**)&gA,  gAB));
    gpuErrchk(cudaMalloc((void**)&pfA, pfAB));

    gpuErrchk(cudaMemcpy( pfA, pfArrays, pfAB, cudaMemcpyHostToDevice));   //copy prefix array for second pass
    
    getIsoIndices <<< nBlocks, blockSize >>> (bitmaskCnt, iA, bitmasks, pfA);

    nBlocks = numIsoSurfaceCubes/blockSize + (numIsoSurfaceCubes%blockSize == 0?0:1);
    processIsoCubes <<< nBlocks, blockSize >>> (numIsoSurfaceCubes, c, iA, pfA, data, g, vA, gA, nVT, tT, bitmasks, dim);

    gpuErrchk(cudaMemcpy(vertexArray, vA, vAB, cudaMemcpyDeviceToHost ));
    gpuErrchk(cudaMemcpy(gradientArray, gA, gAB, cudaMemcpyDeviceToHost ));

    //clock_gettime(CLOCK_MONOTONIC, &end);
    
    //---Cleanup------------------------------------------------------------------
    cudaFree(data);
    cudaFree(g);
    cudaFree(dim);
    cudaFree(bitmasks);
    cudaFree(eT);
    cudaFree(tT);
    cudaFree(nVT);
    cudaFree(pfA);
    cudaFree(iA);
    cudaFree(vA);
    cudaFree(gA);

    delete [] pfArrays;
    delete [] bm;

    *out_verts = vertexArray;
    *out_numVerts = numVertices;

#ifdef CUDA_OUTPUT_STEPS
    cout << "    CUDA: Cuda calculation done." << endl;
#endif
    //cout << "Time used for GPU calculation in milliseconds:" << double(end.tv_nsec-start.tv_nsec)/1000000.0 << endl;
    return true;
}
