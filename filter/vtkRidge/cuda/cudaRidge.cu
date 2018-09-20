// System includes
#include "cudaRidge.h"
#include <vector>
#include <set>
#include <thrust/device_ptr.h>

#include "../mctable.h"
#include "kernels.cu"

#define CUDA_OUTPUT_STEPS
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

typedef double   vec3[3];
typedef  vec3    mat3[3];

bool cudaRidgeData(uint dim[3], const double *data, double** evals_out, uint* numVertices_out, double** vertices_out, double** pointData_out, double featureThreshold, double origin[3], double spacing[3], unsigned int stencilRange, bool valley) {
#ifdef CUDA_OUTPUT_STEPS
    std::cout << "    CUDA: Starting cuda method..." << std::endl;
#endif

    int deviceCount = 0;
    gpuErrchk(cudaGetDeviceCount(&deviceCount));
    std::cout << "    CUDA: number of devices is " << deviceCount << std::endl;

    if (cudaSetDevice(0) != cudaSuccess) {
    std::cout << "    CUDA: Could not set cuda device." << std::endl;
    }
    else {
#ifdef CUDA_OUTPUT_STEPS
    std::cout << "    CUDA: Cuda device set." << std::endl;
#endif
    }
    size_t numPointsTotal=dim[0]*dim[1]*dim[2];
    size_t sgTotal = 3*numPointsTotal*sizeof(double);
    size_t shesseTotal = 3*3*numPointsTotal*sizeof(double);
    size_t seValTotal = 3*numPointsTotal*sizeof(double);
    size_t seVectTotal = 3*numPointsTotal*sizeof(double);
    size_t slargeEVTotal = 3*numPointsTotal*sizeof(double);
    size_t sdTotal = numPointsTotal*sizeof(double);
    size_t sbitmasksTotal = 2*numPointsTotal*sizeof(uint);
    size_t sdim = 3*sizeof(uint);
    size_t snVT = 256*sizeof(uint);
    size_t stT = 256*16*sizeof(uint);
    uint *dim_d;

    size_t memoryRequired = sdTotal + sgTotal + seValTotal + seVectTotal + slargeEVTotal + shesseTotal + sbitmasksTotal + sdim + snVT + stT;
    size_t memoryFree, memoryTotal;
    cudaMemGetInfo(&memoryFree,&memoryTotal);

    size_t memoryKeepFree = 250l*1024l*1024l;
    size_t memoryMin = 1100*1024*1024;
    size_t memoryToUse = min(memoryFree - memoryKeepFree, memoryMin );
    size_t numParts = 1;
    size_t maxNumSlicesInMemory = dim[2];

    std::cout << "Memory on device: " << memoryFree << " byte free (" << memoryTotal << " byte total), " << memoryRequired << " byte required." << std::endl;
    std::cout << "Using at most " << memoryToUse << " bytes." << std::endl;

    unsigned int overlap = 4*stencilRange;

    if (memoryToUse < memoryRequired)
    {
        //calculate memory required for overlapping slice of data for the blocks
        size_t numPointsSlice=dim[0]*dim[1];
        size_t sgSlice = 3*numPointsSlice*sizeof(double);
        size_t shesseSlice = 3*3*numPointsSlice*sizeof(double);
        size_t seValSlice = 3*numPointsSlice*sizeof(double);
        size_t seVectSlice = 3*numPointsSlice*sizeof(double);
        size_t slargeEVSlice = 3*numPointsSlice*sizeof(double);
        size_t sdSlice = numPointsSlice*sizeof(double);
        size_t sbitmasksSlice = 2*numPointsSlice*sizeof(uint);
        size_t memoryRequiredSlice = sdSlice + sgSlice + seValSlice + seVectSlice + slargeEVSlice + shesseSlice + sbitmasksSlice + sdim + snVT + stT;

        //add to required memory per part and test, if data needs to be split into more parts
        maxNumSlicesInMemory = memoryToUse / memoryRequiredSlice;
        numParts = (dim[2] - 1) / maxNumSlicesInMemory + 1;
        numParts += numParts*overlap/maxNumSlicesInMemory;
        while (dim[2] - (numParts-1) * (maxNumSlicesInMemory - overlap) > maxNumSlicesInMemory)
        {
            ++numParts;
        }
//        numParts = (memoryRequired - 1) / (memoryRequiredSlice * (maxNumSlicesInMemory - 1)) + 1; // maxNumSlicesInMemory - 1 to account for overlapping slice
        size_t memoryRequiredPart = memoryRequiredSlice * (maxNumSlicesInMemory - 1);
//        if (memoryToUse < memoryRequiredPart)
//        {
//            ++numParts;
//        }

        std::cout << "Data set needs to be split into " << numParts << " parts of " << memoryRequiredPart << " bytes + " << " bytes kept free." << std::endl;

//        cudaDeviceReset();
//        return false;
    }

    // double minEVal = -0.1;
    double *hesse_d;
    double *g_d, *d_d, *eVect_d, *evals_d, *largeEV_d;
    uint *bitmasks_d, *nVT_d, *tT_d;

    size_t numPointsPart=dim[0]*dim[1]*maxNumSlicesInMemory;
    size_t sgPart = 3*numPointsPart*sizeof(double);
    size_t shessePart = 3*3*numPointsPart*sizeof(double);
    size_t seValPart = 3*numPointsPart*sizeof(double);
    size_t seVectPart = 3*numPointsPart*sizeof(double);
    size_t slargeEVPart = 3*numPointsPart*sizeof(double);
    size_t sdPart = numPointsPart*sizeof(double);
    size_t sbitmasksPart = 2*numPointsPart*sizeof(uint);

//    double *verticesPart[numParts];
//    double *evalsPart[numParts];
//    double *pointDataPart[numParts];
    std::vector< std::vector<double> > verticesPart(numParts);
    std::vector< std::vector<double> > evalsPart(numParts);
    std::vector< std::vector<double> > pointDataPart(numParts);
//    size_t numVerticesPart[numParts];
//    size_t numActualPointsPart[numParts];
//    std::vector<size_t> numVerticesPart(numParts);
//    std::vector<size_t> numActualPointsPart(numParts);

    size_t numVerticesTotal = 0;

    gpuErrchk(cudaMalloc((void**)&d_d, sdPart));
    gpuErrchk(cudaMalloc((void**)&g_d, sgPart));
    gpuErrchk(cudaMalloc((void**)&evals_d, seValPart));
    gpuErrchk(cudaMalloc((void**)&eVect_d, seVectPart));
    gpuErrchk(cudaMalloc((void**)&largeEV_d, slargeEVPart));
    gpuErrchk(cudaMalloc((void**)&hesse_d, shessePart));
    gpuErrchk(cudaMalloc((void**)&bitmasks_d, sbitmasksPart));
    gpuErrchk(cudaMalloc((void**)&dim_d, sdim));
    gpuErrchk(cudaMalloc((void**)&nVT_d, snVT));
    gpuErrchk(cudaMalloc((void**)&tT_d, stT));

#ifdef CUDA_OUTPUT_STEPS
    std::cout << "    CUDA: Memory allocated." << std::endl;
#endif

    for (size_t i = 0; i < numParts; i++)
    {
//        verticesPart[i] = NULL;
//        evalsPart[i] = NULL;
//        pointDataPart[i] = NULL;
        cudaMemGetInfo(&memoryFree,&memoryTotal);
        std::cout << "    CUDA: Available device memory: " << memoryFree << " bytes. " << std::endl;

        size_t sliceOffset = i * (maxNumSlicesInMemory - overlap);
        //size_t sliceOffset = i * maxNumSlicesInMemory;
        size_t offsetIndex = sliceOffset * dim[0] * dim[1];
        size_t numSlices = i < numParts - 1 ? maxNumSlicesInMemory : dim[2] - sliceOffset;
        size_t numPoints = dim[0]*dim[1]*numSlices;
        //numActualPointsPart[i] = numPoints;
        size_t sd = numPoints*sizeof(double);
        //size_t sbitmasks = 2*numPoints*sizeof(uint);
        std::cout << "    CUDA: slice offset: " << sliceOffset << " number of slices " << numSlices << " number of points " << numPoints << std::endl;

        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );

        //uint *bitmasks_h = new uint[2*numPoints];
        gpuErrchk(cudaMemcpy(d_d, data + offsetIndex, sd, cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(nVT_d, numVertsTable, snVT, cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(tT_d, triTable, stT, cudaMemcpyHostToDevice));
        gpuErrchk(cudaMemcpy(dim_d, dim, sdim, cudaMemcpyHostToDevice));

        gpuErrchk(cudaMemset(g_d, 0, sgPart));
//        double *d_debug = new double[numPoints];
//        gpuErrchk(cudaMemcpy(d_d, d + offsetIndex, sd, cudaMemcpyHostToDevice));
//        gpuErrchk(cudaMemcpy(nVT_d, numVertsTable, snVT, cudaMemcpyHostToDevice));
//        gpuErrchk(cudaMemcpy(tT_d, triTable, stT, cudaMemcpyHostToDevice));
//        gpuErrchk(cudaMemcpy(dim_d, dim, sdim, cudaMemcpyHostToDevice));

    #ifdef CUDA_OUTPUT_STEPS
        std::cout << "    CUDA: Memory copied." << std::endl;
    #endif
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );

        uint blockSize = 128;
        int nBlocks = numPoints/blockSize + (numPoints%blockSize == 0?0:1);

        //first calculate gradients
        processPointGradients <<< nBlocks, blockSize >>> (numPoints, d_d, g_d, dim_d, stencilRange, valley);
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );

        //use those to calculate hessian matrix
        processPointHesse <<< nBlocks, blockSize >>> (numPoints, g_d, hesse_d, dim_d, stencilRange);
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );

        //use those to calculate eigenvalues
        processPointEigenVals <<< nBlocks, blockSize >>> (numPoints, hesse_d, evals_d, dim_d);
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );

        //use those to calculate eigenvectors
        processPointEigenVect <<< nBlocks, blockSize >>> (numPoints, hesse_d, evals_d, eVect_d, dim_d);
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );

        //use those to calculate an array of vectors giving the general orientation of eigenvectors per cube
        largestEigenVect <<< nBlocks, blockSize >>> (numPoints, eVect_d, largeEV_d, dim_d);
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );

        //generate bitmasks, arrays for prefix scan
        processEigenData <<< nBlocks, blockSize >>> (numPoints, dim_d, evals_d, eVect_d, g_d, largeEV_d, bitmasks_d, nVT_d, d_d, featureThreshold);
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );

        //gpuErrchk(cudaMemcpy(bitmasks_h, bitmasks_d, sbitmasks, cudaMemcpyDeviceToHost));

        //uint *pfA_h = new uint[2*numPoints];
        uint *pfA_d;
        size_t pfAB = 2*numPoints*sizeof(uint);
        gpuErrchk(cudaMalloc( (void**)&pfA_d, pfAB));

        thrust::device_ptr<uint> prefixArray = thrust::device_pointer_cast(pfA_d);
        thrust::device_ptr<uint> bitmaskArray = thrust::device_pointer_cast(bitmasks_d);

        int addOne = 0;
        int addTris = 0;
        if (bitmaskArray[2*numPoints-1]) {
            addOne = 1;                    //doing exclusive scan, if last element is set, we need one more space
            addTris = bitmaskArray[2*numPoints-1];  //and the tris of course
        }

//        thrust::exclusive_scan(&bitmasks_h[numPoints], bitmasks_h+2*numPoints, &pfA_h[numPoints]);  //index of next cubes
//        thrust::exclusive_scan(&bitmasks_h[0], bitmasks_h+numPoints, &pfA_h[0]);  //index of next cubes
        thrust::exclusive_scan(bitmaskArray + numPoints, bitmaskArray + 2*numPoints, prefixArray + numPoints);  //index of next cubes
        thrust::exclusive_scan(bitmaskArray, bitmaskArray + numPoints, prefixArray);  //index of next cubes

        uint numIsoSurfaceCubes=prefixArray[2*numPoints-1]+addOne;
        uint numVertices=prefixArray[numPoints-1]+addTris;

        if (numIsoSurfaceCubes == 0)
        {
            //numVerticesPart[i] = 0;
            //evalsPart[i] = new double[numPoints*3];
            evalsPart[i].resize(numPoints*3);
            gpuErrchk(cudaMemcpy(evalsPart[i].data(), evals_d, numPoints*3*sizeof(double), cudaMemcpyDeviceToHost));
            continue;
        }
        numVerticesTotal += numVertices;


    #ifdef CUDA_OUTPUT_STEPS
        cout << "    CUDA: numIsoSurfaceCubes: " << numIsoSurfaceCubes << endl;
        cout << "    CUDA: numVertices: " << numVertices << endl;
    #endif

        size_t siA, svA, sgA, sizePointData;
        uint *iA_d;
        double *vA_d;//, *vA_h = new double[numVertices*3];
        verticesPart[i].resize(numVertices*3);
        double *gA_d;
        //double *evals_h = new double[numPoints*3];
        evalsPart[i].resize(numPoints*3);
        double *pointData_d;//, *pointData_h = new double[numVertices];
        pointDataPart[i].resize(numVertices);
        siA = numIsoSurfaceCubes*sizeof(uint);
        svA = numVertices*3*sizeof(double);
        sgA = numVertices*3*sizeof(double);
        sizePointData = numVertices*sizeof(double);

        gpuErrchk(cudaMalloc( (void**)&iA_d,  siA));
        gpuErrchk(cudaMalloc( (void**)&vA_d,  svA));
        gpuErrchk(cudaMalloc( (void**)&gA_d,  sgA));
        gpuErrchk(cudaMalloc( (void**)&pointData_d, sizePointData));

        //gpuErrchk(cudaMemcpy( pfA_d, pfA_h, pfAB, cudaMemcpyHostToDevice));   //copy prefix array for second pass
        getIsoIndices <<< nBlocks, blockSize >>> (numPoints, iA_d, bitmasks_d, pfA_d);
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
//        std::vector<unsigned int> cellIndices(numIsoSurfaceCubes);
//        gpuErrchk(cudaMemcpy(cellIndices.data(), iA_d, siA, cudaMemcpyDeviceToHost));

//        std::set<unsigned int> nodeIndicesTmp;
//        for (std::vector<unsigned int>::iterator it = cellIndices.begin(); it != cellIndices.end(); it++)
//        {
//            nodeIndicesTmp.insert(*it);
//            nodeIndicesTmp.insert(*it + 1);
//            nodeIndicesTmp.insert(*it + dim[0]);
//            nodeIndicesTmp.insert(*it + dim[0] + 1);
//            nodeIndicesTmp.insert(*it + dim[1]*dim[0]);
//            nodeIndicesTmp.insert(*it + dim[1]*dim[0] + 1);
//            nodeIndicesTmp.insert(*it + dim[1]*dim[0] + dim[0]);
//            nodeIndicesTmp.insert(*it + dim[1]*dim[0] + dim[0] + 1);
//        }

        nBlocks = numIsoSurfaceCubes/blockSize + (numIsoSurfaceCubes%blockSize == 0?0:1);
        double3 originTmp;
        originTmp.x = origin[0];
        originTmp.y = origin[1];
        originTmp.z = origin[2] + sliceOffset;
        double3 spacingTmp;
        spacingTmp.x = spacing[0];
        spacingTmp.y = spacing[1];
        spacingTmp.z = spacing[2];
        processRidgeCubesData <<< nBlocks, blockSize >>> (numIsoSurfaceCubes, dim_d, iA_d, pfA_d, eVect_d, g_d, largeEV_d, vA_d, gA_d, nVT_d, tT_d, bitmasks_d, d_d, pointData_d, originTmp, spacingTmp);
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );

        //gpuErrchk(cudaMemcpy(vA_h, vA_d, svA, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(verticesPart[i].data(), vA_d, svA, cudaMemcpyDeviceToHost));
        //gpuErrchk(cudaMemcpy(evals_h, evals_d, sgA, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(evalsPart[i].data(), evals_d, evalsPart[i].size()*sizeof(double), cudaMemcpyDeviceToHost));
        //gpuErrchk(cudaMemcpy(pointData_h, pointData_d, sizePointData, cudaMemcpyDeviceToHost));
        gpuErrchk(cudaMemcpy(pointDataPart[i].data(), pointData_d, sizePointData, cudaMemcpyDeviceToHost));

        gpuErrchk(cudaFree(iA_d));
        gpuErrchk(cudaFree(vA_d));
        gpuErrchk(cudaFree(gA_d));
        gpuErrchk(cudaFree(pfA_d));
        gpuErrchk(cudaFree(pointData_d));
        cout << "vram released" << endl;

        //delete [] bitmasks_h;
        //delete [] pfA_h;
        cout << "host mem released" << endl;

        //numVerticesPart[i] = numVertices;
        //verticesPart[i] = vA_h;
        //evalsPart[i] = evals_h;
        //pointDataPart[i] = pointData_h;
//        *numVertices_out = numVertices;
//        *vertices_out = vA_h;
//        *gradients_out = gA_h;
//        *pointData_out = pointData_h;
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
    }
    gpuErrchk(cudaFree(hesse_d));
    gpuErrchk(cudaFree(evals_d));

    gpuErrchk(cudaFree(dim_d));
    gpuErrchk(cudaFree(d_d));
    gpuErrchk(cudaFree(g_d));
    gpuErrchk(cudaFree(eVect_d));
    gpuErrchk(cudaFree(largeEV_d));


    double *vertices = new double[numVerticesTotal*3];
    double *evals = new double[numPointsTotal*3];
    double *pointData = new double[numVerticesTotal];

    if (!(vertices && evals && pointData))
    {
        return false;
    }

    size_t indexOffset = 0;
    size_t indexOffsetEval = 0;
    //copy everything together
    for (size_t i = 0; i < numParts; i++)
    {
        //copy
        if (!verticesPart[i].empty())
        {
            memcpy(vertices + indexOffset*3, verticesPart[i].data(), verticesPart[i].size()*sizeof(double));
            memcpy(pointData + indexOffset, pointDataPart[i].data(), pointDataPart[i].size()*sizeof(double));
            indexOffset += pointDataPart[i].size();
        }

        //account for overlapping slice
        size_t numEvalsSlice = 3*(dim[0]*dim[1]);
        size_t sliceOffset = i == 0 ? 0 : overlap/2;
        size_t numEvalsToCopy = (i > 0 && i < numParts - 1) ? evalsPart[i].size() - overlap*numEvalsSlice : evalsPart[i].size() - (overlap/2)*numEvalsSlice;
        if (numParts == 1) numEvalsToCopy = evalsPart[i].size();
        memcpy(evals + indexOffsetEval, evalsPart[i].data() + sliceOffset*numEvalsSlice, numEvalsToCopy*sizeof(double));
        indexOffsetEval += numEvalsToCopy;
    }

//    size_t indexOffset = 0;
//    size_t indexOffsetEval = 0;
//    //copy everything together
//    for (size_t i = 0; i < numParts; i++)
//    {
//        //copy
//        if (numVerticesPart[i] > 0)
//        {
//            memcpy(vertices + indexOffset*3, verticesPart[i], 3*numVerticesPart[i]*sizeof(double));
//            memcpy(pointData + indexOffset, pointDataPart[i], numVerticesPart[i]*sizeof(double));
//            indexOffset += numVerticesPart[i];
//            if (verticesPart[i]) delete[] verticesPart[i];
//            if (pointDataPart[i]) delete[] pointDataPart[i];
//        }

//        //account for overlapping slice
//        size_t numEvalsToCopy = i < numParts - 1 ? 3*(numActualPointsPart[i] - dim[0]*dim[1]) : 3*(numActualPointsPart[i]);
//        memcpy(evals + indexOffsetEval*3, evalsPart[i], numEvalsToCopy*sizeof(double));
//        indexOffsetEval += numEvalsToCopy;
//        //delete
//        if (evalsPart[i]) delete[] evalsPart[i];
//    }

    *numVertices_out = numVerticesTotal;
    *vertices_out = vertices;
    *evals_out = evals;
    *pointData_out = pointData;

    cudaDeviceReset();

    return true;
}
/* old stuff
 *
bool cudaRidgeDataOld(uint dim[3], const double* d, double** gradients_out, uint* numVertices_out, double** vertices_out, double** pointData_out, double featureThreshold, double origin[3], double spacing[3]) {
#ifdef CUDA_OUTPUT_STEPS
    std::cout << "    CUDA: Starting cuda method..." << std::endl;
#endif

    int deviceCount = 0;
    gpuErrchk(cudaGetDeviceCount(&deviceCount));
    std::cout << "    CUDA: number of devices is " << deviceCount << std::endl;

    if (cudaSetDevice(0) != cudaSuccess) {
    std::cout << "    CUDA: Could not set cuda device." << std::endl;
    }
    else {
#ifdef CUDA_OUTPUT_STEPS
    std::cout << "    CUDA: Cuda device set." << std::endl;
#endif
    }
    size_t numPoints=dim[0]*dim[1]*dim[2];
    size_t sg = 3*numPoints*sizeof(double);
    size_t shesse = 3*3*numPoints*sizeof(double);
    size_t seVal = numPoints*sizeof(double);
    size_t seVect = 3*numPoints*sizeof(double);
    size_t slargeEV = 3*numPoints*sizeof(double);
    size_t sd = numPoints*sizeof(double);
    size_t sdim = 3*sizeof(uint);
    size_t sdebug = 9*sizeof(double);
    size_t sbitmasks = 2*numPoints*sizeof(uint);
    size_t snVT = 256*sizeof(uint);
    size_t stT = 256*16*sizeof(uint);
    uint *dim_d;

    // double minEVal = -0.1;
    double *g_h = new double[3*numPoints];
    double *eVect_h = new double[3*numPoints];
    double *largeEV_h = new double[3*numPoints];
    double *eVal_h = new double[numPoints];
    double *hesse_h = new double[3*3*numPoints];
    uint *bitmasks_h = new uint[2*numPoints];
    double debug_h[9];
    double *hesse_d;
    double *g_d, *d_d, *eVect_d, *eVal_d, *debug_d, *largeEV_d;
    uint *bitmasks_d, *nVT_d, *tT_d;

    size_t memoryRequired = 128*sizeof(double) + sd + sg + seVal + seVect + slargeEV + shesse + sdim + sbitmasks + snVT + stT + sdebug;
    size_t memoryFree, memoryTotal;
    cudaMemGetInfo(&memoryFree,&memoryTotal);

    size_t memoryKeepFree = 25*1024*1024;
    size_t memoryToUse = memoryFree - memoryKeepFree;

    if (memoryToUse < memoryRequired)
    {
        size_t numParts = (memoryRequired - 1) / memoryToUse + 1;
        std::cout << "Not enough available memory on device. " << memoryFree << " byte free (" << memoryTotal << " byte total), " << memoryRequired << " + " << memoryKeepFree << " byte required." << std::endl;
        std::cout << "Data set needs to be split into " << numParts << " parts." << std::endl;
        return false;
    }

    gpuErrchk(cudaMalloc((void**)&d_d, sd));
    gpuErrchk(cudaMalloc((void**)&g_d, sg));
    gpuErrchk(cudaMalloc((void**)&eVal_d, seVal));
    gpuErrchk(cudaMalloc((void**)&eVect_d, seVect));
    gpuErrchk(cudaMalloc((void**)&largeEV_d, slargeEV));
    gpuErrchk(cudaMalloc((void**)&hesse_d, shesse));
    gpuErrchk(cudaMalloc((void**)&dim_d, sdim));
    gpuErrchk(cudaMalloc((void**)&bitmasks_d, sbitmasks));
    gpuErrchk(cudaMalloc((void**)&nVT_d, snVT));
    gpuErrchk(cudaMalloc((void**)&tT_d, stT));
    gpuErrchk(cudaMalloc((void**)&debug_d, sdebug));

#ifdef CUDA_OUTPUT_STEPS
    std::cout << "    CUDA: Memory allocated." << std::endl;
#endif
    gpuErrchk(cudaMemcpy(d_d, d, sd, cudaMemcpyHostToDevice));
    // std::cout << "    CUDA: Data copied" << std::endl;
    gpuErrchk(cudaMemcpy(nVT_d, numVertsTable, snVT, cudaMemcpyHostToDevice));
    // std::cout << "    CUDA: numVertsTable copied." << std::endl;
    gpuErrchk(cudaMemcpy(tT_d, triTable, stT, cudaMemcpyHostToDevice));
    // std::cout << "    CUDA: TriTable copied." << std::endl;
    gpuErrchk(cudaMemcpy(dim_d, dim, sdim, cudaMemcpyHostToDevice));
    // std::cout << "    CUDA: Dim array copied." << std::endl;
#ifdef CUDA_OUTPUT_STEPS
    std::cout << "    CUDA: Memory copied." << std::endl;
#endif

    uint blockSize = 128;
    int nBlocks = numPoints/blockSize + (numPoints%blockSize == 0?0:1);

    //first calculate gradients
    processPointGradients <<< nBlocks, blockSize >>> (numPoints, d_d, g_d, dim_d);
    gpuErrchk(cudaMemcpy(g_h, g_d, sg, cudaMemcpyDeviceToHost));

    //use those to calculate hessian matrix
    processPointHesse <<< nBlocks, blockSize >>> (numPoints, g_d, hesse_d, dim_d);
    gpuErrchk(cudaMemcpy(hesse_h, hesse_d, shesse, cudaMemcpyDeviceToHost));

    //use those to calculate eigenvalues
    processPointEigenVals <<< nBlocks, blockSize >>> (numPoints, hesse_d, eVal_d, dim_d);
    gpuErrchk(cudaMemcpy(eVal_h, eVal_d, seVal, cudaMemcpyDeviceToHost));

    //use those to calculate eigenvectors
    processPointEigenVect <<< nBlocks, blockSize >>> (numPoints, hesse_d, eVal_d, eVect_d, dim_d);
    gpuErrchk(cudaMemcpy(eVect_h, eVect_d, seVect, cudaMemcpyDeviceToHost));

    //use those to calculate an array of vectors giving the general orientation of eigenvectors per cube
    largestEigenVect <<< nBlocks, blockSize >>> (numPoints, eVect_d, largeEV_d, dim_d);
    gpuErrchk(cudaMemcpy(largeEV_h, largeEV_d, slargeEV, cudaMemcpyDeviceToHost));

    //generate bitmasks, arrays for prefix scan
    processEigenData <<< nBlocks, blockSize >>> (numPoints, dim_d, eVal_d, eVect_d, g_d, largeEV_d, bitmasks_d, nVT_d, d_d, featureThreshold);
    gpuErrchk(cudaMemcpy(bitmasks_h, bitmasks_d, sbitmasks, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(eVal_h, eVal_d, seVal, cudaMemcpyDeviceToHost));

    gpuErrchk(cudaFree(hesse_d));
    gpuErrchk(cudaFree(eVal_d));

    int addOne = 0;
    int addTris = 0;
    if (bitmasks_h[2*numPoints-1]) {
    addOne = 1;                    //doing exclusive scan, if last element is set, we need one more space
    addTris = bitmasks_h[2*numPoints-1];  //and the tris of course
    }

    uint *pfA_h = new uint[2*numPoints];

    thrust::exclusive_scan(&bitmasks_h[numPoints], bitmasks_h+2*numPoints, &pfA_h[numPoints]);  //index of next cubes
    thrust::exclusive_scan(&bitmasks_h[0], bitmasks_h+numPoints, &pfA_h[0]);  //index of next cubes

    uint numIsoSurfaceCubes=pfA_h[2*numPoints-1]+addOne;
    uint numVertices=pfA_h[numPoints-1]+addTris;


#ifdef CUDA_OUTPUT_STEPS
    cout << "    CUDA: numIsoSurfaceCubes: " << numIsoSurfaceCubes << endl;
    cout << "    CUDA: numVertices: " << numVertices << endl;
#endif

    size_t siA, svA, sgA, pfAB, sizePointData;
    uint *pfA_d, *iA_d, *iA_h = new uint[numIsoSurfaceCubes];
    double *vA_d, *vA_h = new double[numVertices*3];
    double *gA_d, *gA_h = new double[numVertices*3];
    double *pointData_d, *pointData_h = new double[numVertices];
    siA = numIsoSurfaceCubes*sizeof(uint);
    svA = numVertices*3*sizeof(double);
    sgA = numVertices*3*sizeof(double);
    sizePointData = numVertices*sizeof(double);
    pfAB = 2*numPoints*sizeof(uint);

    gpuErrchk(cudaMalloc( (void**)&iA_d,  siA));
    gpuErrchk(cudaMalloc( (void**)&vA_d,  svA));
    gpuErrchk(cudaMalloc( (void**)&gA_d,  sgA));
    gpuErrchk(cudaMalloc( (void**)&pfA_d, pfAB));
    gpuErrchk(cudaMalloc( (void**)&pointData_d, sizePointData));

    gpuErrchk(cudaMemcpy( pfA_d, pfA_h, pfAB, cudaMemcpyHostToDevice));   //copy prefix array for second pass
    getIsoIndices <<< nBlocks, blockSize >>> (numPoints, iA_d, bitmasks_d, pfA_d);
    std::vector<unsigned int> cellIndices(numIsoSurfaceCubes);
    gpuErrchk(cudaMemcpy(cellIndices.data(), iA_d, siA, cudaMemcpyDeviceToHost));

    std::set<unsigned int> nodeIndicesTmp;
    for (std::vector<unsigned int>::iterator it = cellIndices.begin(); it != cellIndices.end(); it++)
    {
        nodeIndicesTmp.insert(*it);
        nodeIndicesTmp.insert(*it + 1);
        nodeIndicesTmp.insert(*it + dim[0]);
        nodeIndicesTmp.insert(*it + dim[0] + 1);
        nodeIndicesTmp.insert(*it + dim[1]*dim[0]);
        nodeIndicesTmp.insert(*it + dim[1]*dim[0] + 1);
        nodeIndicesTmp.insert(*it + dim[1]*dim[0] + dim[0]);
        nodeIndicesTmp.insert(*it + dim[1]*dim[0] + dim[0] + 1);
    }

    nBlocks = numIsoSurfaceCubes/blockSize + (numIsoSurfaceCubes%blockSize == 0?0:1);
    double3 originTmp;
    originTmp.x = origin[0];
    originTmp.y = origin[1];
    originTmp.z = origin[2];
    double3 spacingTmp;
    spacingTmp.x = spacing[0];
    spacingTmp.y = spacing[1];
    spacingTmp.z = spacing[2];
    processRidgeCubesData <<< nBlocks, blockSize >>> (numIsoSurfaceCubes, dim_d, iA_d, pfA_d, eVect_d, g_d, largeEV_d, vA_d, gA_d, nVT_d, tT_d, bitmasks_d, d_d, pointData_d, originTmp, spacingTmp);

    gpuErrchk(cudaMemcpy(vA_h, vA_d, svA, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(gA_h, gA_d, sgA, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(debug_h, debug_d, sdebug, cudaMemcpyDeviceToHost));
    gpuErrchk(cudaMemcpy(pointData_h, pointData_d, sizePointData, cudaMemcpyDeviceToHost));

    gpuErrchk(cudaFree(pointData_d));
    gpuErrchk(cudaFree(dim_d));
    gpuErrchk(cudaFree(d_d));
    gpuErrchk(cudaFree(g_d));
    gpuErrchk(cudaFree(eVect_d));
    gpuErrchk(cudaFree(largeEV_d));
    gpuErrchk(cudaFree(iA_d));
    gpuErrchk(cudaFree(vA_d));
    gpuErrchk(cudaFree(gA_d));
    gpuErrchk(cudaFree(pfA_d));
    gpuErrchk(cudaFree(debug_d));
    cout << "vram released" << endl;

    delete [] eVect_h;
    delete [] largeEV_h;
    delete [] eVal_h;
    delete [] hesse_h;
    delete [] bitmasks_h;
    delete [] pfA_h;
    delete [] iA_h;
    cout << "host mem released" << endl;

    *numVertices_out = numVertices;
    *vertices_out = vA_h;
    *gradients_out = gA_h;
    *pointData_out = pointData_h;
    return true;
}

bool cudaGradient(const double* d, uint dim[3], double** gradients) {
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
   
    uint numPoints=dim[0]*dim[1]*dim[2];
    uint sg = 3*numPoints*sizeof(double);
    uint sd = numPoints*sizeof(double);
    uint sdim = 3*sizeof(uint);
    double *g_h = new double[3*numPoints];
    uint *dim_d;
    double *g_d, *d_d;


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

    //calculate gradients
    processPointGradients <<< nBlocks, blockSize >>> (numPoints, d_d, g_d, dim_d);
    gpuErrchk(cudaMemcpy(g_h, g_d, sg, cudaMemcpyDeviceToHost));
    
    gpuErrchk(cudaFree(dim_d));
    gpuErrchk(cudaFree(d_d));
    gpuErrchk(cudaFree(g_d));

    *gradients = g_h;
    return true;
}

bool cudaIsosurface(const double c, const double* in_data, double* in_grads, int in_size[3], uint* out_numVerts, double** out_verts, double** out_grads)
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
    double *vA, *vA_h = new double[numVertices*3];
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

    gpuErrchk(cudaMemcpy(vA_h, vA, vAB, cudaMemcpyDeviceToHost ));
    gpuErrchk(cudaMemcpy(gradientArray, gA, gAB, cudaMemcpyDeviceToHost ));

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

    *out_verts = vA_h;
    *out_grads = gradientArray;
    *out_numVerts = numVertices;

#ifdef CUDA_OUTPUT_STEPS
    cout << "    CUDA: Cuda calculation done." << endl;
#endif
    return true;
}
*/
