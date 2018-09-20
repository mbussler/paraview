#ifndef __CUDARIDGE_H_
#define __CUDARIDGE_H_

#include "undef_atomics_int128.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <thrust/scan.h>
#include <cuda_runtime.h>
#include <cublas.h>

using namespace std;

inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
    if (code != cudaSuccess) {
    std::cout << "GPUassert: " << cudaGetErrorString(code) << "(in " << file << ":" << line << ")" << std::endl;
	if (abort) exit(code);
    }
}

//used to output some information while the cuda host method is run
#define CUDA_OUTPUT_STEPS

#endif
