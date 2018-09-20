#ifndef __KERNELS_CU__
#define __KERNELS_CU__

#include <climits>
#include <cfloat>
#include "std_kernel.cu"

template <typename T>
__device__ void inline myswap(T &a, T &b)
{
    T c(a); a=b; b=c;
}

template <typename T>
inline __device__ void sort3(T &a, T &b, T &c)
{
    if (a < b) {
        if (c < a) myswap(a,c);
    } else {
    if (b < c) myswap(a,b);
        else myswap(a,c);
    }
    if(c<b) myswap(b,c);
}

template <typename T>
inline __device__ void sort3(T l[3])
{
    sort3(l[0], l[1], l[2]);
}

__global__ void processPointEigenVals(uint maxIndex, double* h, double* e, uint dim[3]) {
    uint index=blockIdx.x*blockDim.x + threadIdx.x;
    if (index < maxIndex) {

	uint hmIdx = index*9;
	double polyCoeff[3], lambda[3];
	GetInvariants(&h[hmIdx], polyCoeff);
	
	if (CalculateRoots(polyCoeff, lambda)) {
//	    double smallest;
//	    lambda[0] < lambda[1] ?
//		(lambda[0] < lambda[2] ? smallest=lambda[0] : smallest=lambda[2]):
//		(lambda[1] < lambda[2] ? smallest=lambda[1] : smallest=lambda[2]);
        sort3(lambda);

        e[index*3+0] = lambda[0];
        e[index*3+1] = lambda[1];
        e[index*3+2] = lambda[2];
	}
    else {
        e[index*3+0] = 0.0;
        e[index*3+1] = 0.0;
        e[index*3+2] = 0.0;
	}
    }    
}

__global__ void processPointEigenVect(uint maxIndex, double* h, double* e, double* evect, uint dim[3]) {
    uint index=blockIdx.x*blockDim.x + threadIdx.x;
    if (index < maxIndex) {
	uint hmIdx = index*9;
    double lambda = e[index*3];
	CalculateEigenvector(&h[hmIdx], lambda, &evect[index*3]);
    }
}

__global__ void largestEigenVect(uint maxIndex, double* evect, double* lEV, uint dim[3]) {
    __shared__ uint ofst[8];
    ofst[0] = 0;
    ofst[1] = 1;
    ofst[2] = dim[0];
    ofst[3] = 1+dim[0];
    ofst[4] = 0+dim[0]*dim[1];
    ofst[5] = 1+dim[0]*dim[1];
    ofst[6] = dim[0]+dim[0]*dim[1];
    ofst[7] = 1+dim[0]+dim[0]*dim[1];
    
    uint index=blockIdx.x*blockDim.x + threadIdx.x;
    double C[9];
    if (index < maxIndex) {
	uint i,j,k;
	k=index/(dim[0]*dim[1]);
	j=index%(dim[0]*dim[1])/dim[0];
	i=(index%(dim[0]*dim[1]))%dim[0];
	if (i<dim[0]-1 && j<dim[1]-1 && k<dim[2]-1) {
	    const uint numNodes = 8;
	    double evSet[3][numNodes*2];

	    for (int n=0; n<numNodes; ++n) {
		evSet[0][2*n+0] = evect[(index+ofst[n])*3+0];
		evSet[1][2*n+0] = evect[(index+ofst[n])*3+1];
		evSet[2][2*n+0] = evect[(index+ofst[n])*3+2];
		evSet[0][2*n+1] = -evect[(index+ofst[n])*3+0];
		evSet[1][2*n+1] = -evect[(index+ofst[n])*3+1];
		evSet[2][2*n+1] = -evect[(index+ofst[n])*3+2];
	    }
	
	    //Covariance Matrix
	    for(int j=0; j<3; ++j) {
		for(int i=0; i<3; ++i) {
		    double sum = 0.0;
		    for(int n=0; n<8; ++n){
              sum += evSet[i][n] * evSet[j][n];
		    }
		    C[j*3+i] = sum/16.0;
		}
	    }
    
	    double pqr[3];
	    double lambda[3];
	    double largeEigenV[3];
	    GetInvariants(C, pqr);
	    CalculateRoots(pqr, lambda);

	    //get index of largest eigenvalue
	    int largest = getLargest(lambda);
	    CalculateEigenvector(C, lambda[largest], largeEigenV);
	    lEV[index*3+0]=largeEigenV[0];
	    lEV[index*3+1]=largeEigenV[1];
	    lEV[index*3+2]=largeEigenV[2];
	}
    }
}

__global__ void processEigenData(int maxIndex, uint dim[3], double* eVal, double* ev, double* g, double* lEV, uint* bitmasks, uint* nVT, double* d, double featureThreshold){
    int index=blockIdx.x*blockDim.x + threadIdx.x;

    __shared__ uint ofst[8];
    ofst[0] = 0;
    ofst[1] = 1;
    ofst[2] = 1+dim[0];
    ofst[3] = dim[0];
    ofst[4] = 0+dim[0]*dim[1];
    ofst[5] = 1+dim[0]*dim[1];
    ofst[6] = 1+dim[0]+dim[0]*dim[1];
    ofst[7] = dim[0]+dim[0]*dim[1];

    if (index < maxIndex) {
	uint i,j,k;
	double sumD=0;
	k=index/(dim[0]*dim[1]);
	j=index%(dim[0]*dim[1])/dim[0];
	i=(index%(dim[0]*dim[1]))%dim[0];
	if(i !=  dim[0]-1 && j != dim[1]-1 && k != dim[2]-1) {
	    const uint numNodes = 8;
	    double gradDotEV[numNodes];
	    double orientedEV[numNodes*3];
	    double grads[numNodes*3];
	    for (int n=0; n<numNodes; ++n) {
		sumD += d[(index+ofst[n])];
		orientedEV[n*3+0] = ev[(index+ofst[n])*3+0];
		orientedEV[n*3+1] = ev[(index+ofst[n])*3+1];
		orientedEV[n*3+2] = ev[(index+ofst[n])*3+2];
		grads[n*3+0] = g[(index+ofst[n])*3+0];
		grads[n*3+1] = g[(index+ofst[n])*3+1];
		grads[n*3+2] = g[(index+ofst[n])*3+2];
	        if (dotProduct(&lEV[index*3], &orientedEV[n*3]) < 0) {
		    orientedEV[n*3+0] *= -1;
		    orientedEV[n*3+1] *= -1;
		    orientedEV[n*3+2] *= -1;
		}
		gradDotEV[n] = dotProduct(&grads[n*3], &orientedEV[n*3]);
        }
	    uint bitmask = 0;
	    bitmask |= ((gradDotEV[0] < 0) << 0);
	    bitmask |= ((gradDotEV[1] < 0) << 1);
	    bitmask |= ((gradDotEV[2] < 0) << 2);
	    bitmask |= ((gradDotEV[3] < 0) << 3);
	    bitmask |= ((gradDotEV[4] < 0) << 4);
	    bitmask |= ((gradDotEV[5] < 0) << 5);
	    bitmask |= ((gradDotEV[6] < 0) << 6);
	    bitmask |= ((gradDotEV[7] < 0) << 7);
	    int flagCube = ((bitmask) && (~bitmask & 0xFF));
	    if (flagCube) {
        if (eVal[index*3] < featureThreshold) {
		//if (true) {
		    // bitmasks[index] = nVT[bitmask];
		    // bitmasks[index+maxIndex] = 1;
		    bitmasks[index] = nVT[bitmask];
		    bitmasks[index+maxIndex] = 1;
		    // double eValInterp = 0;
		    // //s 1.5 p 5.0 -> -0.03
		    // //s 1.25 p 3.0 -> -0.0045
		    // for (int edg=0; edg<12; ++edg) {
		    // 	if (nVT[bitmask] & (0x1 << edg)) {
		    // 	    eValInterp = max(eVal[index+eI[edg*2]], eVal[index+eI[edg*2+1]]);
		    // 	    // gradInterp[0] = g[(index+eI[edg*2])*3+0]
		    // 	    if (!(eValInterp < -0.015 && eValInterp > -0.03)) {
		    // 		bitmasks[index] = nVT[bitmask];
		    // 		bitmasks[index+maxIndex] = 1;
		    // 	    }
		    // 	    else {
		    // 		bitmasks[index] = 0;
		    // 		bitmasks[index+maxIndex] = 0;
		    // 	    }
		    // 	}
		    // }
		}
		else {
		    bitmasks[index] = 0;
		    bitmasks[index+maxIndex] = 0;
		}
	    }
	    else {
		bitmasks[index] = 0;
		bitmasks[index+maxIndex] = 0;
	    }
	}
	else {
	    bitmasks[index] = 0;
	    bitmasks[index+maxIndex] = 0;
	}
    }
}

__global__ void processPointHesseOld(uint maxIndex, double* g, double* h, uint dim[3]) {
    uint index=blockIdx.x*blockDim.x + threadIdx.x;
    if (index < maxIndex) {
	uint k=index/(dim[0]*dim[1]);
	uint j=index%(dim[0]*dim[1])/dim[0];
	uint i=(index%(dim[0]*dim[1]))%dim[0];
	if (i>0 && i<dim[0]-1) {         //central differences for inner cubes
	    uint hrow =(index)*9+0;
	    uint neg = (index-1)*3;
	    uint pos = (index+1)*3;
	    h[hrow+0]=(g[pos+0]-g[neg+0])/2.0;
	    h[hrow+1]=(g[pos+1]-g[neg+1])/2.0;
	    h[hrow+2]=(g[pos+2]-g[neg+2])/2.0;
	}
	else {
	    if (i == 0) {                //forward differences for bounding cubes
		uint hrow =(index)*9+0;
		uint neg = (index+0)*3;
		uint pos = (index+1)*3;
		h[hrow+0]=(g[pos+0]-g[neg+0]);
		h[hrow+1]=(g[pos+1]-g[neg+1]);
		h[hrow+2]=(g[pos+2]-g[neg+2]);
	    }
	    if (i == dim[0]-1) {         //forward differences for bounding cubes
		uint hrow =(index)*9+0;
		uint neg = (index-1)*3;
		uint pos = (index+0)*3;
		h[hrow+0]=(g[pos+0]-g[neg+0]);
		h[hrow+1]=(g[pos+1]-g[neg+1]);
		h[hrow+2]=(g[pos+2]-g[neg+2]);
	    }
	}
	if (j>0 && j<dim[1]-1) {
	    uint hrow =(index)*9+3;
	    uint neg = (index-dim[0])*3;
	    uint pos = (index+dim[0])*3;
	    h[hrow+0]=(g[pos+0]-g[neg+0])/2.0;
	    h[hrow+1]=(g[pos+1]-g[neg+1])/2.0;
	    h[hrow+2]=(g[pos+2]-g[neg+2])/2.0;
	}
	else {
	    if (j == 0) {                //forward differences for bounding cubes
		uint hrow =(index)*9+3;
		uint neg = (index-  0   )*3;
		uint pos = (index+dim[0])*3;
		h[hrow+0]=(g[pos+0]-g[neg+0]);
		h[hrow+1]=(g[pos+1]-g[neg+1]);
		h[hrow+2]=(g[pos+2]-g[neg+2]);
	    }
	    if (j == dim[1]-1) {         //forward differences for bounding cubes
		uint hrow =(index)*9+3;
		uint pos = (index-dim[0])*3;
		uint neg = (index+  0   )*3;
		h[hrow+0]=(g[pos+0]-g[neg+0]);
		h[hrow+1]=(g[pos+1]-g[neg+1]);
		h[hrow+2]=(g[pos+2]-g[neg+2]);
	    }
	}
	if (k>0 && k<dim[2]-1) {
	    uint hrow =(index)*9+6;
	    uint neg = (index-dim[0]*dim[1])*3;
	    uint pos = (index+dim[0]*dim[1])*3;
	    h[hrow+0]=(g[pos+0]-g[neg+0])/2.0;
	    h[hrow+1]=(g[pos+1]-g[neg+1])/2.0;
	    h[hrow+2]=(g[pos+2]-g[neg+2])/2.0;
	}
	else {
	    if (k == 0) {         //forward differences for bounding cubes
		uint hrow =(index)*9+6;
		uint neg = (index-      0      )*3;
		uint pos = (index+dim[0]*dim[1])*3;
		h[hrow+0]=(g[pos+0]-g[neg+0]);
		h[hrow+1]=(g[pos+1]-g[neg+1]);
		h[hrow+2]=(g[pos+2]-g[neg+2]);
	    }
	    if (k == dim[2]-1) {                //forward differences for bounding cubes
		uint hrow =(index)*9+6;
		uint neg = (index-dim[0]*dim[1])*3;
		uint pos = (index+      0      )*3;
		h[hrow+0]=(g[pos+0]-g[neg+0]);
		h[hrow+1]=(g[pos+1]-g[neg+1]);
		h[hrow+2]=(g[pos+2]-g[neg+2]);
	    }
	}
    }    
}

__global__ void processPointHesse(uint maxIndex, double* gradient, double* hessian, uint dim[3], unsigned int stencilRange)
{
    uint index=blockIdx.x*blockDim.x + threadIdx.x;
    if (index < maxIndex) {
        uint z=index/(dim[0]*dim[1]);
        uint y=index%(dim[0]*dim[1])/dim[0];
        uint x=(index%(dim[0]*dim[1]))%dim[0];

        mat3 m = { {0}, {0}, {0} };
        vec3 rhsX = { 0 };
        vec3 rhsY = { 0 };
        vec3 rhsZ = { 0 };

        //vec3 xyzi = {x, y, z}; // ###, vi;
        double gX = gradient[index*3+0];
        double gY = gradient[index*3+1];
        double gZ = gradient[index*3+2];

        // compute all neighbors inside level 'range'
        int width = 2*stencilRange+1;
        int neighCnt = width*width*width;//computeNodeNeighborsN(i, range, neighborsN);

        int effNeighCnt = 0;
        for (int k = 0; k < neighCnt; k++) {
            int offsetX = (k%width) - stencilRange;
            int offsetY = ((k/width)%width) - stencilRange;
            int offsetZ = (k/(width*width)) - stencilRange;

            double neighborX = x + offsetX;
            double neighborY = y + offsetY;
            double neighborZ = z + offsetZ;

            if (neighborX >= 0 && neighborX < dim[0] &&
                neighborY >= 0 && neighborY < dim[1] &&
                neighborZ >= 0 && neighborZ < dim[2])
            {
                int neighborIndex = (neighborX) + (neighborY)*dim[0] + (neighborZ)*dim[0]*dim[1];
                effNeighCnt++;

                m[0][0] += offsetX*offsetX;  m[0][1] += offsetX*offsetY;  m[0][2] += offsetX*offsetZ;
                m[1][0] += offsetY*offsetX;  m[1][1] += offsetY*offsetY;  m[1][2] += offsetY*offsetZ;
                m[2][0] += offsetZ*offsetX;  m[2][1] += offsetZ*offsetY;  m[2][2] += offsetZ*offsetZ;

                // Relative function values of neighbors
                double gradientX = gradient[neighborIndex*3+0] - gX;
                double gradientY = gradient[neighborIndex*3+1] - gY;
                double gradientZ = gradient[neighborIndex*3+2] - gZ;

                rhsX[0] += gradientX*offsetX;  rhsX[1] += gradientX*offsetY;  rhsX[2] += gradientX*offsetZ;
                rhsY[0] += gradientY*offsetY;  rhsY[1] += gradientY*offsetY;  rhsY[2] += gradientY*offsetZ;
                rhsZ[0] += gradientZ*offsetX;  rhsZ[1] += gradientZ*offsetY;  rhsZ[2] += gradientZ*offsetZ;
            }
        }

        if (effNeighCnt < 3) { // TODO: how many neighbors needed?
            //if (effNeighCnt < 4) { // TODO: how many neighbors needed? think 4 but still can get singular ...
            //if (gradDefault) out->setVector3(i, defaultG);
        }

        double det = mat3det(m);
        if (det == 0) det = 1;
        mat3 h;

        h[0][0] = vec3det(rhsX, m[1], m[2]) / det;	// dV0/dx
        h[0][1] = vec3det(m[0], rhsX, m[2]) / det;	// dV0/dy
        h[0][2] = vec3det(m[0], m[1], rhsX) / det;	// dV0/dz

        h[1][0] = vec3det(rhsY, m[1], m[2]) / det;	// dV0/dx
        h[1][1] = vec3det(m[0], rhsY, m[2]) / det;	// dV0/dy
        h[1][2] = vec3det(m[0], m[1], rhsY) / det;	// dV0/dz

        h[2][0] = vec3det(rhsZ, m[1], m[2]) / det;	// dV0/dx
        h[2][1] = vec3det(m[0], rhsZ, m[2]) / det;	// dV0/dy
        h[2][2] = vec3det(m[0], m[1], rhsZ) / det;	// dV0/dz

//        for (int v=0; v<9; v++) {
//            if (grad[v] > FLT_MAX) {
//                grad[v] = FLT_MAX;
//            }
//            else if (grad[v] < -FLT_MAX) {
//                grad[v] = -FLT_MAX;
//            }
//        }

        hessian[index*9+0]=h[0][0];
        hessian[index*9+1]=h[0][1];
        hessian[index*9+2]=h[0][2];
        hessian[index*9+3]=h[1][0];
        hessian[index*9+4]=h[1][1];
        hessian[index*9+5]=h[1][2];
        hessian[index*9+6]=h[2][0];
        hessian[index*9+7]=h[2][1];
        hessian[index*9+8]=h[2][2];
    }
}

//__global__ void processPointGradientsOld(uint maxIndex, double* data, double* gradient, uint dim[3]) {
//    uint index=blockIdx.x*blockDim.x + threadIdx.x;
//    if (index < maxIndex) {
//    uint k=index/(dim[0]*dim[1]);
//    uint j=index%(dim[0]*dim[1])/dim[0];
//    uint i=(index%(dim[0]*dim[1]))%dim[0];
//    if (i>0 && i<dim[0]-1) {         //central differences for inner cubes
//        uint neg = index-1;
//        uint pos = index+1;
//        gradient[3*index+0]=(data[pos]-data[neg])/2.0;
//    }
//    else {
//        if (i == 0) {                //forward differences for bounding cubes
//            uint neg = index;
//            uint pos = index+1;
//            gradient[3*index+0]=(data[pos]-data[neg]);
//        }
//        if (i == dim[0]-1) {         //forward differences for bounding cubes
//            uint neg = index-1;
//            uint pos = index;
//            gradient[3*index+0]=(data[pos]-data[neg]);
//        }
//    }
//    if (j>0 && j<dim[1]-1) {
//        uint neg = index-dim[0];
//        uint pos = index+dim[0];
//        gradient[3*index+1]=(data[pos]-data[neg])/2.0;
//    }
//    else {
//        if (j == 0) {
//            uint neg = index;
//            uint pos = index+dim[0];
//            gradient[3*index+1]=(data[pos]-data[neg]);
//        }
//        if (j == dim[1]-1) {
//            uint neg = index-dim[0];
//            uint pos = index;
//            gradient[3*index+1]=(data[pos]-data[neg]);
//        }
//    }
//    if (k>0 && k<dim[2]-1) {
//        uint neg = index-dim[0]*dim[1];
//        uint pos = index+dim[0]*dim[1];
//        gradient[3*index+2]=(data[pos]-data[neg])/2.0;
//    }
//    else {
//        if (k == 0) {
//            uint neg = index;
//            uint pos = index+dim[0]*dim[1];
//            gradient[3*index+2]=(data[pos]-data[neg]);
//        }
//        if (k == dim[2]-1) {
//            uint neg = index-dim[0]*dim[1];
//            uint pos = index;
//            gradient[3*index+2]=(data[pos]-data[neg]);
//        }
//    }
//    }
//}

__global__ void processPointGradients(uint maxIndex, double* data, double* gradient, uint dim[3], unsigned int stencilRange, bool valley)
{
    uint index=blockIdx.x*blockDim.x + threadIdx.x;
    if (index < maxIndex) {
        uint z=index/(dim[0]*dim[1]);
        uint y=index%(dim[0]*dim[1])/dim[0];
        uint x=(index%(dim[0]*dim[1]))%dim[0];

        mat3 m = { {0}, {0}, {0} };
        vec3 rhs = { 0 };

        //vec3 xyzi = {x, y, z}; // ###, vi;
        double fi = data[index];

        // compute all neighbors inside level 'range'
        int width = 2*stencilRange+1;
        int neighCnt = width*width*width;//computeNodeNeighborsN(i, range, neighborsN);

        int effNeighCnt = 0;
        for (int k = 0; k < neighCnt; k++) {
            int offsetX = (k%width) - stencilRange;
            int offsetY = ((k/width)%width) - stencilRange;
            int offsetZ = (k/(width*width)) - stencilRange;

            double neighborX = x + offsetX;
            double neighborY = y + offsetY;
            double neighborZ = z + offsetZ;

            if (neighborX >= 0 && neighborX < dim[0] &&
                neighborY >= 0 && neighborY < dim[1] &&
                neighborZ >= 0 && neighborZ < dim[2])
            {
                int neighborIndex = (neighborX) + (neighborY)*dim[0] + (neighborZ)*dim[0]*dim[1];
                effNeighCnt++;

                m[0][0] += offsetX*offsetX;  m[0][1] += offsetX*offsetY;  m[0][2] += offsetX*offsetZ;
                m[1][0] += offsetY*offsetX;  m[1][1] += offsetY*offsetY;  m[1][2] += offsetY*offsetZ;
                m[2][0] += offsetZ*offsetX;  m[2][1] += offsetZ*offsetY;  m[2][2] += offsetZ*offsetZ;

                // Relative function values of neighbors
                double f = data[neighborIndex] - fi;

                rhs[0] += f*offsetX;  rhs[1] += f*offsetY;  rhs[2] += f*offsetZ;
            }
        }

        if (effNeighCnt < 3) { // TODO: how many neighbors needed?
            //if (effNeighCnt < 4) { // TODO: how many neighbors needed? think 4 but still can get singular ...
            //if (gradDefault) out->setVector3(i, defaultG);
        }

        double det = mat3det(m);
        if (det == 0) det = 1;
        vec3 grad;

        grad[0] = vec3det(rhs, m[1], m[2]) / det;	// dV0/dx
        grad[1] = vec3det(m[0], rhs, m[2]) / det;	// dV0/dy
        grad[2] = vec3det(m[0], m[1], rhs) / det;	// dV0/dz

        for (int v=0; v<3; v++) {
            if (grad[v] > FLT_MAX) {
                grad[v] = FLT_MAX;
            }
            else if (grad[v] < -FLT_MAX) {
                grad[v] = -FLT_MAX;
            }
        }

        if (valley){
            gradient[3*index+0]=-grad[0];
            gradient[3*index+1]=-grad[1];
            gradient[3*index+2]=-grad[2];
        }
        else {
            gradient[3*index+0]=grad[0];
            gradient[3*index+1]=grad[1];
            gradient[3*index+2]=grad[2];
        }
    }
}

__global__ void processIsoCubes(int maxIndex, double c, uint* iA, uint* pfA, double* d, double* g, double* vA, double* gA, uint* nVT, uint* tT, uint* bitmasks, int dim[3]) {
    uint index=blockIdx.x*blockDim.x + threadIdx.x;
    if (index < maxIndex) {
	uint curCubeIdx = iA[index];
	uint vtxArrayOffset = pfA[curCubeIdx]*3;
	uint bitmask = 0;
	bitmask |= ((d[curCubeIdx]                         < c) << 0);               //lower left front
	bitmask |= ((d[curCubeIdx+1]                       < c) << 1);               //lower right front
	bitmask |= ((d[curCubeIdx+1+dim[0]]                < c) << 2);               //upper right front
	bitmask |= ((d[curCubeIdx+dim[0]]                  < c) << 3);               //upper left front
	bitmask |= ((d[curCubeIdx+dim[0]*dim[1]]           < c) << 4);               //lower left back
	bitmask |= ((d[curCubeIdx+1+dim[0]*dim[1]]         < c) << 5);               //lower right back
	bitmask |= ((d[curCubeIdx+1+dim[0]+dim[0]*dim[1] ] < c) << 6);               //upper right back
	bitmask |= ((d[curCubeIdx+dim[0]+dim[0]*dim[1]]    < c) << 7);               //upper left back

	double i,j,k;
	k=double(curCubeIdx/(dim[0]*dim[1]));
	j=double(curCubeIdx%(dim[0]*dim[1])/dim[0]);
	i=double((curCubeIdx%(dim[0]*dim[1]))%dim[0]);

	__shared__ uint eI[24];
	eI[0*2+0] = 0;
	eI[0*2+1] = 1;
	eI[1*2+0] = 1;
	eI[1*2+1] = 1+dim[0];
	eI[2*2+0] = dim[0];
	eI[2*2+1] = 1+dim[0];
	eI[3*2+0] = 0;
	eI[3*2+1] = dim[0];
	eI[4*2+0] = dim[0]*dim[1];
	eI[4*2+1] = 1+dim[0]*dim[1];
	eI[5*2+0] = 1+dim[0]*dim[1];
	eI[5*2+1] = 1+dim[0]+dim[0]*dim[1];
	eI[6*2+0] = dim[0]+dim[0]*dim[1];
	eI[6*2+1] = 1+dim[0]+dim[0]*dim[1];
	eI[7*2+0] = dim[0]*dim[1];
	eI[7*2+1] = dim[0]+dim[0]*dim[1];
	eI[8*2+0] = 0;
	eI[8*2+1] = dim[0]*dim[1];
	eI[9*2+0] = 1;
	eI[9*2+1] = 1+dim[0]*dim[1];
	eI[10*2+0] = 1+dim[0];
	eI[10*2+1] = 1+dim[0]+dim[0]*dim[1];
	eI[11*2+0] = dim[0];
	eI[11*2+1] = dim[0]+dim[0]*dim[1];

	double eV[36];            //vertices on edges to be interpolated
	double eG[36];            //gradients on edges to be interpolated
	int pIE[12];              //pointer to interpolated edge coordinate
	eV[0*3+0]=double(0.5); pIE[0]=0*3+0;
	eV[0*3+1]=double(0);
	eV[0*3+2]=double(0);
	eV[1*3+0]=double(1);
	eV[1*3+1]=double(0.5); pIE[1]=1*3+1;
	eV[1*3+2]=double(0);
	eV[2*3+0]=double(0.5); pIE[2]=2*3+0;
	eV[2*3+1]=double(1);
	eV[2*3+2]=double(0); 
	eV[3*3+0]=double(0);
	eV[3*3+1]=double(0.5); pIE[3]=3*3+1;
	eV[3*3+2]=double(0);
	eV[4*3+0]=double(0.5); pIE[4]=4*3+0;
	eV[4*3+1]=double(0);
	eV[4*3+2]=double(1);
	eV[5*3+0]=double(1);
	eV[5*3+1]=double(0.5); pIE[5]=5*3+1;
	eV[5*3+2]=double(1);
	eV[6*3+0]=double(0.5); pIE[6]=6*3+0;
	eV[6*3+1]=double(1);
	eV[6*3+2]=double(1);
	eV[7*3+0]=double(0);
	eV[7*3+1]=double(0.5); pIE[7]=7*3+1;
	eV[7*3+2]=double(1);
	eV[8*3+0]=double(0);
	eV[8*3+1]=double(0);
	eV[8*3+2]=double(0.5); pIE[8]=8*3+2;
	eV[9*3+0]=double(1);
	eV[9*3+1]=double(0);
	eV[9*3+2]=double(0.5); pIE[9]=9*3+2;
	eV[10*3+0]=double(1);
	eV[10*3+1]=double(1);
	eV[10*3+2]=double(0.5); pIE[10]=10*3+2;
	eV[11*3+0]=double(0);
	eV[11*3+1]=double(1);
	eV[11*3+2]=double(0.5); pIE[11]=11*3+2;

	for (uint edgeIdx=0; edgeIdx<12; ++edgeIdx) {
	    double f1 = d[eI[edgeIdx*2+0]+curCubeIdx];
	    double f2 = d[eI[edgeIdx*2+1]+curCubeIdx];
	    double* g1 = &g[(eI[edgeIdx*2+0]+curCubeIdx)*3];
	    double* g2 = &g[(eI[edgeIdx*2+1]+curCubeIdx)*3];
	    double x1 = 0.0;
	    double x2 = 1.0;
	    double grad[3];
	    double x_c = findValueOE(f1,x1,f2,x2,c);
	    interpolateVector(g1,g2,x_c,grad);
	    eV[pIE[edgeIdx]] = x_c;
	    eG[edgeIdx*3+0] = grad[0];
	    eG[edgeIdx*3+1] = grad[1];
	    eG[edgeIdx*3+2] = grad[2];
	}
	
	for (uint triNr=0; triNr<nVT[bitmask]/3; ++triNr) {
	    uint edge0 = tT[bitmask*16+triNr*3+0];
	    uint edge1 = tT[bitmask*16+triNr*3+1];
	    uint edge2 = tT[bitmask*16+triNr*3+2];
	    
	    vA[vtxArrayOffset+9*triNr+0] = eV[edge0*3+0]+i;
	    vA[vtxArrayOffset+9*triNr+1] = eV[edge0*3+1]+j;
	    vA[vtxArrayOffset+9*triNr+2] = eV[edge0*3+2]+k;

	    gA[vtxArrayOffset+9*triNr+0] = eG[edge0*3+0];
	    gA[vtxArrayOffset+9*triNr+1] = eG[edge0*3+1];
	    gA[vtxArrayOffset+9*triNr+2] = eG[edge0*3+2];

	    vA[vtxArrayOffset+9*triNr+3] = eV[edge1*3+0]+i;
	    vA[vtxArrayOffset+9*triNr+4] = eV[edge1*3+1]+j;
	    vA[vtxArrayOffset+9*triNr+5] = eV[edge1*3+2]+k;

	    gA[vtxArrayOffset+9*triNr+3] = eG[edge1*3+0];
	    gA[vtxArrayOffset+9*triNr+4] = eG[edge1*3+1];
	    gA[vtxArrayOffset+9*triNr+5] = eG[edge1*3+2];

	    vA[vtxArrayOffset+9*triNr+6] = eV[edge2*3+0]+i;
	    vA[vtxArrayOffset+9*triNr+7] = eV[edge2*3+1]+j;
	    vA[vtxArrayOffset+9*triNr+8] = eV[edge2*3+2]+k;

	    gA[vtxArrayOffset+9*triNr+6] = eG[edge2*3+0];
	    gA[vtxArrayOffset+9*triNr+7] = eG[edge2*3+1];
	    gA[vtxArrayOffset+9*triNr+8] = eG[edge2*3+2];
	}
    }
}


__global__ void processRidgeCubes(int maxIndex, uint dim[3], uint* iA, uint* pfA, double* ev, double* g, double* lEV, double* vA, double* gA, uint* nVT, uint* tT, uint* bitmasks) {
    uint index=blockIdx.x*blockDim.x + threadIdx.x;

    __shared__ uint ofst[8];
    ofst[0] = 0;                       
    ofst[1] = 1;                       
    ofst[2] = 1+dim[0];                  
    ofst[3] = dim[0];                
    ofst[4] = 0+dim[0]*dim[1];
    ofst[5] = 1+dim[0]*dim[1];
    ofst[6] = 1+dim[0]+dim[0]*dim[1];
    ofst[7] = dim[0]+dim[0]*dim[1];

    if (index < maxIndex) {
	const uint numNodes = 8;
	uint curCubeIdx = iA[index];
	uint vtxArrayOffset = pfA[curCubeIdx]*3;
	uint bitmask = 0;
	uint k_int=curCubeIdx/(dim[0]*dim[1]);
	uint j_int=curCubeIdx%(dim[0]*dim[1])/dim[0];
	uint i_int=(curCubeIdx%(dim[0]*dim[1]))%dim[0];
	double gradDotEV[numNodes];
	double grads[numNodes*3];
	if(i_int !=  dim[0]-1 && j_int != dim[1]-1 && k_int != dim[2]-1) {
	    double orientedEV[numNodes*3];
	    double grads[numNodes*3];
	    for (int n=0; n<numNodes; ++n) {
		orientedEV[n*3+0] = ev[(curCubeIdx+ofst[n])*3+0];
		orientedEV[n*3+1] = ev[(curCubeIdx+ofst[n])*3+1];
		orientedEV[n*3+2] = ev[(curCubeIdx+ofst[n])*3+2];
		
		grads[n*3+0] = g[(curCubeIdx+ofst[n])*3+0];
		grads[n*3+1] = g[(curCubeIdx+ofst[n])*3+1];
		grads[n*3+2] = g[(curCubeIdx+ofst[n])*3+2];

	        if (dotProduct(&lEV[curCubeIdx*3], &orientedEV[n*3]) < 0) {
		    orientedEV[n*3+0] *= -1;
		    orientedEV[n*3+1] *= -1;
		    orientedEV[n*3+2] *= -1;
		}
		gradDotEV[n] = dotProduct(&grads[n*3], &orientedEV[n*3]);
	    }
	    bitmask |= ((gradDotEV[0] < 0) << 0);               //lower left front
	    bitmask |= ((gradDotEV[1] < 0) << 1);               //lower right front
	    bitmask |= ((gradDotEV[2] < 0) << 2);               //upper right front
	    bitmask |= ((gradDotEV[3] < 0) << 3);               //upper left front
	    bitmask |= ((gradDotEV[4] < 0) << 4);               //lower left back
	    bitmask |= ((gradDotEV[5] < 0) << 5);               //lower right back
	    bitmask |= ((gradDotEV[6] < 0) << 6);               //upper right back
	    bitmask |= ((gradDotEV[7] < 0) << 7);               //upper left back
	}
	double i,j,k;
	k=double(k_int);
	j=double(j_int);
	i=double(i_int);

	__shared__ uint eI[24];
	eI[0*2+0] = 0;
	eI[0*2+1] = 1;
	eI[2*2+0] = 3;
	eI[2*2+1] = 2;
	eI[4*2+0] = 4;
	eI[4*2+1] = 5;
	eI[6*2+0] = 7;
	eI[6*2+1] = 6;

	eI[1*2+0] = 1;
	eI[1*2+1] = 2;
	eI[3*2+0] = 0;
	eI[3*2+1] = 3;
	eI[5*2+0] = 5;
	eI[5*2+1] = 6;
	eI[7*2+0] = 4;
	eI[7*2+1] = 7;

	eI[8*2+0] = 0;
	eI[8*2+1] = 4;
	eI[9*2+0] = 1;
	eI[9*2+1] = 5;
	eI[10*2+0] = 2;
	eI[10*2+1] = 6;
	eI[11*2+0] = 3;
	eI[11*2+1] = 7;

	double eV[36];            //vertices on edges to be interpolated
	double eG[36];            //gradients on edges to be interpolated
	int pIE[12];              //pointer to interpolated edge coordinate
	eV[0*3+0]=double(0.5); pIE[0]=0*3+0;
	eV[0*3+1]=double(0);
	eV[0*3+2]=double(0);
	eV[1*3+0]=double(1);
	eV[1*3+1]=double(0.5); pIE[1]=1*3+1;
	eV[1*3+2]=double(0);
	eV[2*3+0]=double(0.5); pIE[2]=2*3+0;
	eV[2*3+1]=double(1);
	eV[2*3+2]=double(0); 
	eV[3*3+0]=double(0);
	eV[3*3+1]=double(0.5); pIE[3]=3*3+1;
	eV[3*3+2]=double(0);
	eV[4*3+0]=double(0.5); pIE[4]=4*3+0;
	eV[4*3+1]=double(0);
	eV[4*3+2]=double(1);
	eV[5*3+0]=double(1);
	eV[5*3+1]=double(0.5); pIE[5]=5*3+1;
	eV[5*3+2]=double(1);
	eV[6*3+0]=double(0.5); pIE[6]=6*3+0;
	eV[6*3+1]=double(1);
	eV[6*3+2]=double(1);
	eV[7*3+0]=double(0);
	eV[7*3+1]=double(0.5); pIE[7]=7*3+1;
	eV[7*3+2]=double(1);
	eV[8*3+0]=double(0);
	eV[8*3+1]=double(0);
	eV[8*3+2]=double(0.5); pIE[8]=8*3+2;
	eV[9*3+0]=double(1);
	eV[9*3+1]=double(0);
	eV[9*3+2]=double(0.5); pIE[9]=9*3+2;
	eV[10*3+0]=double(1);
	eV[10*3+1]=double(1);
	eV[10*3+2]=double(0.5); pIE[10]=10*3+2;
	eV[11*3+0]=double(0);
	eV[11*3+1]=double(1);
	eV[11*3+2]=double(0.5); pIE[11]=11*3+2;

	for (uint edgeIdx=0; edgeIdx<12; ++edgeIdx) {
	    double f1 = gradDotEV[eI[edgeIdx*2+0]];
	    double f2 = gradDotEV[eI[edgeIdx*2+1]];
	    double* g1 = &grads[(eI[edgeIdx*2+0])*3];
	    double* g2 = &grads[(eI[edgeIdx*2+1])*3];
	    double x1 = 0.0;
	    double x2 = 1.0;
	    double igrad[3];
	    double x_c = findValueOE(f1,x1,f2,x2,0.0);	
	    interpolateVector(g1,g2,x_c,igrad);
	    // interpolateVector(g1,g2,x_c,grad);
	    // if (edgeIdx ==  0 || edgeIdx ==  2 || edgeIdx ==  4 || edgeIdx ==  6 ||
	    // 	edgeIdx ==  1 || edgeIdx ==  3 || edgeIdx ==  5 || edgeIdx ==  7 ||
	    // 	edgeIdx ==  8 || edgeIdx ==  9 || edgeIdx == 10 || edgeIdx == 11 )
	    	eV[pIE[edgeIdx]] = x_c;
	    // else
	    // 	eV[pIE[edgeIdx]] = 0.5;

	    eG[edgeIdx*3+0] = igrad[0];
	    eG[edgeIdx*3+1] = igrad[1];
	    eG[edgeIdx*3+2] = igrad[2];
	}
	for (uint triNr=0; triNr<nVT[bitmask]/3; ++triNr) {
	    uint edge0 = tT[bitmask*16+triNr*3+0];
	    uint edge1 = tT[bitmask*16+triNr*3+1];
	    uint edge2 = tT[bitmask*16+triNr*3+2];
	    
	    vA[vtxArrayOffset+9*triNr+0] = eV[edge0*3+0]+i;
	    vA[vtxArrayOffset+9*triNr+1] = eV[edge0*3+1]+j;
	    vA[vtxArrayOffset+9*triNr+2] = eV[edge0*3+2]+k;

	    gA[vtxArrayOffset+9*triNr+0] = eG[edge0*3+0];
	    gA[vtxArrayOffset+9*triNr+1] = eG[edge0*3+1];
	    gA[vtxArrayOffset+9*triNr+2] = eG[edge0*3+2];

	    vA[vtxArrayOffset+9*triNr+3] = eV[edge1*3+0]+i;
	    vA[vtxArrayOffset+9*triNr+4] = eV[edge1*3+1]+j;
	    vA[vtxArrayOffset+9*triNr+5] = eV[edge1*3+2]+k;

	    gA[vtxArrayOffset+9*triNr+3] = eG[edge1*3+0];
	    gA[vtxArrayOffset+9*triNr+4] = eG[edge1*3+1];
	    gA[vtxArrayOffset+9*triNr+5] = eG[edge1*3+2];

	    vA[vtxArrayOffset+9*triNr+6] = eV[edge2*3+0]+i;
	    vA[vtxArrayOffset+9*triNr+7] = eV[edge2*3+1]+j;
	    vA[vtxArrayOffset+9*triNr+8] = eV[edge2*3+2]+k;

	    gA[vtxArrayOffset+9*triNr+6] = eG[edge2*3+0];
	    gA[vtxArrayOffset+9*triNr+7] = eG[edge2*3+1];
	    gA[vtxArrayOffset+9*triNr+8] = eG[edge2*3+2];
	}
    }
}

#define IDX 12345
__global__ void processRidgeCubesData(int maxIndex,
                                      uint dim[3],
                                      uint* iA,
                                      uint* pfA,
                                      double* ev,
                                      double* g,
                                      double* lEV,
                                      double* vA,
                                      double* gA,
                                      uint* nVT,
                                      uint* tT,
                                      uint* bitmasks,
                                      double *data,
                                      double *pointData,
                                      double3 origin,
                                      double3 spacing) {

    uint index=blockIdx.x*blockDim.x + threadIdx.x;

    __shared__ uint ofst[8];
    ofst[0] = 0;                       
    ofst[1] = 1;                       
    ofst[2] = 1+dim[0];                  
    ofst[3] = dim[0];                
    ofst[4] = 0+dim[0]*dim[1];
    ofst[5] = 1+dim[0]*dim[1];
    ofst[6] = 1+dim[0]+dim[0]*dim[1];
    ofst[7] = dim[0]+dim[0]*dim[1];

    if (index < maxIndex) {
	const uint numNodes = 8;
    uint curCubeIdx = iA[index];
    uint vtxArrayOffset = pfA[curCubeIdx]*3;
    uint pointDataArrayOffset = pfA[curCubeIdx];
	uint bitmask = 0;
	uint k_int=curCubeIdx/(dim[0]*dim[1]);
	uint j_int=curCubeIdx%(dim[0]*dim[1])/dim[0];
	uint i_int=(curCubeIdx%(dim[0]*dim[1]))%dim[0];
	double gradDotEV[numNodes];
    double grads[numNodes*3];
    double nodeData[numNodes];
    
	if(i_int !=  dim[0]-1 && j_int != dim[1]-1 && k_int != dim[2]-1) {
	    double orientedEV[numNodes*3];
        //double grads[numNodes*3];
	    for (int n=0; n<numNodes; ++n) {
			orientedEV[n*3+0] = ev[(curCubeIdx+ofst[n])*3+0];
			orientedEV[n*3+1] = ev[(curCubeIdx+ofst[n])*3+1];
			orientedEV[n*3+2] = ev[(curCubeIdx+ofst[n])*3+2];
			
            grads[n*3+0] = g[(curCubeIdx+ofst[n])*3+0];
            grads[n*3+1] = g[(curCubeIdx+ofst[n])*3+1];
            grads[n*3+2] = g[(curCubeIdx+ofst[n])*3+2];

			if (dotProduct(&lEV[curCubeIdx*3], &orientedEV[n*3]) < 0) {
				orientedEV[n*3+0] *= -1;
				orientedEV[n*3+1] *= -1;
				orientedEV[n*3+2] *= -1;
			}
			gradDotEV[n] = dotProduct(&grads[n*3], &orientedEV[n*3]);
            nodeData[n] = data[curCubeIdx+ofst[n]];
	    }
	    bitmask |= ((gradDotEV[0] < 0) << 0);               //lower left front
	    bitmask |= ((gradDotEV[1] < 0) << 1);               //lower right front
	    bitmask |= ((gradDotEV[2] < 0) << 2);               //upper right front
	    bitmask |= ((gradDotEV[3] < 0) << 3);               //upper left front
	    bitmask |= ((gradDotEV[4] < 0) << 4);               //lower left back
	    bitmask |= ((gradDotEV[5] < 0) << 5);               //lower right back
	    bitmask |= ((gradDotEV[6] < 0) << 6);               //upper right back
	    bitmask |= ((gradDotEV[7] < 0) << 7);               //upper left back
	}
	double i,j,k;
	k=double(k_int);
	j=double(j_int);
	i=double(i_int);
	
	
	//if (index%1000 == 3) printf("%f %f %f %f   %f %f %f %f\n",
	//							nodeData[0], nodeData[1], nodeData[2], nodeData[3],
	//							nodeData[4], nodeData[5], nodeData[6], nodeData[7]);

	__shared__ uint eI[24];
	eI[0*2+0] = 0;
	eI[0*2+1] = 1;
	eI[2*2+0] = 3;
	eI[2*2+1] = 2;
	eI[4*2+0] = 4;
	eI[4*2+1] = 5;
	eI[6*2+0] = 7;
	eI[6*2+1] = 6;

	eI[1*2+0] = 1;
	eI[1*2+1] = 2;
	eI[3*2+0] = 0;
	eI[3*2+1] = 3;
	eI[5*2+0] = 5;
	eI[5*2+1] = 6;
	eI[7*2+0] = 4;
	eI[7*2+1] = 7;

	eI[8*2+0] = 0;
	eI[8*2+1] = 4;
	eI[9*2+0] = 1;
	eI[9*2+1] = 5;
	eI[10*2+0] = 2;
	eI[10*2+1] = 6;
	eI[11*2+0] = 3;
	eI[11*2+1] = 7;

	double eV[36];            //vertices on edges to be interpolated
//    double eG[36];            //gradients on edges to be interpolated
	int pIE[12];              //pointer to interpolated edge coordinate
	eV[0*3+0]=double(0.5); pIE[0]=0*3+0;
	eV[0*3+1]=double(0);
	eV[0*3+2]=double(0);
	eV[1*3+0]=double(1);
	eV[1*3+1]=double(0.5); pIE[1]=1*3+1;
	eV[1*3+2]=double(0);
	eV[2*3+0]=double(0.5); pIE[2]=2*3+0;
	eV[2*3+1]=double(1);
	eV[2*3+2]=double(0); 
	eV[3*3+0]=double(0);
	eV[3*3+1]=double(0.5); pIE[3]=3*3+1;
	eV[3*3+2]=double(0);
	eV[4*3+0]=double(0.5); pIE[4]=4*3+0;
	eV[4*3+1]=double(0);
	eV[4*3+2]=double(1);
	eV[5*3+0]=double(1);
	eV[5*3+1]=double(0.5); pIE[5]=5*3+1;
	eV[5*3+2]=double(1);
	eV[6*3+0]=double(0.5); pIE[6]=6*3+0;
	eV[6*3+1]=double(1);
	eV[6*3+2]=double(1);
	eV[7*3+0]=double(0);
	eV[7*3+1]=double(0.5); pIE[7]=7*3+1;
	eV[7*3+2]=double(1);
	eV[8*3+0]=double(0);
	eV[8*3+1]=double(0);
	eV[8*3+2]=double(0.5); pIE[8]=8*3+2;
	eV[9*3+0]=double(1);
	eV[9*3+1]=double(0);
	eV[9*3+2]=double(0.5); pIE[9]=9*3+2;
	eV[10*3+0]=double(1);
	eV[10*3+1]=double(1);
	eV[10*3+2]=double(0.5); pIE[10]=10*3+2;
	eV[11*3+0]=double(0);
	eV[11*3+1]=double(1);
	eV[11*3+2]=double(0.5); pIE[11]=11*3+2;

	for (uint edgeIdx=0; edgeIdx<12; ++edgeIdx) {
        double f1 = gradDotEV[eI[edgeIdx*2+0]];
        double f2 = gradDotEV[eI[edgeIdx*2+1]];
//        double* g1 = &grads[(eI[edgeIdx*2+0])*3];
//        double* g2 = &grads[(eI[edgeIdx*2+1])*3];
	    double x1 = 0.0;
	    double x2 = 1.0;
//        double igrad[3];
	    double x_c = findValueOE(f1,x1,f2,x2,0.0);	
//        interpolateVector(g1,g2,x_c,igrad);
	    // interpolateVector(g1,g2,x_c,grad);
	    // if (edgeIdx ==  0 || edgeIdx ==  2 || edgeIdx ==  4 || edgeIdx ==  6 ||
	    // 	edgeIdx ==  1 || edgeIdx ==  3 || edgeIdx ==  5 || edgeIdx ==  7 ||
	    // 	edgeIdx ==  8 || edgeIdx ==  9 || edgeIdx == 10 || edgeIdx == 11 )
	    	eV[pIE[edgeIdx]] = x_c;
	    // else
	    // 	eV[pIE[edgeIdx]] = 0.5;

//        eG[edgeIdx*3+0] = igrad[0];
//        eG[edgeIdx*3+1] = igrad[1];
//        eG[edgeIdx*3+2] = igrad[2];
	}
	    
	for (uint triNr=0; triNr<nVT[bitmask]/3; ++triNr) {
        uint edge0 = tT[bitmask*16+triNr*3+0];
        uint edge1 = tT[bitmask*16+triNr*3+1];
        uint edge2 = tT[bitmask*16+triNr*3+2];
	    
        vA[vtxArrayOffset+9*triNr+0] = spacing.x*(eV[edge0*3+0]+i) + origin.x;
        vA[vtxArrayOffset+9*triNr+1] = spacing.y*(eV[edge0*3+1]+j) + origin.y;
        vA[vtxArrayOffset+9*triNr+2] = spacing.z*(eV[edge0*3+2]+k) + origin.z;

//        gA[vtxArrayOffset+9*triNr+0] = eG[edge0*3+0];
//        gA[vtxArrayOffset+9*triNr+1] = eG[edge0*3+1];
//        gA[vtxArrayOffset+9*triNr+2] = eG[edge0*3+2];
	    
        double c0 = eV[pIE[edge0]];
        pointData[pointDataArrayOffset + triNr*3] = nodeData[eI[edge0*2]] * (1 - c0) + nodeData[eI[edge0*2+1]] * c0;

        vA[vtxArrayOffset+9*triNr+3] = spacing.x*(eV[edge1*3+0]+i) + origin.x;
        vA[vtxArrayOffset+9*triNr+4] = spacing.y*(eV[edge1*3+1]+j) + origin.y;
        vA[vtxArrayOffset+9*triNr+5] = spacing.z*(eV[edge1*3+2]+k) + origin.z;

//        gA[vtxArrayOffset+9*triNr+3] = eG[edge1*3+0];
//        gA[vtxArrayOffset+9*triNr+4] = eG[edge1*3+1];
//        gA[vtxArrayOffset+9*triNr+5] = eG[edge1*3+2];
	    
        double c1 = eV[pIE[edge1]];
        pointData[pointDataArrayOffset + triNr*3 + 1] = nodeData[eI[edge1*2]] * (1 - c1) + nodeData[eI[edge1*2+1]] * c1;

        vA[vtxArrayOffset+9*triNr+6] = spacing.x*(eV[edge2*3+0]+i) + origin.x;
        vA[vtxArrayOffset+9*triNr+7] = spacing.y*(eV[edge2*3+1]+j) + origin.y;
        vA[vtxArrayOffset+9*triNr+8] = spacing.z*(eV[edge2*3+2]+k) + origin.z;

//        gA[vtxArrayOffset+9*triNr+6] = eG[edge2*3+0];
//        gA[vtxArrayOffset+9*triNr+7] = eG[edge2*3+1];
//        gA[vtxArrayOffset+9*triNr+8] = eG[edge2*3+2];
	    
        double c2 = eV[pIE[edge2]];
        pointData[pointDataArrayOffset + triNr*3 + 2] = nodeData[eI[edge2*2]] * (1 - c2) + nodeData[eI[edge2*2+1]] * c2;
	}
    }
}

__global__ void getIsoIndices(int maxIndex, uint* iA, uint* bitmasks, uint* pfA) {
    uint index=blockIdx.x*blockDim.x + threadIdx.x;
    if (index < maxIndex) {
	if (bitmasks[index+maxIndex]) {
	    iA[pfA[index+maxIndex]] = index;
	}
    }
}

__global__ void processDataArray(int maxIndex, double c, int dim[3], double* d, uint* bitmasks, uint* nVT){
    int index=blockIdx.x*blockDim.x + threadIdx.x;
    if (index < maxIndex) {
	uint bitmask = 0;
	uint i,j,k;
	k=index/(dim[0]*dim[1]);
	j=index%(dim[0]*dim[1])/dim[0];
	i=(index%(dim[0]*dim[1]))%dim[0];
	if(i !=  dim[0]-1 && j != dim[1]-1 && k != dim[2]-1) {
	    bitmask |= ((d[index]                         < c) << 0);               //lower left front
	    bitmask |= ((d[index+1]                       < c) << 1);               //lower right front
	    bitmask |= ((d[index+1+dim[0]]                < c) << 2);               //upper right front
	    bitmask |= ((d[index+dim[0]]                  < c) << 3);               //upper left front
	    bitmask |= ((d[index+dim[0]*dim[1]]           < c) << 4);               //lower left back
	    bitmask |= ((d[index+1+dim[0]*dim[1]]         < c) << 5);               //lower right back
	    bitmask |= ((d[index+1+dim[0]+dim[0]*dim[1] ] < c) << 6);               //upper right back
	    bitmask |= ((d[index+dim[0]+dim[0]*dim[1]]    < c) << 7);               //upper left back
	    bitmasks[index] = nVT[bitmask];
	    bitmasks[index+maxIndex] = ((bitmask) && (~bitmask & 0xFF));
	}
	else {
	    bitmasks[index] = 0;
	    bitmasks[index+maxIndex] = 0;//((bitmask) && (~bitmask & 0xFF) & 0x1);
	}
    }
}

#endif
