#ifndef __KERNELS_CU__
#define __KERNELS_CU__

__global__ void processPoints(uint maxIndex, double* d, double* g, uint dim[3]) {
    uint index=blockIdx.x*blockDim.x + threadIdx.x;
    if (index < maxIndex) {
	uint k=index/(dim[0]*dim[1]);
	uint j=index%(dim[0]*dim[1])/dim[0];
	uint i=(index%(dim[0]*dim[1]))%dim[0];
	if (i>0 && i<dim[0]-1) {         //central differences for inner cubes
	    uint neg = index-1;
	    uint pos = index+1;
	    g[3*index+0]=(d[pos]-d[neg])/2.0;
	}
	else {
	    if (i == 0) {                //forward differences for bounding cubes
	    	uint neg = index;
	    	uint pos = index+1;
	    	g[3*index+0]=(d[pos]-d[neg]);
	    }
	    if (i == dim[0]-1) {         //forward differences for bounding cubes
	    	uint neg = index-1;
	    	uint pos = index;
	    	g[3*index+0]=(d[pos]-d[neg]);
	    }
	}
	if (j>0 && j<dim[1]-1) {
	    uint neg = index-dim[0];
	    uint pos = index+dim[0];
	    g[3*index+1]=(d[pos]-d[neg])/2.0;
	}
	else {
	    if (j == 0) {
	    	uint neg = index;
	    	uint pos = index+dim[0];
	    	g[3*index+1]=(d[pos]-d[neg]);
	    }
	    if (j == dim[1]-1) {
	    	uint neg = index-dim[0];
	    	uint pos = index;
	    	g[3*index+1]=(d[pos]-d[neg]);
	    }
	}
	if (k>0 && k<dim[2]-1) {
	    uint neg = index-dim[0]*dim[1];
	    uint pos = index+dim[0]*dim[1];
	    g[3*index+2]=(d[pos]-d[neg])/2.0;
	}
	else {
	    if (k == 0) {
	    	uint neg = index;
	    	uint pos = index+dim[0]*dim[1];
	    	g[3*index+2]=(d[pos]-d[neg]);
	    }
	    if (k == dim[2]-1) {
	    	uint neg = index-dim[0]*dim[1];
	    	uint pos = index;
	    	g[3*index+2]=(d[pos]-d[neg]);
	    }
	}
    }
}

__device__ double findValueOE(double f1, double x1, double f2, double x2, double c){
    return ((f2-c)*x1+(c-f1)*x2)/(f2-f1);
}

__device__ void findGradientOE(double f1[3], double x1, double f2[3], double x2, double c, double* out){
    out[0] = ((f2[0]-c)*x1+(c-f1[0])*x2)/(f2[0]-f1[0]);
    out[1] = ((f2[1]-c)*x1+(c-f1[1])*x2)/(f2[1]-f1[1]);
    out[2] = ((f2[2]-c)*x1+(c-f1[2])*x2)/(f2[2]-f1[2]);
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
	    double* g1 = &g[eI[edgeIdx*2+0]+curCubeIdx];
	    double* g2 = &g[eI[edgeIdx*2+1]+curCubeIdx];
	    double x1 = 0.0;
	    double x2 = 1.0;
	    double grad[3];
	    eV[pIE[edgeIdx]] = findValueOE(f1,x1,f2,x2,c);
	    findGradientOE(g1,x1,g2,x2,c,grad);
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

	    gA[vtxArrayOffset+9*triNr+3] = eG[edge0*3+0];
	    gA[vtxArrayOffset+9*triNr+4] = eG[edge0*3+1];
	    gA[vtxArrayOffset+9*triNr+5] = eG[edge0*3+2];

	    vA[vtxArrayOffset+9*triNr+6] = eV[edge2*3+0]+i;
	    vA[vtxArrayOffset+9*triNr+7] = eV[edge2*3+1]+j;
	    vA[vtxArrayOffset+9*triNr+8] = eV[edge2*3+2]+k;

	    gA[vtxArrayOffset+9*triNr+6] = eG[edge0*3+0];
	    gA[vtxArrayOffset+9*triNr+7] = eG[edge0*3+1];
	    gA[vtxArrayOffset+9*triNr+8] = eG[edge0*3+2];
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
