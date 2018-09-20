#include "isosurface.h"

#include "vtkPolyDataAlgorithm.h" //superclass
#include "vtkExecutive.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkTriangle.h"
#include "vtkFloatArray.h"

//#include "mctable.h" // <-- use extern defines here
extern uint edgeTable[256];
extern uint triTable[256][16];
extern uint numVertsTable[256];



//---Isosurface---------------------------------------------------------------
Isosurface::Isosurface(){
    m_points = vtkSmartPointer<vtkPoints>::New();    
    m_tris = vtkSmartPointer<vtkCellArray>::New();    
    m_polydata = vtkSmartPointer<vtkPolyData>::New();
    m_grads = vtkSmartPointer<vtkFloatArray>::New();
}

Isosurface::~Isosurface(){
    m_dim[0] = 0;
    m_dim[1] = 0;
    m_dim[2] = 0;
    m_data = NULL;
    m_grad_data = NULL;
    m_hasGradients = false;
}

double Isosurface::findValueOnEdge(double f1, double x1, double f2, double x2, double c){
    return ((f2-c)*x1+(c-f1)*x2)/(f2-f1);
}

void Isosurface::generateTriangles(double c) {

    // double *vertsCPU = NULL;
    double *vertsGPU = NULL;
    double *gradsGPU = NULL;
    double *verts = NULL;
    double *trigrads = NULL;

    uint numVertsGPU;

    // generateTrianglesCPU(c, &numVertsCPU, &vertsCPU);
    // grads are the gradients per point of dataset, gradsGPU will be the calculated interpolated ones for the triangles
    generateTrianglesGPU(c, &numVertsGPU, &vertsGPU, &gradsGPU);

    m_points->Reset();
    m_tris->Reset();
    
    verts = vertsGPU;
    trigrads = gradsGPU;

    for(int i=0; i<numVertsGPU; i+=3) {
	vtkIdType id[1];
	vtkSmartPointer<vtkTriangle> tri = vtkSmartPointer<vtkTriangle>::New();

	id[0] = m_points->InsertNextPoint(&verts[(i+0)*3]);
	tri->GetPointIds()->SetId(0, id[0]);

	id[0] = m_points->InsertNextPoint(&verts[(i+1)*3]);
	tri->GetPointIds()->SetId(1, id[0]);

	id[0] = m_points->InsertNextPoint(&verts[(i+2)*3]);
	tri->GetPointIds()->SetId(2, id[0]);
	m_tris->InsertNextCell(tri);
    }

    if (vertsGPU) delete [] vertsGPU; vertsGPU = NULL;

    cout << "Storing " << m_points->GetNumberOfPoints() << " points." << endl;
    m_polydata->SetPoints(m_points);
    cout << "Storing " << m_tris->GetNumberOfCells() << " triangles." << endl;
    m_polydata->SetPolys(m_tris);

    cout << "Triangles generated. Exiting method." << endl;
}

void Isosurface::generateTrianglesCPU(double c, uint* numVerts, double** verts) {
    cout << "\nIsosurface with interpolation on CPU" << endl;

    double eV[36];
    double* pIE[12];
    int edgeDataIndices[12][2];
    //timespec start, end;
    uint *bm = new uint[2*m_dim[0]*m_dim[1]*m_dim[2]];
    std::vector<double> vA;

    GetEdgeVertices(eV, pIE);

    //clock_gettime(CLOCK_MONOTONIC, &start);

    for (int k=0; k<m_dim[2]-1; ++k) {
	for (int j=0; j<m_dim[1]-1; ++j) {
	    for (int i=0; i<m_dim[0]-1; ++i) {	
		uint bitmask = 0;
		int index = i + j*m_dim[0]+ k*(m_dim[0]*m_dim[1]);
		if (m_data[index]                              < c) bitmask |= 0x01;               //lower left front
		if (m_data[index+1]                            < c) bitmask |= 0x02;               //lower right front
		if (m_data[index+1+m_dim[0]]                   < c) bitmask |= 0x04;               //upper right front
		if (m_data[index+m_dim[0]]                     < c) bitmask |= 0x08;               //upper left front
		if (m_data[index+m_dim[0]*m_dim[1]]            < c) bitmask |= 0x10;               //lower left back
		if (m_data[index+1+m_dim[0]*m_dim[1]]          < c) bitmask |= 0x20;               //lower right back
		if (m_data[index+1+m_dim[0]+m_dim[0]*m_dim[1]] < c) bitmask |= 0x40;               //upper right back
		if (m_data[index+m_dim[0]+m_dim[0]*m_dim[1]]   < c) bitmask |= 0x80;               //upper left back
		bm[index] = numVertsTable[bitmask];
		bm[index+m_dim[0]*m_dim[1]*m_dim[2]] = 0;
		uint edgebit=0x1;
		GetEdgeIndices(index, (int*)edgeDataIndices);
		for(int edgeIdx = 0; edgeIdx<12; ++edgeIdx) {
		    if (edgeTable[bitmask] & edgebit) {
			double f1 = m_data[edgeDataIndices[edgeIdx][0]];
			double x1 = 0.0;
			double f2 = m_data[edgeDataIndices[edgeIdx][1]];
			double x2 = 1.0;
			*pIE[edgeIdx] = findValueOnEdge(f1,x1,f2,x2,c);
		    }
		    edgebit = edgebit << 1;
		}
		for (int triNr=0; triNr<numVertsTable[bitmask]/3; ++triNr) {
		    vtkSmartPointer<vtkTriangle> tri = vtkSmartPointer<vtkTriangle>::New();
		    vtkIdType id[1];
		    double p0[3], p1[3], p2[3];
		    int edge0 = triTable[bitmask][0+triNr*3];
		    int edge1 = triTable[bitmask][1+triNr*3];
		    int edge2 = triTable[bitmask][2+triNr*3];
		    
		    p0[0] = eV[edge0*3+0]+double(i);
		    p0[1] = eV[edge0*3+1]+double(k);
		    p0[2] = eV[edge0*3+2]+double(j);
		    vA.push_back(p0[0]);
		    vA.push_back(p0[1]);
		    vA.push_back(p0[2]);
		    
		    p1[0] = eV[edge1*3+0]+double(i);
		    p1[1] = eV[edge1*3+1]+double(k);
		    p1[2] = eV[edge1*3+2]+double(j);
		    vA.push_back(p1[0]);
		    vA.push_back(p1[1]);
		    vA.push_back(p1[2]);
		    
		    p2[0] = eV[edge2*3+0]+double(i);
		    p2[1] = eV[edge2*3+1]+double(k);
		    p2[2] = eV[edge2*3+2]+double(j);
		    vA.push_back(p2[0]);
		    vA.push_back(p2[1]);
		    vA.push_back(p2[2]);
		}
	    }
	}
    }
    delete [] bm;

    //clock_gettime(CLOCK_MONOTONIC, &end);
    //cout << "Time used for CPU calculation in milliseconds:" << double(end.tv_nsec-start.tv_nsec)/1000000.0 << endl;

    *numVerts = vA.size()/3;
    *verts = new double[vA.size()];
    memcpy(*verts, &vA[0], vA.size()*sizeof(double));
}

void Isosurface::generateTrianglesGPU(double c, uint* numVerts, double** verts, double** triGrads) {
    cout << "\nIsosurface with interpolation on GPU" << endl;
    cudaIsosurface(c, m_data, m_grad_data, m_dim, numVerts, verts, triGrads);
}

//----Gradients---------------------------------------------------------------------------------
void Isosurface::generateGradients() {
    // double *gradsGPU = NULL;
    generateGradientsGPU();
    m_hasGradients = true;
}

void Isosurface::generateGradientsGPU() {
#ifdef ISO_SHOW_GRADIENTS
    m_grads->SetNumberOfComponents(3);
    m_grads->SetName("Gradients");
    vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
#endif

    cudaGradient(m_data, m_dim, &m_grad_data);

#ifdef ISO_SHOW_GRADIENTS
    for (uint idx=0; idx<m_dim[0]*m_dim[1]*m_dim[2]; ++idx) {
    	vtkIdType id[1];
    	id[0]=idx;
    	double k=double(idx/(m_dim[0]*m_dim[1]));
    	double j=double(idx%(m_dim[0]*m_dim[1])/m_dim[0]);
    	double i=double((idx%(m_dim[0]*m_dim[1]))%m_dim[0]);
	m_points->InsertPoint(id[0],i,j,k);
	verts->InsertNextCell(VTK_VERTEX, id);
	m_grads->InsertNextTuple(&m_grad_data[idx*3]);
    }

    // cout << "numInner:" << numInner << endl;
    m_polydata->SetPoints(m_points);
    m_polydata->SetVerts(verts);
    m_polydata->GetPointData()->AddArray(m_grads);
#endif
}

void Isosurface::GetEdgeIndices(int index, int* outArray ){
    outArray[0*2+0] = index;
    outArray[0*2+1] = index+1;

    outArray[1*2+0] = index+1;
    outArray[1*2+1] = index+1+m_dim[0];

    outArray[2*2+0] = index+m_dim[0];
    outArray[2*2+1] = index+1+m_dim[0];

    outArray[3*2+0] = index;
    outArray[3*2+1] = index+m_dim[0];

    outArray[4*2+0] = index+m_dim[0]*m_dim[1];
    outArray[4*2+1] = index+1+m_dim[0]*m_dim[1];

    outArray[5*2+0] = index+1+m_dim[0]*m_dim[1];
    outArray[5*2+1] = index+1+m_dim[0]+m_dim[0]*m_dim[1];

    outArray[6*2+0] = index+m_dim[0]+m_dim[0]*m_dim[1];
    outArray[6*2+1] = index+1+m_dim[0]+m_dim[0]*m_dim[1];

    outArray[7*2+0] = index+m_dim[0]*m_dim[1];
    outArray[7*2+1] = index+m_dim[0]+m_dim[0]*m_dim[1];

    outArray[8*2+0] = index;
    outArray[8*2+1] = index+m_dim[0]*m_dim[1];

    outArray[9*2+0] = index+1;
    outArray[9*2+1] = index+1+m_dim[0]*m_dim[1];

    outArray[10*2+0] = index+1+m_dim[0];
    outArray[10*2+1] = index+1+m_dim[0]+m_dim[0]*m_dim[1];

    outArray[11*2+0] = index+m_dim[0];
    outArray[11*2+1] = index+m_dim[0]+m_dim[0]*m_dim[1];

}

void Isosurface::GetEdgeVertices(double* verts, double* pVerts[12] ) {
    verts[0*3+0]=double(0.5); pVerts[0]=&verts[0*3+0];
    verts[0*3+1]=double(0);
    verts[0*3+2]=double(0);

    verts[1*3+0]=double(1);
    verts[1*3+1]=double(0);
    verts[1*3+2]=double(0.5); pVerts[1]=&verts[1*3+2];

    verts[2*3+0]=double(0.5); pVerts[2]=&verts[2*3+0];
    verts[2*3+1]=double(0);
    verts[2*3+2]=double(1); 

    verts[3*3+0]=double(0);
    verts[3*3+1]=double(0);
    verts[3*3+2]=double(0.5); pVerts[3]=&verts[3*3+2];

    verts[4*3+0]=double(0.5); pVerts[4]=&verts[4*3+0];
    verts[4*3+1]=double(1);
    verts[4*3+2]=double(0);

    verts[5*3+0]=double(1);
    verts[5*3+1]=double(1);
    verts[5*3+2]=double(0.5); pVerts[5]=&verts[5*3+2];

    verts[6*3+0]=double(0.5); pVerts[6]=&verts[6*3+0];
    verts[6*3+1]=double(1);
    verts[6*3+2]=double(1);

    verts[7*3+0]=double(0);
    verts[7*3+1]=double(1);
    verts[7*3+2]=double(0.5); pVerts[7]=&verts[7*3+2];

    verts[8*3+0]=double(0);
    verts[8*3+1]=double(0.5); pVerts[8]=&verts[8*3+1];
    verts[8*3+2]=double(0);

    verts[9*3+0]=double(1);
    verts[9*3+1]=double(0.5); pVerts[9]=&verts[9*3+1];
    verts[9*3+2]=double(0);

    verts[10*3+0]=double(1);
    verts[10*3+1]=double(0.5); pVerts[10]=&verts[10*3+1];
    verts[10*3+2]=double(1);

    verts[11*3+0]=double(0);
    verts[11*3+1]=double(0.5); pVerts[11]=&verts[11*3+1];
    verts[11*3+2]=double(1);
}
