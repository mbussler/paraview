#ifndef _ISOSURFACE_H_
#define _ISOSURFACE_H_

//only dump components with size greater than
#define DISPLAY_MIN_SIZE 16
typedef unsigned int uint;

#include <vector>
/* #include <helper_math.h> */
#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <algorithm>
//#include <values.h>
#include <cstring>
#include "cuda_runtime.h"

#include "vtkSmartPointer.h" // compiler errors if this is forward declared
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkPolyDataAlgorithm.h" //superclass
#include "vtkImageData.h"
#include "vtkExecutive.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkTriangle.h"
#include "vtkUnsignedIntArray.h"
#include "linalg.h"

extern uint edgeTable[256];
extern uint triTable[256][16];
extern uint numVertsTable[256];

extern bool cudaRidgeData( uint dim[3], 
                           const double* data, 
                           const double* grad, 
                           const double* hesse, 
                           double** evals, 
                           uint* numVertices, 
                           double** vertices, 
                           double **pointData,
                           double featureThreshold, 
                           double origin[3], 
                           double spacing[3], 
                           unsigned 
                           int stencilRange, 
                           bool valley = false);

class Isosurface {
 public:
    Isosurface();
    ~Isosurface();
    uint m_dim[3];
    double *m_data;
    double *m_grad_data;
    double *m_hesse_data;
    void generateTriangles(int regThreshold, double angleThreshold, double minDataValue, double maxDataValue, double featureThreshold, double waveSpeed, unsigned int stencilRange);
    void generateGradients();
    vtkPolyData* GetPolyDataPointer() {return m_VTKPolyData;}
    //vtkPolyData* GetRigdeNodePointer() {return m_RidgeNodesPolyData;}
    vtkImageData* GetEigenValueGrid() {return m_EigenValues;}
    void CleanUp(void);
    void SetOrigin(double origin[3]) { m_origin[0] = origin[0]; m_origin[1] = origin[1]; m_origin[2] = origin[2]; }
    void GetOrigin(double origin[3]) { origin[0] = m_origin[0]; origin[1] = m_origin[1]; origin[2] = m_origin[2]; }
    void SetSpacing(double spacing[3]) { m_spacing[0] = spacing[0]; m_spacing[1] = spacing[1]; m_spacing[2] = spacing[2]; }
    void GetSpacing(double spacing[3]) { spacing[0] = m_spacing[0]; spacing[1] = m_spacing[1]; spacing[2] = m_spacing[2]; }
    void SetExtent(int extent[6]) { m_extent[0] = extent[0]; m_extent[1] = extent[1]; m_extent[2] = extent[2];
                                    m_extent[3] = extent[3]; m_extent[4] = extent[4]; m_extent[5] = extent[5];}
    void GetExtent(int extent[6]) { extent[0] = m_extent[0]; extent[1] = m_extent[1]; extent[2] = m_extent[2];
                                    extent[3] = m_extent[3]; extent[4] = m_extent[4]; extent[5] = m_extent[5];}
    void SetValley(bool valley){m_valley = valley;}
    bool GetValley(bool valley){ return valley;}

    void SetManifold(bool manifold){m_manifold = manifold;}
    bool GetManifold(bool manifold){ return manifold;}


 private:
    
     inline static double* normalize(double* v) 
     {
         double l=-sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); 
         v[0]=v[0]/l; v[1]=v[1]/l; v[2]=v[2]/l; 
         return v;
     }

    void mergeTriangles( float3* pos, float* data, 
                         int totalVerts, 
                         std::vector<unsigned int>& indices, 
                         std::vector<float3>& vertices, 
                         std::vector<float>& scalars );

    int labelComponents_rek( std::vector<unsigned int>* indices, 
                             bool manifold, int triangle, 
                             std::map<std::pair<int, int>, std::vector<int> > *edges, 
                             bool *triaVisited, 
                             int *flippedCnt,	
                             int compLabel,	
                             std::vector<int> *triangleComponents );  
    
    void labelComponents( std::vector<unsigned int>* indices, 
                          bool manifold, 
                          std::vector<int> *triangleLabels, 
                          std::vector<int> *objectSizes);

    bool m_hasGradients;
    vtkSmartPointer<vtkFloatArray> m_VTKdataGrads;
    vtkSmartPointer<vtkPoints> m_VTKdataPoints;
    vtkSmartPointer<vtkCellArray> m_VTKisoTriangles;
    vtkSmartPointer<vtkFloatArray> m_VTKisoGrads;
    vtkSmartPointer<vtkPolyData> m_VTKPolyData;
    vtkSmartPointer<vtkDoubleArray> m_VTKPointData;

    vtkSmartPointer<vtkImageData> m_EigenValues;
    bool generateRidgeSurface(uint* numVerts, double** verts, double **evals, double **pointData, double featureThreshold, unsigned int stencilRange);
    //void generateGradientsGPU();
    //void generateTrianglesCPU(double c, uint* numVerts, double** verts);
    //void generateTrianglesGPU(double c, uint* numVerts, double** verts, double** triGrads);
    double findValueOnEdge(double f1, double x1, double f2, double x2, double c);
    void GetEdgeIndices(int, int*);
    void GetEdgeVertices(double*, double*[12]);
    double m_origin[3];
    double m_spacing[3];
    int m_extent[6];
    bool m_valley;
    bool m_manifold;
};

class FKeys3 {
 public:
 FKeys3(float k1, float k2, float k3) : key1(k1), key2(k2), key3(k3) { }
    bool operator<(const FKeys3 &right) const {
	return (key1 < right.key1 || 
		(key1 == right.key1 && 
		 (key2 < right.key2 || 
		  (key2 == right.key2 && 
		   (key3 < right.key3))))
		);
    }
    float key1, key2, key3;
};

#endif // _ISOSURFACE_H_
