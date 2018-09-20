#ifndef _ISOSURFACE_H_
#define _ISOSURFACE_H_

#include "vtkSmartPointer.h" // compiler errors if this is forward declared

class vtkFloatArray;
class vtkPoints;
class vtkCellArray;
class vtkPolyData;

typedef unsigned int uint;

///////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Creates isosurface                                                                       ///
/// Creates an isosurface with a threshold passed in ParaView 'Create isosurface at' variable and   ///
/// scalar data on GPU using Marching Cubes algorithm.                                              ///
/// @param[in] c Isosurface threshold value                                                         ///
/// @param[in] data Scalar data set                                                                 ///
/// @param[in] dim Dimensions of data set in x,y,z                                                  ///
/// @param[out] numVerts Address of uint to store number of generated vertices                       ///
/// @param[out] vertices Address of pointer to receive array of doubles containing vertices x,y,z    ///
/// @param[out] gradients Address of pointer to receive array of doubles containing vertices x,y,z    ///
///////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" bool cudaIsosurface(const double c, const double* data, double* grads, uint dim[3], uint* numVerts, double** vertices, double** gradients);

///////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Creates isosurface                                                                       ///
/// Creates an array of gradient vectors based on the scalar field using GPGPU                      ///
/// @param[in] data Scalar data set                                                                 ///
/// @param[in] dim Dimensions of data set in x,y,z                                                  ///
/// @param[out] gradients Adress of pointer to receive array of doubles containing gradients x,y,z  ///
///////////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" bool cudaGradient(const double* data, uint dim[3], double** gradients);

class Isosurface {
 public:
    Isosurface();
    ~Isosurface();
    uint m_dim[3];
    double *m_data;
    double *m_grad_data;
    void generateTriangles(double c);
    void generateGradients();
    vtkPolyData* GetPolyDataPointer() {return m_polydata;}
 private:
    bool m_hasGradients;
    vtkSmartPointer<vtkFloatArray> m_grads;
    vtkSmartPointer<vtkPoints> m_points;
    vtkSmartPointer<vtkCellArray> m_tris;
    vtkSmartPointer<vtkPolyData> m_polydata;
    void generateGradientsGPU();
    void generateTrianglesCPU(double c, uint* numVerts, double** verts);
    void generateTrianglesGPU(double c, uint* numVerts, double** verts, double** triGrads);
    double findValueOnEdge(double f1, double x1, double f2, double x2, double c);
    void GetEdgeIndices(int, int*);
    void GetEdgeVertices(double*, double*[12]);
};

#endif // _ISOSURFACE_H_
