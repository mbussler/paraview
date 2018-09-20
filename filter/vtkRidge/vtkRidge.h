#ifndef __vtkRidge_h
#define __vtkRidge_h


#include "vtkPolyDataAlgorithm.h" //superclass
#include "vtkSmartPointer.h" // compiler errors if this is forward declared

typedef unsigned int uint;

class vtkPolyData;
class vtkImageData;
class vtkTransform;
class vtkInformation;
class vtkInformationVector;
class vtkIterativeClosestPointTransform;

#define DUMPVECTOR(x) (cout << "(" << x[0] << "," << x[1] << "," << x[2] << ")" << endl)

class vtkRidge : public vtkPolyDataAlgorithm
{
 public:
  static vtkRidge *New();
  vtkTypeMacro(vtkRidge, vtkPolyDataAlgorithm)
  void PrintSelf(ostream &os, vtkIndent indent);
  //void SetIsoThreshold(double iT);

  vtkSetClampMacro(AngleThreshold, float, 0.0, 90.0)
  vtkGetMacro(AngleThreshold, float)

  vtkSetMacro(WaveSpeed, float)
  vtkGetMacro(WaveSpeed, float)

  vtkSetMacro(MinDataValue, float)
  vtkGetMacro(MinDataValue, float)

  vtkSetMacro(MaxDataValue, float)
  vtkGetMacro(MaxDataValue, float)

  vtkSetMacro(EigenValueThreshold, float)
  vtkGetMacro(EigenValueThreshold, float)

  vtkSetMacro(RegionThreshold, int)
  vtkGetMacro(RegionThreshold, int)

  vtkSetMacro(StencilRange, int)
  vtkGetMacro(StencilRange, int)

  vtkSetMacro(Valley, int)
  vtkGetMacro(Valley, int)

  vtkSetMacro(UseInput, bool);
  vtkGetMacro(UseInput, bool);

  vtkSetMacro(Manifold, bool);
  vtkGetMacro(Manifold, bool);

  //double isoThreshold;
 protected:
  vtkRidge();
  ~vtkRidge();
  uint dim[3];
  // Make sure the pipeline knows what type we expect as input
  int FillInputPortInformation( int port, vtkInformation* info );
  int FillOutputPortInformation( int port, vtkInformation* info );
  // Generate output
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *); //the function that makes this class work with the vtk pipeline

  void getData(vtkSmartPointer<vtkImageData> in_grid, double* &in_data, int numComponents, bool &allocated_in_data);

  float AngleThreshold;
  float WaveSpeed;
  float MaxDataValue;
  float MinDataValue;
  float EigenValueThreshold;
  int RegionThreshold;
  int StencilRange;
  int Valley;
  bool UseInput;
  bool Manifold;
  vtkSmartPointer<vtkImageData> m_evalImageData;
};
#endif
