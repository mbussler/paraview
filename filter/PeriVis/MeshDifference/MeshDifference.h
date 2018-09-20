
#ifndef __MeshDifference_h
#define __MeshDifference_h

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkImplicitFunction;
class vtkIncrementalPointLocator;
class vtkPolyData;
class vtkIdList;

class MeshDifference : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(MeshDifference,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // 
  static MeshDifference *New();


protected:
  MeshDifference();
  ~MeshDifference();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:

  MeshDifference(const MeshDifference&);  // Not implemented.
  void operator=(const MeshDifference&);  // Not implemented.
};

#endif
