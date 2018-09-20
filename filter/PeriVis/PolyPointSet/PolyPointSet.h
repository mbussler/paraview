
#ifndef __PolyPointSet_h
#define __PolyPointSet_h

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkUnstructuredGridAlgorithm.h"

class vtkImplicitFunction;
class vtkIncrementalPointLocator;

class vtkIdList;

class PolyPointSet : public vtkUnstructuredGridAlgorithm
{
public:
  vtkTypeMacro(PolyPointSet,vtkUnstructuredGridAlgorithm);

  vtkGetMacro(MaxClipDistance, double);
  vtkSetMacro(MaxClipDistance, double);
  vtkGetMacro(MinClipDistance, double);
  vtkSetMacro(MinClipDistance, double);
  // Description:
  // Construct with user-specified implicit function; InsideOut turned off;
  // value set to 0.0; and generate clip scalars turned off.
  static PolyPointSet *New();

protected:
  PolyPointSet(vtkImplicitFunction *cf=NULL);
  ~PolyPointSet();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  PolyPointSet(const PolyPointSet&);  // Not implemented.
  void operator=(const PolyPointSet&);  // Not implemented.

  virtual int FillInputPortInformation(int port, vtkInformation* info);

  double MaxClipDistance;
  double MinClipDistance;

};

#endif
