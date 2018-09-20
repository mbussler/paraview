
#ifndef __ImplicitTetClip_h
#define __ImplicitTetClip_h

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkUnstructuredGridAlgorithm.h"
#include <vtkSmartPointer.h>
#include <vtkPolygon.h>
#include <vector>

class vtkImplicitFunction;
class vtkIncrementalPointLocator;

class vtkIdList;

class ImplicitTetClip : public vtkUnstructuredGridAlgorithm
{
public:
  vtkTypeMacro(ImplicitTetClip,vtkUnstructuredGridAlgorithm);

  vtkGetMacro(InsideOut, bool);
  vtkSetMacro(InsideOut, bool);
  vtkGetMacro(ClipDistance, double);
  vtkSetMacro(ClipDistance, double);
  vtkGetMacro(Tolerance, double);
  vtkSetMacro(Tolerance, double);
  vtkGetMacro(ClipByCellEdges, bool);
  vtkSetMacro(ClipByCellEdges, bool);
  

  // Description:
  // Construct with user-specified implicit function; InsideOut turned off;
  // value set to 0.0; and generate clip scalars turned off.
  static ImplicitTetClip *New();

protected:
  ImplicitTetClip(vtkImplicitFunction *cf=NULL);
  ~ImplicitTetClip();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  bool IntersectWithPolys(vtkUnstructuredGrid * output, vtkIdList* cellPtIds, int ptId1, int ptId2, std::vector<vtkSmartPointer<vtkPolygon> > &polys);

private:
  ImplicitTetClip(const ImplicitTetClip&);  // Not implemented.
  void operator=(const ImplicitTetClip&);  // Not implemented.

  virtual int FillInputPortInformation(int port, vtkInformation* info);

  bool InsideOut;
  double ClipDistance;
  double Tolerance;
  bool ClipByCellEdges;

};

#endif
