
#ifndef __ClipLines_h
#define __ClipLines_h

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkImplicitFunction;
class vtkIncrementalPointLocator;
class vtkPolyData;
class vtkIdList;

class ClipLines : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(ClipLines,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct with user-specified implicit function; InsideOut turned off;
  // value set to 0.0; and generate clip scalars turned off.
  static ClipLines *New();

  // Description:
  // Set/Get the InsideOut flag. When off, a vertex is considered
  // inside the implicit function if its value is greater than the
  // Value ivar. When InsideOutside is turned on, a vertex is
  // considered inside the implicit function if its implicit function
  // value is less than or equal to the Value ivar.  InsideOut is off
  // by default.
  vtkSetMacro(InsideOut,int);
  vtkGetMacro(InsideOut,int);
  vtkBooleanMacro(InsideOut,int);

  vtkSetMacro(CrincleClip,int);
  vtkGetMacro(CrincleClip,int);
  vtkBooleanMacro(CrincleClip,int);

  // Description
  // Specify the implicit function with which to perform the
  // clipping. If you do not define an implicit function, then the input
  // scalar data will be used for clipping.
  virtual void SetClipFunction(vtkImplicitFunction*);
  vtkGetObjectMacro(ClipFunction,vtkImplicitFunction);

  // Description:
  // Return the mtime also considering the locator and clip function.
  unsigned long GetMTime();

protected:
  ClipLines(vtkImplicitFunction *cf=NULL);
  ~ClipLines();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  vtkImplicitFunction *ClipFunction;

  int InsideOut;
  int CrincleClip;

private:

  vtkSmartPointer<vtkIdList> getConnectedPoints( vtkPolyData * pd, vtkIdType ptId, vtkSmartPointer<vtkIdList> bondIds = NULL);

  ClipLines(const ClipLines&);  // Not implemented.
  void operator=(const ClipLines&);  // Not implemented.
};

#endif
