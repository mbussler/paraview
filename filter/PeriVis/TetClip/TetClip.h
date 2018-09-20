
#ifndef __TetClip_h
#define __TetClip_h

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkImplicitFunction;
class vtkIncrementalPointLocator;
class vtkIdList;

class TetClip : public vtkUnstructuredGridAlgorithm
{
public:
  vtkTypeMacro(TetClip,vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct with user-specified implicit function; InsideOut turned off;
  // value set to 0.0; and generate clip scalars turned off.
  static TetClip *New();

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
  TetClip(vtkImplicitFunction *cf=NULL);
  ~TetClip();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  vtkImplicitFunction *ClipFunction;

  int InsideOut;
  int CrincleClip;

private:

  TetClip(const TetClip&);  // Not implemented.
  void operator=(const TetClip&);  // Not implemented.
};

#endif
