
#ifndef __ImplicitPolydataClip_h
#define __ImplicitPolydataClip_h

#include "vtkFiltersCoreModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"

class vtkImplicitFunction;
class vtkIncrementalPointLocator;

class ImplicitPolydataClip : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(ImplicitPolydataClip,vtkPolyDataAlgorithm);

  vtkGetMacro(InsideOut, bool);
  vtkSetMacro(InsideOut, bool);

  // Description:
  // Construct with user-specified implicit function; InsideOut turned off;
  // value set to 0.0; and generate clip scalars turned off.
  static ImplicitPolydataClip *New();

protected:
  ImplicitPolydataClip(vtkImplicitFunction *cf=NULL);
  ~ImplicitPolydataClip();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  ImplicitPolydataClip(const ImplicitPolydataClip&);  // Not implemented.
  void operator=(const ImplicitPolydataClip&);  // Not implemented.

  bool InsideOut;

};

#endif
