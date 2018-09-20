
#ifndef __vtkTetTrimSurface_h
#define __vtkTetTrimSurface_h

#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkTetTrimSurface : public vtkUnstructuredGridAlgorithm
{
public:
    vtkTypeMacro(vtkTetTrimSurface,vtkUnstructuredGridAlgorithm);

    // Description:
    static vtkTetTrimSurface *New();

    vtkSetMacro(TrimDepth, int);
    vtkGetMacro(TrimDepth, int);

protected:
    vtkTetTrimSurface();
    ~vtkTetTrimSurface();

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:

    vtkTetTrimSurface(const vtkTetTrimSurface&);  // Not implemented.
    void operator=(const vtkTetTrimSurface&);  // Not implemented.

    int TrimDepth;
};

#endif
