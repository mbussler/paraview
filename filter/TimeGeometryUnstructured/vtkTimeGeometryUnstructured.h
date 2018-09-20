#ifndef __vtkTimeGeometryUnstructured_h
#define __vtkTimeGeometryUnstructured_h

#include "vtkUnstructuredGridAlgorithm.h" //superclass
#include <vtkSmartPointer.h>

class vtkPolyData;
class vtkDoubleArray;

class vtkTimeGeometryUnstructured : public vtkUnstructuredGridAlgorithm
{
public:
    static vtkTimeGeometryUnstructured *New();
    vtkTypeMacro(vtkTimeGeometryUnstructured, vtkUnstructuredGridAlgorithm)

    vtkSetMacro( ResetOnTimestep, int );
    vtkGetMacro( ResetOnTimestep, int );

protected:
    vtkTimeGeometryUnstructured();
    ~vtkTimeGeometryUnstructured();

    virtual int RequestData(vtkInformation *,
        vtkInformationVector **,
        vtkInformationVector *);

    int ResetOnTimestep;
    double currentTime;
    vtkSmartPointer<vtkUnstructuredGrid> growMesh;
};

#endif