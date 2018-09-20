#ifndef __vtkTimeGeometry_h
#define __vtkTimeGeometry_h

#include "vtkPolyDataAlgorithm.h" //superclass
#include <vtkSmartPointer.h>

class vtkPolyData;
class vtkDoubleArray;

class vtkTimeGeometry : public vtkPolyDataAlgorithm
{
public:
    static vtkTimeGeometry *New();
    vtkTypeMacro(vtkTimeGeometry, vtkPolyDataAlgorithm)

    vtkSetMacro( ResetOnTimestep, int );
    vtkGetMacro( ResetOnTimestep, int );

protected:
    vtkTimeGeometry();
    ~vtkTimeGeometry();

    virtual int RequestData(vtkInformation *,
        vtkInformationVector **,
        vtkInformationVector *);

    int ResetOnTimestep;
    double currentTime;
    vtkSmartPointer<vtkPolyData> growMesh;
};

#endif