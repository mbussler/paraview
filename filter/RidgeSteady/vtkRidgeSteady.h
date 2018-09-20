#ifndef __vtkRidgeSteady_h
#define __vtkRidgeSteady_h

#include "vtkPolyDataAlgorithm.h" //superclass
#include "vtkSmartPointer.h"

#include <vector>

class vtkPolyData;
class vtkDoubleArray;

class vtkRidgeSteady : public vtkPolyDataAlgorithm
{
public:
    static vtkRidgeSteady *New();
    vtkTypeMacro(vtkRidgeSteady, vtkPolyDataAlgorithm)

    vtkSetMacro( MinDistance, double );
    vtkGetMacro( MinDistance, double );
    vtkSetMacro( ResetOnTimestep, int );
    vtkGetMacro( ResetOnTimestep, int );
    vtkSetMacro( MergeRidges, bool );
    vtkGetMacro( MergeRidges, bool );

protected:
    vtkRidgeSteady();
    ~vtkRidgeSteady();

    // Make sure the pipeline knows what type we expect as input
    int FillInputPortInformation( int port, vtkInformation* info );
    int FillOutputPortInformation( int port, vtkInformation* info );

    virtual int RequestData(vtkInformation *,
        vtkInformationVector **,
        vtkInformationVector *);

    double MinDistance;
    int ResetOnTimestep;
    bool MergeRidges;

    vtkSmartPointer<vtkPolyData> oldMesh;
};

#endif