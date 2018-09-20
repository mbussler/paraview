#ifndef __vtkRidgeGrow_h
#define __vtkRidgeGrow_h

#include "vtkPolyDataAlgorithm.h" //superclass
#include "vtkSmartPointer.h"

#include <vector>
#include <set>

class vtkPolyData;
class vtkDoubleArray;

class vtkRidgeGrow : public vtkPolyDataAlgorithm
{
public:
    static vtkRidgeGrow *New();
    vtkTypeMacro(vtkRidgeGrow, vtkPolyDataAlgorithm)

    vtkSetMacro( MinDistance, double );
    vtkGetMacro( MinDistance, double );
    vtkSetMacro( ResetOnTimestep, int );
    vtkGetMacro( ResetOnTimestep, int );
    vtkSetMacro( ClearCellGrowth, bool );
    vtkGetMacro( ClearCellGrowth, bool );
    vtkSetMacro( TargetDistanceMethod, int );
    vtkGetMacro( TargetDistanceMethod, int );

protected:
    vtkRidgeGrow();
    ~vtkRidgeGrow();

    // Make sure the pipeline knows what type we expect as input
    int FillInputPortInformation( int port, vtkInformation* info );
    int FillOutputPortInformation( int port, vtkInformation* info );

    virtual int RequestData(vtkInformation *,
        vtkInformationVector **,
        vtkInformationVector *);

    void getCellNeighbors( vtkPolyData* mesh, vtkIdType cellId, std::vector<vtkIdType>& cellNeighbors);

    double MinDistance;
    int ResetOnTimestep;
    bool ClearCellGrowth;
    int TargetDistanceMethod; //!< point-to-point if 0, point-to-cell if 1

    std::vector<int> refPtIds;
    std::vector<bool> refCellIds;
    std::vector<bool> isolatedPt;
    vtkSmartPointer<vtkPolyData> growMesh;
};

#endif