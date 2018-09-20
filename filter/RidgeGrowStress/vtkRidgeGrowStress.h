#ifndef __vtkRidgeGrowStress_h
#define __vtkRidgeGrowStress_h

#include "vtkPolyDataAlgorithm.h" //superclass
#include "vtkSmartPointer.h"

#include <vector>
#include <set>

class vtkPolyData;
class vtkDoubleArray;

class vtkRidgeGrowStress : public vtkPolyDataAlgorithm
{
public:
    static vtkRidgeGrowStress *New();
    vtkTypeMacro(vtkRidgeGrowStress, vtkPolyDataAlgorithm)

    vtkSetMacro( MinDistance, double );
    vtkGetMacro( MinDistance, double );
    vtkSetMacro( MaxGrowthRate, double );
    vtkGetMacro( MaxGrowthRate, double );
    vtkSetMacro( ResetOnTimestep, int );
    vtkGetMacro( ResetOnTimestep, int );
    vtkSetMacro( ClearCellGrowth, bool );
    vtkGetMacro( ClearCellGrowth, bool );
    vtkSetMacro( InterpolateAdditionalData, bool );
    vtkGetMacro( InterpolateAdditionalData, bool );
    vtkSetStringMacro(ArrayName);
    vtkGetStringMacro(ArrayName);

protected:
    vtkRidgeGrowStress();
    ~vtkRidgeGrowStress();

    // Make sure the pipeline knows what type we expect as input
    int FillInputPortInformation( int port, vtkInformation* info );
    int FillOutputPortInformation( int port, vtkInformation* info );

    virtual int RequestData(vtkInformation *,
        vtkInformationVector **,
        vtkInformationVector *);

    void getCellNeighbors( vtkPolyData* mesh, vtkIdType cellId, std::vector<vtkIdType>& cellNeighbors);

    double MinDistance;
    double MaxGrowthRate;
    int ResetOnTimestep;
    bool ClearCellGrowth;
    bool InterpolateAdditionalData;
    char* ArrayName;

    std::vector<int> refPtIds;
    std::vector<bool> refCellIds;
    vtkSmartPointer<vtkPolyData> growMesh;
};

#endif