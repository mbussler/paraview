#ifndef __vtkMeshComponents_h
#define __vtkMeshComponents_h

#include "vtkPolyDataAlgorithm.h" //superclass
#include "vtkSmartPointer.h"

#include <vector>
#include <set>

class vtkPolyData;
class vtkDoubleArray;

class vtkMeshComponents : public vtkPolyDataAlgorithm
{
public:
    static vtkMeshComponents *New();
    vtkTypeMacro(vtkMeshComponents, vtkPolyDataAlgorithm)

protected:
    vtkMeshComponents();
    ~vtkMeshComponents();

    virtual int RequestData(vtkInformation *,
        vtkInformationVector **,
        vtkInformationVector *);

    void getCellNeighbors( vtkPolyData* mesh, vtkIdType cellId, std::set<vtkIdType>& cellNeighbors);

};

#endif