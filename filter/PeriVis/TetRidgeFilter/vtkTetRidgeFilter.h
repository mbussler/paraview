
#ifndef __vtkTetRidgeFilter_h
#define __vtkTetRidgeFilter_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkSmartPointer.h"
#include <map>
#include <utility>
#include <vector>

class vtkPolyData;
class vtkIdList;

typedef struct {
    double x,y,z;
} XYZ;

class vtkTetRidgeFilter : public vtkPolyDataAlgorithm
{
public:
    vtkTypeMacro(vtkTetRidgeFilter,vtkPolyDataAlgorithm);

    // Description:
    static vtkTetRidgeFilter *New();

    vtkSetMacro( MinDataValue, double);
    vtkGetMacro( MinDataValue, double);
    vtkSetMacro( MaxDataValue, double);
    vtkGetMacro( MaxDataValue, double);
    vtkSetMacro( Manifold, bool);
    vtkGetMacro( Manifold, bool);
    vtkSetMacro(RegionThreshold, int);
    vtkGetMacro(RegionThreshold, int);

protected:
    vtkTetRidgeFilter();
    ~vtkTetRidgeFilter();

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    void getCellNeighbors( vtkPolyData* mesh, vtkIdType cellId, std::vector<vtkIdType>& cellNeighbors);

private:

    vtkTetRidgeFilter(const vtkTetRidgeFilter&);  // Not implemented.
    void operator=(const vtkTetRidgeFilter&);  // Not implemented.

    double MaxDataValue;
    double MinDataValue;
    int RegionThreshold;
    bool Manifold;


};


#endif
