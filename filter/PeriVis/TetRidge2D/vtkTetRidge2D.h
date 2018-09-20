
#ifndef __vtkTetRidge2D_h
#define __vtkTetRidge2D_h

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

class vtkTetRidge2D : public vtkPolyDataAlgorithm
{
public:
    vtkTypeMacro(vtkTetRidge2D,vtkPolyDataAlgorithm);

    // Description:
    static vtkTetRidge2D *New();

    vtkSetMacro( IsoValue, double);
    vtkGetMacro( IsoValue, double);
    vtkSetMacro( Valley, int);
    vtkGetMacro( Valley, int);
    vtkSetMacro( MinDataValue, double);
    vtkGetMacro( MinDataValue, double);
    vtkSetMacro( MaxDataValue, double);
    vtkGetMacro( MaxDataValue, double);
    vtkSetMacro( Manifold, bool);
    vtkGetMacro( Manifold, bool);
    vtkSetMacro(RegionThreshold, int);
    vtkGetMacro(RegionThreshold, int);

protected:
    vtkTetRidge2D();
    ~vtkTetRidge2D();

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    void getLargestEV(double *ev, double lev[3]);

private:

    XYZ VertexInterpolate( double isovalue, XYZ p1, XYZ p2, double valp1, double valp2, double data1, double data2, double& val);

    vtkTetRidge2D(const vtkTetRidge2D&);  // Not implemented.
    void operator=(const vtkTetRidge2D&);  // Not implemented.
    void InsertLine(XYZ np0, XYZ np1, vtkPoints * newPoints, vtkCellArray * newLines);

    virtual int FillInputPortInformation(int port, vtkInformation* info);

    void mergeTriangles( vtkPoints* points, vtkDataArray* values, 
        std::vector<unsigned int>& indices, 
        std::vector<XYZ>& vertices, 
        std::vector<double>& scalars );

    void labelComponents( std::vector<unsigned int>* indices, 
                          bool manifold, 
                          std::vector<int> *triangleLabels, 
                          std::vector<int> *objectSizes);

    int labelComponents_rek( std::vector<unsigned int>* indices, 
                             bool manifold, 
                             int triangle, 
                             std::map<std::pair<int, int>, std::vector<int> > *edges, 
                             bool *triaVisited, 
                             int *flippedCnt, 
                             int compLabel, 
                             std::vector<int> *triangleComponents);

    double MaxDataValue;
    double MinDataValue;
    double IsoValue;
    int Valley;
    int RegionThreshold;
    bool Manifold;


};

class FKeys3 {
public:
    FKeys3(float k1, float k2, float k3) : key1(k1), key2(k2), key3(k3) { }
    bool operator<(const FKeys3 &right) const {
        return (key1 < right.key1 || 
            (key1 == right.key1 && 
            (key2 < right.key2 || 
            (key2 == right.key2 && 
            (key3 < right.key3))))
            );
    }
    float key1, key2, key3;
};

#endif
