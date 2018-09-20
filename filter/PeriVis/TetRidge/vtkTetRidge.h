
#ifndef __vtkTetRidge_h
#define __vtkTetRidge_h

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

class vtkTetRidge : public vtkPolyDataAlgorithm
{
public:
    vtkTypeMacro(vtkTetRidge,vtkPolyDataAlgorithm);

    // Description:
    static vtkTetRidge *New();

    vtkSetMacro( MergeTolerance, double);
    vtkGetMacro( MergeTolerance, double);
    vtkSetMacro( Valley, int);
    vtkGetMacro( Valley, int);
    vtkSetMacro( MinDataValue, double);
    vtkGetMacro( MinDataValue, double);
    vtkSetMacro( MaxDataValue, double);
    vtkGetMacro( MaxDataValue, double);
    vtkSetMacro( Manifold, bool);
    vtkGetMacro( Manifold, bool);
    vtkSetMacro( StrictFilter, bool);
    vtkGetMacro( StrictFilter, bool);
    vtkSetMacro(RegionThreshold, int);
    vtkGetMacro(RegionThreshold, int);
    vtkSetMacro(EigenValueThreshold, double);
    vtkGetMacro(EigenValueThreshold, double);

protected:
    vtkTetRidge();
    ~vtkTetRidge();

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    bool applyFilter(double val0, double val1, double val2, double ew0, double t, double ew1, double ew2);

    void getLargestEV(double *ev, double lev[3]);

private:

    XYZ VertexInterpolate( double isovalue, XYZ p1, XYZ p2, 
                           double valp1, double valp2, 
                           double data1, double data2, double& val, 
                           double ew1, double ew2, double& ew,
                           double d1[9], double d2[9], double d[9] );

    vtkTetRidge(const vtkTetRidge&);  // Not implemented.
    void operator=(const vtkTetRidge&);  // Not implemented.
    void InsertTriangle(XYZ np0, XYZ np1, XYZ np2, vtkPoints * newPoints, vtkCellArray * newTris);

    virtual int FillInputPortInformation(int port, vtkInformation* info);
    virtual int FillOutputPortInformation(int port, vtkInformation* info);
    void mergeTriangles( vtkPoints* points, 
        vtkDataArray* values, 
        vtkDataArray* data_in, 
        std::vector<unsigned int>& indices, 
        std::vector<XYZ>& vertices, 
        std::vector<double>& scalars,
        vtkDataArray* data_merged );

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
    double EigenValueThreshold;
    double MergeTolerance;
    int Valley;
    int RegionThreshold;
    bool Manifold;
    bool StrictFilter;


};

class FKeys3 {
public:
    FKeys3(float k1, float k2, float k3) 
        : key1(k1), key2(k2), key3(k3) { }
    bool operator<(const FKeys3 &right) const {
        return (key1 < right.key1 || 
            (key1 == right.key1 && 
            (key2 < right.key2 || 
            (key2 == right.key2 && 
            (key3 < right.key3))))
            );
    }
    //bool operator<(const FKeys3 &right) const {
    //    return (key1 < (right.key1+tolerance) || 
    //           (abs(key1-right.key1) <= tolerance && 
    //           (key2 < (right.key2+tolerance) || 
    //           (abs(key2-right.key2) <= tolerance && 
    //           (key3 < (right.key3+tolerance)))))
    //        );
    //}
    //bool operator<(const FKeys3 &right) const {
    //    return (key1 < (right.key1-tolerance) || 
    //        (abs(key1 - right.key1) < tolerance && 
    //        (key2 < (right.key2-tolerance) || 
    //        (abs(key2 - right.key2) < tolerance && 
    //        (key3 < right.key3-tolerance)))));
    //}
    //bool operator==(const FKeys3 &right) const {
    //    return (
    //        abs( key1 - right.key1 ) < tolerance &&
    //        abs( key2 - right.key2 ) < tolerance &&
    //        abs( key3 - right.key3 ) < tolerance  );
    //}
    float key1, key2, key3;

private:
    //double tolerance;
};

#endif
