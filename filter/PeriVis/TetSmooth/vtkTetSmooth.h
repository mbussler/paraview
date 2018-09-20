
#ifndef __vtkTetSmooth_h
#define __vtkTetSmooth_h

#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkTetSmooth : public vtkUnstructuredGridAlgorithm
{
public:
    vtkTypeMacro(vtkTetSmooth,vtkUnstructuredGridAlgorithm);

    // Description:
    static vtkTetSmooth *New();

    vtkSetMacro(CenterWeight, int);
    vtkGetMacro(CenterWeight, int);
    vtkSetMacro(NumberOfIterations, int);
    vtkGetMacro(NumberOfIterations, int);
    vtkSetMacro(SmoothingMethod, int);
    vtkGetMacro(SmoothingMethod, int);
    vtkSetMacro(SmoothRadius, double);
    vtkGetMacro(SmoothRadius, double);
    

protected:
    vtkTetSmooth();
    ~vtkTetSmooth();

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:

    vtkTetSmooth(const vtkTetSmooth&);  // Not implemented.
    void operator=(const vtkTetSmooth&);  // Not implemented.

    virtual int FillInputPortInformation(int port, vtkInformation* info);

    virtual int FillOutputPortInformation(int port, vtkInformation* info);

    int SmoothingMethod;
    int CenterWeight;
    int NumberOfIterations;
    double SmoothRadius;
};

#endif
