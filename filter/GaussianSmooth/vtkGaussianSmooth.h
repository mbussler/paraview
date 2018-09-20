
#ifndef __vtkGaussianSmooth_h
#define __vtkGaussianSmooth_h

#include "vtkPointSetAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkGaussianSmooth : public vtkPointSetAlgorithm
{
public:
    vtkTypeMacro(vtkGaussianSmooth,vtkPointSetAlgorithm);

    // Description:
    static vtkGaussianSmooth *New();

    vtkSetMacro(SmoothRadius, double);
    vtkGetMacro(SmoothRadius, double);
    

protected:
    vtkGaussianSmooth();
    ~vtkGaussianSmooth();

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:

    vtkGaussianSmooth(const vtkGaussianSmooth&);  // Not implemented.
    void operator=(const vtkGaussianSmooth&);  // Not implemented.

    double SmoothRadius;
};

#endif
