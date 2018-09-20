
#ifndef __vtkLeastSquaresGradients_h
#define __vtkLeastSquaresGradients_h

#include "vtkSmartPointer.h"
#include "vtkPointSetAlgorithm.h"

class vtkLeastSquaresGradients : public vtkPointSetAlgorithm
{
public:
    vtkTypeMacro(vtkLeastSquaresGradients, vtkPointSetAlgorithm);

    // Description:
    static vtkLeastSquaresGradients *New();

    vtkSetMacro(Radius, double);
    vtkGetMacro(Radius, double);
    

protected:
    vtkLeastSquaresGradients();
    ~vtkLeastSquaresGradients();

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:

    vtkLeastSquaresGradients(const vtkLeastSquaresGradients&);  // Not implemented.
    void operator=(const vtkLeastSquaresGradients&);  // Not implemented.

    double Radius;
};

#endif
