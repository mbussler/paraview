
#ifndef __vtkTetGen_h
#define __vtkTetGen_h

#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkTetGen : public vtkUnstructuredGridAlgorithm
{
public:
    vtkTypeMacro(vtkTetGen,vtkUnstructuredGridAlgorithm);

    // Description:
    static vtkTetGen *New();

protected:
    vtkTetGen();
    ~vtkTetGen();

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:

    vtkTetGen(const vtkTetGen&);  // Not implemented.
    void operator=(const vtkTetGen&);  // Not implemented.

    virtual int FillInputPortInformation(int port, vtkInformation* info);
    virtual int FillOutputPortInformation(int port, vtkInformation* info);
};

#endif
