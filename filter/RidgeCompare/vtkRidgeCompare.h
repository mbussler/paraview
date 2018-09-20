#ifndef __vtkRidgeCompare_h
#define __vtkRidgeCompare_h

#include "vtkMultiTimeStepAlgorithm.h" //superclass

class vtkRidgeCompare : public vtkMultiTimeStepAlgorithm
{
public:
    static vtkRidgeCompare *New();
    vtkTypeMacro(vtkRidgeCompare, vtkMultiTimeStepAlgorithm)

    vtkSetMacro( TargetDistanceMethod, int );
    vtkGetMacro( TargetDistanceMethod, int );

protected:
    vtkRidgeCompare();
    ~vtkRidgeCompare();

    // Make sure the pipeline knows what type we expect as input
    int FillInputPortInformation( int port, vtkInformation* info );
    int FillOutputPortInformation( int port, vtkInformation* info );

    // Generate output
    virtual int RequestDataObject(vtkInformation *,
        vtkInformationVector **,
        vtkInformationVector *);

    virtual int RequestUpdateExtent(vtkInformation *,
        vtkInformationVector **,
        vtkInformationVector *);
    virtual int RequestInformation(vtkInformation *,
                                    vtkInformationVector **,
                                    vtkInformationVector *);

    virtual int RequestData(vtkInformation *,
        vtkInformationVector **,
        vtkInformationVector *);

    int TargetDistanceMethod; //!< point-to-point if 0, point-to-cell if 1

};

#endif