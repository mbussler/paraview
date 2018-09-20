#ifndef __vtkRidgeMerge_h
#define __vtkRidgeMerge_h

#include "vtkMultiTimeStepAlgorithm.h" //superclass

class vtkRidgeMerge : public vtkMultiTimeStepAlgorithm
{
public:
    static vtkRidgeMerge *New();
    vtkTypeMacro(vtkRidgeMerge, vtkMultiTimeStepAlgorithm)

    vtkSetMacro( MergeRange, int );
    vtkGetMacro( MergeRange, int );
    vtkSetMacro( StartTimestep, int);
    vtkGetMacro( StartTimestep, int);
    vtkSetMacro( MinDistance, double );
    vtkGetMacro( MinDistance, double );

protected:
    vtkRidgeMerge();
    ~vtkRidgeMerge();

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

    int MergeRange;
    int StartTimestep;
    double MinDistance;

};

#endif