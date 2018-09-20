#ifndef __vtkRidgeTime_h
#define __vtkRidgeTime_h

#include "vtkMultiTimeStepAlgorithm.h" //superclass

class vtkRidgeTime : public vtkMultiTimeStepAlgorithm
{
public:
    static vtkRidgeTime *New();
    vtkTypeMacro(vtkRidgeTime, vtkMultiTimeStepAlgorithm)

    vtkGetVectorMacro( MergeRange, int, 2 );
    void SetMergeRange( int in, int out);
    void SetMergeRange( int range[2]);

    vtkGetVectorMacro(Translation,double,3);
    void SetTranslation(double x, double y, double z);
    void SetTranslation(double xyz[3]);

protected:
    vtkRidgeTime();
    ~vtkRidgeTime();

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

    int MergeRange[2];
    double Translation[3];

};

#endif