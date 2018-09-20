/*=========================================================================

  Program:   Visualization Toolkit
  Module:    PeriMultiFitTaylor.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __PeriMultiFitTaylor_h
#define __PeriMultiFitTaylor_h

#include "vtkImageAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkIdList;
class vtkPoints;
class vtkPolyData;
class vtkDoubleArray;

#ifdef _WIN32
    #ifdef PeriMultiFitTaylor_EXPORTS
        #define PeriMultiFitTaylor_API __declspec(dllexport)
    #else
        #define PeriMultiFitTaylor_API __declspec(dllimport)
    #endif
#else
    #define PeriMultiFitTaylor_API 
#endif

class PeriMultiFitTaylor_API PeriMultiFitTaylor : public vtkImageAlgorithm
{
  
  public:
   
    vtkTypeMacro(PeriMultiFitTaylor, vtkImageAlgorithm);
    void PrintSelf(ostream &os, vtkIndent indent);
    
    static PeriMultiFitTaylor *New();

    // Description:
    // Specify i-j-k dimensions on which to sample input points.
    vtkGetVectorMacro(SampleDimensions,int,3);

    // Description:
    // Set the i-j-k dimensions on which to sample the distance function.
    void SetSampleDimensions(int i, int j, int k);

    // Description:
    // Set the i-j-k dimensions on which to sample the distance function.
    void SetSampleDimensions(int dim[3]);

    //vtkSetMacro(FittingFunction, int);
    //vtkGetMacro(FittingFunction, int);

    vtkSetMacro( FitConstantTerm       , bool );
    vtkSetMacro( FitLinearTerms        , bool );
    vtkSetMacro( FitHyperbolicTerms    , bool );
    vtkSetMacro( FitQuadraticTerms     , bool );

    vtkGetMacro( FitConstantTerm       , bool );
    vtkGetMacro( FitLinearTerms        , bool );
    vtkGetMacro( FitHyperbolicTerms    , bool );
    vtkGetMacro( FitQuadraticTerms     , bool );

    vtkSetMacro( Interpolate           , bool );
    vtkGetMacro( Interpolate           , bool );

    vtkSetMacro( OutputGradient        , bool );
    vtkGetMacro( OutputGradient        , bool );

    vtkSetMacro( OutputHessian         , bool );
    vtkGetMacro( OutputHessian         , bool );

    vtkSetMacro( FitAlternativeFunction, bool );
    vtkGetMacro( FitAlternativeFunction, bool );

    

   // Description:
  // These are basically a convenience method that calls SetInputArrayToProcess
  // to set the array used as the input scalars.  The fieldAssociation comes
  // from the vtkDataObject::FieldAssocations enum.  The fieldAttributeType
  // comes from the vtkDataSetAttributes::AttributeTypes enum.
  virtual void SetInputScalars(int fieldAssociation, const char *name);
  virtual void SetInputScalars(int fieldAssociation, int fieldAttributeType);

  protected:
    PeriMultiFitTaylor();
    ~PeriMultiFitTaylor();

    virtual int RequestInformation ( vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

    virtual int ProcessRequest(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector);

    // see algorithm for more info
    virtual int FillInputPortInformation(int port, vtkInformation* info);

    double ModelBounds[6];
    int SampleDimensions[3];
    //int FittingFunction;

    bool FitConstantTerm;
    bool FitLinearTerms;
    bool FitHyperbolicTerms;
    bool FitQuadraticTerms;
    bool Interpolate;
    bool OutputHessian;
    bool OutputGradient;

    bool FitAlternativeFunction;

  private:
    PeriMultiFitTaylor(const PeriMultiFitTaylor&);  // Not implemented.
    void operator=(const PeriMultiFitTaylor&);  // Not implemented.

    vtkSmartPointer<vtkIdList> GetConnectedVertices( 
        int id, vtkSmartPointer<vtkPolyData> mesh, 
        int fieldAssocation, vtkDataArray *array,
        vtkSmartPointer<vtkDoubleArray> values );

    void LeastSquaresFitting(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkDoubleArray> values, double coeff[]);
    double evaluateFunction(double pos[], double coeff[]);
    void evaluateGradient(double pos[], double coeff[], double grad[]);
    void evaluateHessian(double a[], double coeff[], double H_f_a[]);

    void ComputeModelBounds(double origin[3], double spacing[3]);

};

#endif
;
