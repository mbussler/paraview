/*=========================================================================

  Program:   Visualization Toolkit
  Module:    PeriMultiFitHesse.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __PeriMultiFitHesse_h
#define __PeriMultiFitHesse_h

#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkIdList;
class vtkPoints;
class vtkPolyData;
class vtkDoubleArray;

class PeriMultiFitHesse : public vtkUnstructuredGridAlgorithm
{
  
  public:
   
    vtkTypeMacro(PeriMultiFitHesse, vtkUnstructuredGridAlgorithm);
    void PrintSelf(ostream &os, vtkIndent indent);
    
    static PeriMultiFitHesse *New();

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

   // Description:
  // These are basically a convenience method that calls SetInputArrayToProcess
  // to set the array used as the input scalars.  The fieldAssociation comes
  // from the vtkDataObject::FieldAssocations enum.  The fieldAttributeType
  // comes from the vtkDataSetAttributes::AttributeTypes enum.
  virtual void SetInputScalars(int fieldAssociation, const char *name);
  virtual void SetInputScalars(int fieldAssociation, int fieldAttributeType);

  protected:
    PeriMultiFitHesse();
    ~PeriMultiFitHesse();

    //virtual int RequestInformation ( vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
//     int FillInputPortInformation(int port, vtkInformation *info);

    virtual int ProcessRequest(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector);

    // see algorithm for more info
    virtual int FillInputPortInformation(int port, vtkInformation* info);

    bool FitConstantTerm;
    bool FitLinearTerms;
    bool FitHyperbolicTerms;
    bool FitQuadraticTerms;
    bool Interpolate;

  private:
    PeriMultiFitHesse(const PeriMultiFitHesse&);  // Not implemented.
    void operator=(const PeriMultiFitHesse&);  // Not implemented.

    vtkSmartPointer<vtkIdList> GetConnectedVertices( 
        int id, vtkSmartPointer<vtkPolyData> mesh, 
        int fieldAssocation, vtkDataArray *array,
        vtkSmartPointer<vtkDoubleArray> values );

    void LeastSquaresFitting(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkDoubleArray> values, double coeff[]);
    double evaluateFunction(double pos[], double coeff[]);
    void evaluateGradient(double pos[], double coeff[], double grad[]);
    void evaluateHessian(double a[], double coeff[], double H_f_a[]);
};

#endif
;