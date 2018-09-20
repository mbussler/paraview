/*=========================================================================

  Program:   Visualization Toolkit
  Module:    PeriMultiFit.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __PeriMultiFit_h
#define __PeriMultiFit_h

#include "vtkImageAlgorithm.h"
#include "vtkSmartPointer.h"

class vtkIdList;
class vtkPoints;
class vtkPolyData;
class vtkDoubleArray;

class PeriMultiFit : public vtkImageAlgorithm
{
  
  public:
   
    vtkTypeMacro(PeriMultiFit, vtkImageAlgorithm);
    void PrintSelf(ostream &os, vtkIndent indent);
    
    static PeriMultiFit *New();

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

	vtkSetMacro( FitConstantTerm	   , bool );
	vtkSetMacro( FitLinearTerms		   , bool );
	vtkSetMacro( FitHyperbolicTerms    , bool );
	vtkSetMacro( FitQuadraticTerms	   , bool );

	vtkGetMacro( FitConstantTerm	   , bool );
	vtkGetMacro( FitLinearTerms		   , bool );
	vtkGetMacro( FitHyperbolicTerms    , bool );
	vtkGetMacro( FitQuadraticTerms	   , bool );

   // Description:
  // These are basically a convenience method that calls SetInputArrayToProcess
  // to set the array used as the input scalars.  The fieldAssociation comes
  // from the vtkDataObject::FieldAssocations enum.  The fieldAttributeType
  // comes from the vtkDataSetAttributes::AttributeTypes enum.
  virtual void SetInputScalars(int fieldAssociation, const char *name);
  virtual void SetInputScalars(int fieldAssociation, int fieldAttributeType);

  protected:
    PeriMultiFit();
    ~PeriMultiFit();

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

  private:
    PeriMultiFit(const PeriMultiFit&);  // Not implemented.
    void operator=(const PeriMultiFit&);  // Not implemented.

    vtkSmartPointer<vtkIdList> GetConnectedVertices( 
        int id, vtkSmartPointer<vtkPolyData> mesh, 
        int fieldAssocation, vtkDataArray *array,
        vtkSmartPointer<vtkDoubleArray> values );

	void LeastSquaresFitting(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkDoubleArray> values, double coeff[]);
	double evaluateFunction(double pos[], double coeff[]);
	void evaluateGradient(double pos[], double coeff[], double grad[]);
};

#endif
;