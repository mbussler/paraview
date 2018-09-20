
#ifndef __vtkInterpolatedVelocityFieldPCA_h
#define __vtkInterpolatedVelocityFieldPCA_h


#include "vtkAbstractInterpolatedVelocityField.h"

class vtkInterpolatedVelocityFieldPCA : public vtkAbstractInterpolatedVelocityField
{
public:
  vtkTypeMacro( vtkInterpolatedVelocityFieldPCA, vtkFunctionSet );
  void PrintSelf( ostream & os, vtkIndent indent );

  // Description:
  // Set/Get the caching flag. If this flag is turned ON, there are two levels
  // of caching for derived concrete class vtkInterpolatedVelocityField and one
  // level of caching for derived concrete class vtkCellLocatorInterpolatedVelocityField.
  // Otherwise a global cell location is always invoked for evaluating the
  // function values at any point.
  vtkSetMacro( Caching, bool );
  vtkGetMacro( Caching, bool );

  // Description:
  // Get the caching statistics. CacheHit refers to the number of level #0 cache
  // hits while CacheMiss is the number of level #0 cache misses.
  vtkGetMacro( CacheHit, int );
  vtkGetMacro( CacheMiss, int );

  vtkGetObjectMacro( LastDataSet, vtkDataSet );

  // Description:
  // Get/Set the id of the cell cached from last evaluation.
  vtkGetMacro( LastCellId, vtkIdType );
  virtual void SetLastCellId( vtkIdType c ) { this->LastCellId = c; }

  // Description:
  // Set the id of the most recently visited cell of a dataset.
  virtual void SetLastCellId( vtkIdType c, int dataindex ) = 0;

  // Description:
  // Get/Set the name of a spcified vector array. By default it is NULL, with
  // the active vector array for use.
  vtkGetStringMacro( VectorsSelection );
  vtkGetMacro(VectorsType,int);

 // Description:
  // the association type (see vtkDataObject::FieldAssociations)
  // and the name of the velocity data field
  void SelectVectors(int fieldAssociation, const char * fieldName );


  // Description:
  // Set/Get the flag indicating vector post-normalization (following vector
  // interpolation). Vector post-normalization is required to avoid the
  // 'curve-overshooting' problem (caused by high velocity magnitude) that
  // occurs when Cell-Length is used as the step size unit (particularly the
  // Minimum step size unit). Furthermore, it is required by RK45 to achieve,
  // as expected, high numerical accuracy (or high smoothness of flow lines)
  // through adaptive step sizing. Note this operation is performed (when
  // NormalizeVector TRUE) right after vector interpolation such that the
  // differing amount of contribution of each node (of a cell) to the
  // resulting direction of the interpolated vector, due to the possibly
  // significantly-differing velocity magnitude values at the nodes (which is
  // the case with large cells), can be reflected as is. Also note that this
  // flag needs to be turned to FALSE after vtkInitialValueProblemSolver::
  // ComputeNextStep() as subsequent operations, e.g., vorticity computation,
  // may need non-normalized vectors.
  vtkSetMacro( NormalizeVector, bool );
  vtkGetMacro( NormalizeVector, bool );

  // Description:
  // Import parameters. Sub-classes can add more after chaining.
  virtual void CopyParameters( vtkInterpolatedVelocityFieldPCA * from )
    { this->Caching = from->Caching; }


  // Description:
  // Evaluate the velocity field f at point (x, y, z).
  virtual int FunctionValues( double * x, double * f ) = 0;

  // Description:
  // Set the last cell id to -1 to incur a global cell search for the next point.
  void ClearLastCellId() { this->LastCellId = -1; }

  // Description:
  // Get the interpolation weights cached from last evaluation. Return 1 if the
  // cached cell is valid and 0 otherwise.
  int GetLastWeights( double * w );
  int GetLastLocalCoordinates( double pcoords[3] );

protected:
  vtkInterpolatedVelocityFieldPCA();
  ~vtkInterpolatedVelocityFieldPCA();

  static const double TOLERANCE_SCALE;

  int       CacheHit;
  int       CacheMiss;
  int       WeightsSize;
  bool      Caching;
  bool      NormalizeVector;
  int       VectorsType;
  char *    VectorsSelection;
  double *  Weights;
  double    LastPCoords[3];
  vtkIdType LastCellId;
  vtkDataSet *     LastDataSet;
  vtkGenericCell * Cell;
  vtkGenericCell * GenCell; // the current cell


  // Description:
  // Set the name of a specific vector to be interpolated.
  vtkSetStringMacro( VectorsSelection );

  // Description:
  // Evaluate the velocity field f at point (x, y, z) in a specified dataset
  // by invoking vtkDataSet::FindCell() to locate the next cell if the given
  // point is outside the current cell. To address vtkPointSet, vtkPointLocator
  // is involved via vtkPointSet::FindCell() in vtkInterpolatedVelocityField
  // for cell location. In vtkCellLocatorInterpolatedVelocityField, this function
  // is invoked just to handle vtkImageData and vtkRectilinearGrid that are not
  // assigned with any vtkAbstractCellLocatot-type cell locator.
  virtual int FunctionValues( vtkDataSet * ds, double * x, double * f );

//BTX
  friend class vtkTemporalInterpolatedVelocityField;
  // Description:
  // If all weights have been computed (parametric coords etc all valid), a
  // scalar/vector can be quickly interpolated using the known weights and
  // the cached generic cell. This function is primarily reserved for use by
  // vtkTemporalInterpolatedVelocityField
  void FastCompute( vtkDataArray * vectors, double f[3] );
  bool InterpolatePoint( vtkPointData * outPD, vtkIdType outIndex );
  vtkGenericCell * GetLastCell()
    { return ( this->LastCellId != -1 ) ? this->GenCell : NULL; }
//ETX

private:
  vtkInterpolatedVelocityFieldPCA
    ( const vtkInterpolatedVelocityFieldPCA & );  // Not implemented.
  void operator = ( const vtkInterpolatedVelocityFieldPCA & );  // Not implemented.
};



#endif
