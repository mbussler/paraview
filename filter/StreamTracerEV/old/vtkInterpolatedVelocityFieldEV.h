/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkInterpolatedVelocityFieldEV.h

=========================================================================*/
// .NAME vtkInterpolatedVelocityFieldEV - A concrete class for obtaining
//  the interpolated velocity values at a point in Eigenvector data.
//
//
// .SECTION Description
// The main functionality supported here is the point location inside vtkOverlappingAMR data set.


#ifndef __vtkInterpolatedVelocityFieldEV_h
#define __vtkInterpolatedVelocityFieldEV_h

#include <vtkAbstractInterpolatedVelocityField.h>
//BTX
#include <vector> // STL Header; Required for vector
//ETX


class vtkInterpolatedVelocityFieldEV : public vtkAbstractInterpolatedVelocityField
{
public:
  vtkTypeMacro( vtkInterpolatedVelocityFieldEV, vtkAbstractInterpolatedVelocityField );

  // Description:
  // Evaluate the velocity field f at point (x, y, z).
  virtual int FunctionValues( double * x, double * f )
  {
    vtkErrorMacro( <<"Not Implemented!" );
  };

  // Description:
  // Set the last cell id to -1 to incur a global cell search for the next point.
  void ClearLastCellId() { this->LastCellId = -1; }

  // Description:
  // Get the interpolation weights cached from last evaluation. Return 1 if the
  // cached cell is valid and 0 otherwise.
  int GetLastWeights( double * w );
  int GetLastLocalCoordinates( double pcoords[3] );

  void SetLastDirection( double dir[3]);
  
protected:
  vtkInterpolatedVelocityFieldEV();
  ~vtkInterpolatedVelocityFieldEV();

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
  
  double lastStepDirection[3];


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
  vtkInterpolatedVelocityFieldEV
    ( const vtkInterpolatedVelocityFieldEV & );  // Not implemented.
  void operator = ( const vtkInterpolatedVelocityFieldEV & );  // Not implemented.
};



#endif
