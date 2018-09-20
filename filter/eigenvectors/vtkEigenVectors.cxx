/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkEigenVectors.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkEigenVectors.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkDataSet.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkTransform.h"

#include "linalg.h"

//#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>
//#include <boost/concept_check.hpp>

/* performance measure */
#include "timer.h"
#include <QElapsedTimer>

//using Eigen::MatrixXf;
//using Eigen::VectorXf;

vtkStandardNewMacro(vtkEigenVectors);

inline bool SameSign(float a, float b) {
    return a*b >= 0.0f;
}

inline void swap( int& a, int& b) {
    int swp=a; a=b; b=swp;
};

template <class T>
inline void swap( T& a, T& b) {
    T swp=a; a=b; b=swp;
};

//==============================================================================
vtkEigenVectors::vtkEigenVectors()
{
  this->SetNumberOfInputPorts(1);
  
  // by default, process active point tensors
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::TENSORS);

  // by default, process active point scalars
  this->SetInputArrayToProcess(1, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::SCALARS);

  this->Output3EV = true;
  
  this->EigenvalueMethod = 0;
  
  this->FixEVs = false;
  
  this->AnisotrophyMeasure = false;
}
//==============================================================================
vtkEigenVectors::~vtkEigenVectors()
{
  // clean up
}

//==============================================================================
int vtkEigenVectors::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkDataSet *input = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataSet *output = vtkDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDataArray *inTensors;
  double tensor[9];
  vtkIdType numPts, inPtId, ptIncr;
  
  double *m[3], w[3], *v[3];
  double m0[3], m1[3], m2[3];
  double v0[3], v1[3], v2[3];
  double xv[3], yv[3], zv[3];

  // set up working matrices
  m[0] = m0; m[1] = m1; m[2] = m2;
  v[0] = v0; v[1] = v1; v[2] = v2;

  float m11,m12,m21,m22;
  
  vtkDebugMacro(<<"Calculating Eigen vectors");

  vtkPointData *outPD = output->GetPointData();
  inTensors = this->GetInputArrayToProcess(0, inputVector);
  numPts = input->GetNumberOfPoints();

  if ( !inTensors || numPts < 1 )
    {
    vtkErrorMacro(<<"No data!");
    return 1;
    }

  // check tensor dimensions
  bool b2DTensor = inTensors->GetNumberOfComponents() == 4;
  bool b3DTensor = inTensors->GetNumberOfComponents() == 9;
  
    
  if ( !b2DTensor && !b3DTensor  )
    {
    vtkErrorMacro(<<"Input data array has "<<inTensors->GetNumberOfComponents()<<" components. Only 2D (4 components) and 3D (9 components) tensors are supported!");
    return 1;
    }

  int n = b3DTensor || this->Output3EV ? 3 : 2;

  // only copy scalar data through
  vtkDoubleArray* ev0 = vtkDoubleArray::New();
  vtkDoubleArray* ev1 = vtkDoubleArray::New();
  vtkDoubleArray* ev2 = vtkDoubleArray::New();
  
  ev0->SetName("ev0");
  ev0->SetNumberOfComponents(n);
  ev0->SetNumberOfTuples(numPts);

  ev1->SetName("ev1");
  ev1->SetNumberOfComponents(n);
  ev1->SetNumberOfTuples(numPts);

  if( b3DTensor || this->Output3EV ) 
  {
    ev2->SetName("ev2");
    ev2->SetNumberOfComponents(n);
    ev2->SetNumberOfTuples(numPts);

  }
  
  vtkDoubleArray* ew0 = vtkDoubleArray::New();
  vtkDoubleArray* ew1 = vtkDoubleArray::New();
  vtkDoubleArray* ew2 = vtkDoubleArray::New();

  ew0->SetName("ew0");
  ew0->SetNumberOfComponents(1);
  ew0->SetNumberOfTuples(numPts);
  ew1->SetName("ew1");
  ew1->SetNumberOfComponents(1);
  ew1->SetNumberOfTuples(numPts);

  vtkFloatArray* ani_linear;
  vtkFloatArray* ani_planar;
  vtkFloatArray* isotrophy; 

  if( b3DTensor || this->Output3EV ) 
  {
    ew2->SetName("ew2");
    ew2->SetNumberOfComponents(1);
    ew2->SetNumberOfTuples(numPts);

    if( this->AnisotrophyMeasure){
      ani_linear = vtkFloatArray::New();
      ani_linear->SetName("linear anisotrophy measure");
      ani_linear->SetNumberOfComponents(1);
      ani_linear->SetNumberOfTuples(numPts);

      ani_planar = vtkFloatArray::New();
      ani_planar->SetName("planar anisotrophy measure");
      ani_planar->SetNumberOfComponents(1);
      ani_planar->SetNumberOfTuples(numPts);

      isotrophy = vtkFloatArray::New();
      isotrophy->SetName("isotrophy measure");
      isotrophy->SetNumberOfComponents(1);
      isotrophy->SetNumberOfTuples(numPts);
    }
  }
 
  QElapsedTimer timer;
  timer.start();
 
  //
  // Traverse all Input points, calculate and store eigen vectors
  //

  //Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> es;
  
  n = b3DTensor ? 3 : 2; // evaluate according to input data
  for (inPtId=0; inPtId < numPts; inPtId++)
    {
    // Translation is postponed

    inTensors->GetTuple(inPtId, tensor);

    if( this->EigenvalueMethod == 0) // Use linalg
    {
      if( n == 2)
      {
        mat2 m;

        for( int i=0; i<n; i++) { // row index 
          for( int j=0; j<n; j++) { // col index
            m[i][j] = tensor[ n*i+j];
          }
        }
        vec2 lambda;
        mat2eigenvalues( m, lambda);
        
        mat2 ev;
        mat2realEigenvector(m, lambda[0], ev[0]);
        mat2realEigenvector(m, lambda[1], ev[1]);

        // find largest EW
        int sorted[2];
        sorted[0]=0; sorted[1]=1;
        
        if( lambda[sorted[0]] < lambda[sorted[1]]) {
          swap(sorted[0], sorted[1]);
        }

        // copy result
        for( int i=0; i<n; i++) {
          w[i] = lambda[sorted[i]];
          for( int j=0; j<n; j++) {
            v[i][j] = ev[sorted[i]][j];        
          }
        }
      }

      else if( n == 3)
      {
        mat3 m;

        for( int i=0; i<n; i++) { // row index 
          for( int j=0; j<n; j++) { // col index
            m[i][j] = tensor[ n*i+j];
          }
        }
        vec3 lambda;
        mat3eigenvalues( m, lambda);
        
        mat3 ev;
        mat3realEigenvector(m, lambda[0], ev[0]);
        mat3realEigenvector(m, lambda[1], ev[1]);
        mat3realEigenvector(m, lambda[2], ev[2]);
      
        // find largest EW
        int sorted[3];
        sorted[0]=0; sorted[1]=1; sorted[2]=2;

        if( lambda[sorted[0]] < lambda[sorted[1]])
          swap(sorted[0],sorted[1]);
        if( lambda[sorted[1]] < lambda[sorted[2]])
          swap(sorted[1], sorted[2]);
        if( lambda[sorted[0]] < lambda[sorted[1]])
          swap(sorted[0],sorted[1]);

        // copy result
        for( int i=0; i<n; i++) {
          w[i] = lambda[sorted[i]];
          for( int j=0; j<n; j++) {
            v[i][j] = ev[sorted[i]][j];
          }
        }
      }      
    }
  
    else if( this->EigenvalueMethod == 1) // Use Eigen
    {
    //  
    //  MatrixXf m(n,n);
    //  for( int i=0; i<n; i++) // row index
    //  {
    //    for( int j=0; j<n; j++) // col index
    //    {
    //      m(i,j) = tensor[ n*i+j];
    //    }
    //  }

    //  es.compute(m,Eigen::ComputeEigenvectors);
    //  
    //  VectorXf ew = es.eigenvalues();
    //  MatrixXf ev = es.eigenvectors();

    //  if( this->FixEVs ) 
    //  {
    //    for( int i=0; i<2; i++) 
    //    {
    //      float c11 = m(0,0)-ew(i);
    //      float c12 = m(0,1);
    //      float c21 = m(0,1);
    //      float c22 = m(1,1)-ew(i);
    //      
    //      // just to be sure. 
    //      assert( SameSign(c11, c12) ==  SameSign(c21, c22));

    //      // now check EVs and flip component if necessary
    //      float ev1 = ev(i,0);
    //      float ev2 = ev(i,1);
    //      
    //      if(  (  SameSign(c11, c12) &&  SameSign(ev1, ev2)) ||
    //            ( !SameSign(c11, c12) && !SameSign(ev1, ev2))  )
    //      {
    //          ev(i,0) *=-1;
    //      }
    //    }
    //  }

    //  // find largest EW
    //  int sorted[3];
    //  for( int i=0; i<n; i++) {
    //    sorted[i] = i;
    //  }
    //  
    //  // sort by descending order
    //  bool swapped=true;
    //  while(swapped) 
    //  {
    //    swapped=false;
    //    for( int i=0; i<(n-1); i++) {
    //      if( ew(sorted[i]) < ew(sorted[i+1]) ) {
    //        int swp = sorted[i];
    //        sorted[i] = sorted[i+1];
    //        sorted[i+1] = swp;
    //        swapped = true;
    //      }
    //    }
    //  }

    //  // copy result
    //  for( int i=0; i<n; i++)
    //  {
    //    for( int j=0; j<n; j++)
    //    {
    //      v[i][j] = ev(sorted[i],j);        
    //    }
    //    w[i] = ew(sorted[i]);
    //  }
    } 
    else // use vtkMath::JacobiN
    {
      for (int i=0; i<n; i++) // row index
      {
        for (int j=0; j<n; j++) // col index
        {
          m[i][j] = tensor[n*i+j];
        }
      }
      
      m11 = m[0][0];
      m12 = m[0][1];
      m21 = m[1][0];
      m22 = m[1][1];
      
      vtkMath::JacobiN(m, n, w, v);

      if( this->FixEVs ) 
      {
        for( int n=0; n<2; n++) 
        {
          float c11 = m11-w[n];
          float c12 = m12;
          float c21 = m21;
          float c22 = m22-w[n];
          
          // just to be sure. 
          assert( SameSign(c11, c12) ==  SameSign(c21, c22));

          // now check EVs and flip component if necessary
          float ev1 = v[n][0];
          float ev2 = v[n][1];
          
          if(  (  SameSign(c11, c12) &&  SameSign(ev1, ev2)) ||
                ( !SameSign(c11, c12) && !SameSign(ev1, ev2))  )
          {
              v[n][0] *=-1;
          }
        }
      }
    }
    
    if( b2DTensor && this->Output3EV )
    {
      // add third component
      v0[2]=v1[2]=0;                     // set last component of EVs to 0
      v2[0]=0; v2[1]=0; v2[2]=1; w[2]=1; // set third EV to (0,0,1), EW to 1
    }

    //copy eigenvectors
    ev0->SetTupleValue(inPtId, v0);
    ew0->SetTupleValue(inPtId, &w[0]);
    ev1->SetTupleValue(inPtId, v1);
    ew1->SetTupleValue(inPtId, &w[1]);
  
    if( b3DTensor || this->Output3EV )
      {
      ev2->SetTupleValue(inPtId, v2);
      ew2->SetTupleValue(inPtId, &w[2]);

        if( this->AnisotrophyMeasure ) 
        {
          //sort by absolute values in descending order
          //float w_0 = std::abs(w[0]);
          //float w_1 = std::abs(w[1]);
          //float w_2 = std::abs(w[2]);
          //if( w0 < w1)  swap<float>(w0, w1);
          //if( w1 < w2)  swap<float>(w1, w2);
          //if( w0 < w1)  swap<float>(w0, w1);

          //float sum = w_0+w_1+w_2;
          //float linear =      (w_0-w_1) / sum;
          //float planar = 2.0f*(w_1-w_2) / sum;
          //float sphere = 3.0f*(w_2     ) / sum;
          //ani_linear->SetTuple1( inPtId, linear);
          //ani_planar->SetTuple1( inPtId, planar);
          //isotrophy->SetTuple1( inPtId, sphere);

          float linear = (w[0]-w[1]) / w[0]+w[1];
          float planar = 2.0f*(w[1]) / w[0]+w[1];

          float sphere = 3.0f*( std::abs(w[2]) ) / std::abs(w[0]) + std::abs(w[1]) + std::abs(w[2]);

          ani_linear->SetTuple1( inPtId, linear);
          ani_planar->SetTuple1( inPtId, planar);

          isotrophy->SetTuple1( inPtId, sphere);
          
        }
      }
    }
  
  write_timer("EigenVectors", "calculate", timer.elapsed());

  vtkDebugMacro(<<"Generated " << numPts <<" Eigenvectors");

  //
  // Update output and release memory
  //

  output->ShallowCopy(input);

  vtkPointData *pd = output->GetPointData();
  pd->AddArray(ev0);
  pd->AddArray(ev1);
  pd->AddArray(ew0);
  pd->AddArray(ew1);

  if( b3DTensor || this->Output3EV )
  {
    pd->AddArray(ev2);
    pd->AddArray(ew2);

    if( this->AnisotrophyMeasure) {
      //pd->SetTCoords(anisotrophy);
      pd->AddArray(ani_linear);
      pd->AddArray(ani_planar);
      pd->AddArray(isotrophy);

      ani_linear->Delete();
      ani_planar->Delete();
      isotrophy->Delete();
    }
  }
  
  ev0->Delete();
  ev1->Delete();
  ev2->Delete();

  ew0->Delete();
  ew1->Delete();
  ew2->Delete();
  
  output->Squeeze();

  return 1;
}

//==============================================================================
void vtkEigenVectors::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
