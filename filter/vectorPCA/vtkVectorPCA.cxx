/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkVectorPCA.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkVectorPCA.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkCell.h"
#include "vtkCellData.h"
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
#include "vtkSmartPointer.h"
#include "vtkTriangle.h"

#include "linalg.h"

#include <stdlib.h>

inline bool SameSign(float a, float b) {
    return a*b >= 0.0f;
}

inline void swap( int&a, int&b ){
  int swp=a; a=b; b=swp;
};

vtkStandardNewMacro(vtkVectorPCA);

//==============================================================================
vtkVectorPCA::vtkVectorPCA()
{
  this->SetNumberOfInputPorts(1);
  
  // by default, process active point tensors
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::VECTORS);

}
//==============================================================================
vtkVectorPCA::~vtkVectorPCA()
{
  // clean up
}

//==============================================================================
int vtkVectorPCA::RequestData(
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

  vtkIdType cellId, ptId, faceId, facePtId;
  vtkIdType numCells, numPts, numFaces;
  vtkPointData *inPD=input->GetPointData();
  vtkCellData *outCD=output->GetCellData();
  int maxCellSize=input->GetMaxCellSize();
  vtkIdList *cellPts, *facePts;
  vtkCell *cell, *face;
  vtkDataArray *inVelocities;
  double vecIn[3], faceVel[3], fcoords[3], x[3], lEV[3];
  double *weights;
  int subId;
  double C[9];

  double pointCoordsA[3], pointCoordsB[3], AC[3];
  double base, height, area, cellFlux;
  
  double *m[3], w[3], *v[3];
  double m0[3], m1[3], m2[3];
  double v0[3], v1[3], v2[3];
  double xv[3], yv[3], zv[3];

  // set up working matrices
  m[0] = m0; m[1] = m1; m[2] = m2;
  v[0] = v0; v[1] = v1; v[2] = v2;

  float m11,m12,m21,m22;

  if ( (numCells=input->GetNumberOfCells()) < 1 )
    {
    vtkDebugMacro(<<"No input cell data!");
    return 1;
    }

  inVelocities = this->GetInputArrayToProcess(0, inputVector);
  numPts = input->GetNumberOfPoints();

  if ( !inVelocities || numPts < 1 )
  {
    vtkErrorMacro(<<"No data!");
    return 1;
  }

  vtkSmartPointer<vtkDoubleArray> mainComponent = vtkSmartPointer<vtkDoubleArray>::New();
  mainComponent->SetName("vec_pca");
  mainComponent->SetNumberOfComponents(3);
  mainComponent->SetNumberOfTuples(numCells);
  
  vtkDebugMacro(<<"Calculating main component via PCA");

  int abort=0;
  vtkIdType progressInterval = numCells/10 + 1;
  int hasEmptyCells = 0;
  for (cellId=0; cellId < numCells && !abort; cellId++)
  {
    if ( ! (cellId % progressInterval) )
    {
      vtkDebugMacro(<<"Processing #" << cellId);
      this->UpdateProgress (0.5*cellId/numCells);
      abort = this->GetAbortExecute();
    }
    
    cell = input->GetCell(cellId);

    if (cell->GetCellType() != VTK_EMPTY_CELL)
    {  
      // do PCA of cell velocities and store in mainComponent array
      vtkIdType numNodes = cell->GetPointIds()->GetNumberOfIds();
      double* evSet[3];
      evSet[0] = new double(numNodes*2);
      evSet[1] = new double(numNodes*2);
      evSet[2] = new double(numNodes*2);

      for( vtkIdType cellPtId = 0; cellPtId < numNodes; cellPtId++ )
      {
          vtkIdType ptId = cell->GetPointId(cellPtId);
          inVelocities->GetTuple( ptId, vecIn);
          evSet[0][2*cellPtId+0] =  vecIn[0];
          evSet[1][2*cellPtId+0] =  vecIn[1];
          evSet[2][2*cellPtId+0] =  vecIn[2];
          evSet[0][2*cellPtId+1] = -vecIn[0];
          evSet[1][2*cellPtId+1] = -vecIn[1];
          evSet[2][2*cellPtId+1] = -vecIn[2];
      }        
      
      //Covariance Matrix
      mat3 m;
      for(int j=0; j<3; ++j) {
        for(int i=0; i<3; ++i) {
            double sum = 0.0;
            for(int n=0; n<2*numNodes; ++n){
              sum += evSet[i][n] * evSet[j][n];
            }
            m[i][j] = sum/2.0 * (double)numNodes;
        }
      }

      vec3 lambda;
      mat3eigenvalues(m, lambda);

      mat3 ev;
      mat3realEigenvector(m, lambda[0], ev[0]);
      mat3realEigenvector(m, lambda[1], ev[1]);
      mat3realEigenvector(m, lambda[2], ev[2]);
      
      // find largest EW
      int n=3;
      int sorted[3];
      for( int i=0; i<n; i++) {
        sorted[i] = i;
      }
      
      // sort by descending order
      bool swapped=true;
      while(swapped) 
      {
        swapped=false;
        for( int i=0; i<(n-1); i++) {
          if( lambda[sorted[i]] < lambda[sorted[i+1]]) {
            swap( sorted[i], sorted[i+1]);
            swapped = true;
          }
        }
      }
      
      lEV[0] = ev[sorted[0]][0];
      lEV[1] = ev[sorted[0]][1];
      lEV[2] = ev[sorted[0]][2];

      delete[] evSet[0];
      delete[] evSet[1];
      delete[] evSet[2];

    }
    else
    {
      hasEmptyCells = 1;
    }
    mainComponent->SetTuple(cellId, lEV);
  }  

  vtkDebugMacro(<<"Calculated " << numPts <<" VectorPCA values");

  //
  // Update output and release memory
  //

  output->ShallowCopy(input);

  outCD->AddArray(mainComponent);
  
  output->Squeeze();

  return 1;
}

//==============================================================================
void vtkVectorPCA::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void vtkVectorPCA::GetInvariants(double* m, double* pqr) 
{
    pqr[0] = -m[0*3+0]*m[1*3+1]*m[2*3+2]
    -m[0*3+1]*m[1*3+2]*m[2*3+0]
    -m[0*3+2]*m[1*3+0]*m[2*3+1]
    +m[2*3+0]*m[1*3+1]*m[0*3+2]
    +m[2*3+1]*m[1*3+2]*m[0*3+0]
    +m[2*3+2]*m[1*3+0]*m[0*3+1];


    pqr[1] = m[1*3+1]*m[2*3+2] - m[2*3+1]*m[1*3+2]
    + m[0*3+0]*m[1*3+1] - m[1*3+0]*m[0*3+1]
    + m[0*3+0]*m[2*3+2] - m[0*3+2]*m[2*3+0];

    pqr[2] = -(m[0*3+0]+m[1*3+1]+m[2*3+2]);
}


int vtkVectorPCA::CalculateRoots(double* a, double* r) {
    double c1 = a[1] - a[2]*a[2]/3.;
    double c0 = a[0] - a[1]*a[2]/3. + 2./27.*a[2]*a[2]*a[2];
    // Make cubic coefficient 4 and linear coefficient +- 3
    // by substituting y = z*k and multiplying with 4/k^3
    if (c1 == 0) {
        if (c0 == 0)     r[0] =  0;
        else if (c0 > 0) r[0] = -pow(c0, 1./3.);
        else             r[0] =  pow(-c0, 1./3.);
    }
    else {
        bool negc1 = c1 < 0;
        double absc1 = negc1 ? -c1 : c1;
        double k = sqrt(4./3.*absc1);
        double d0 = c0 * 4./(k*k*k);
        // Find the first solution
        if (negc1) {
            if (d0 > 1)       r[0] = -cosh(acosh(d0)/3);
            else if (d0 > -1) r[0] = -cos(acos(d0)/3);
            else              r[0] =  cosh(acosh(-d0)/3);
        }
        else {
            r[0] = -sinh(asinh(d0)/3);
        }
        // Transform back
        r[0] *= k;
    }
    r[0] -= a[2]/3;
    // Other two solutions
    double p = r[0] + a[2];
    double q = r[0]*p + a[1];
    double discrim = p*p - 4*q;
    //    if (forceReal && discrim < 0.0) discrim = 0.0;
    if (discrim >= 0) {
        double root = sqrt(discrim);
        r[1] = (-p - root)/2.;
        r[2] = (-p + root)/2.;
        return 3;
    }
    else {
        double root = sqrt(-discrim);
        r[1] = -p/2;
        r[2] = root/2.;
        return 1;
    }    
}


void vtkVectorPCA::CalculateEigenvector(double* m, double lambda, double* evect) {
    double norm[3];
    double cross[9];
    double red[9];

    red[0] =  m[0]-lambda;
    red[1] =  m[1];
    red[2] =  m[2];
    red[3] =  m[3];
    red[4] =  m[4]-lambda;
    red[5] =  m[5];
    red[6] =  m[6];
    red[7] =  m[7];
    red[8] =  m[8]-lambda;

    CalculateCross(&red[3], &red[6], &cross[0]);
    CalculateCross(&red[6], &red[0], &cross[3]);
    CalculateCross(&red[0], &red[3], &cross[6]);

    norm[0] = NormSquared(&cross[0]);
    norm[1] = NormSquared(&cross[3]);
    norm[2] = NormSquared(&cross[6]);
    
    int best = getLargest(norm);
    double len = sqrt(norm[best]);

    if (len > 0) {
        evect[0] = cross[best*3+0] / len;
        evect[1] = cross[best*3+1] / len;
        evect[2] = cross[best*3+2] / len;
    }
    else {
        evect[0] = 0.0;
        evect[1] = 0.0;
        evect[2] = 0.0;
    }
}


void vtkVectorPCA::CalculateCross(double* v1, double* v2, double* res) {
    res[0] = v1[1]*v2[2] - v1[2]*v2[1];
    res[1] = v1[2]*v2[0] - v1[0]*v2[2];
    res[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

double vtkVectorPCA::NormSquared(double* v) {
    return (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}


int vtkVectorPCA::getLargest(double* v) {
    if (v[0]>v[1]) {
    if (v[0]>v[2]) {
        return 0;
    }
    else {
        return 2;
    }
    }
    else {
    if (v[1]>v[2]) {
        return 1;
    }
    else {
        return 2;
    }
    }
}
