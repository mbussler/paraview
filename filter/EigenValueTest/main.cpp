
#include <stdio.h>
//#include <vtkMath.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "linalg.h"

using Eigen::Matrix2f;
using Eigen::Vector2f;

void GetInvariants(double* m, double* pqr);
int CalculateRoots(double* a, double* r);
void CalculateEigenvector(double* m, double lambda, double* evect);
void CalculateCross(double* v1, double* v2, double* res);
double NormSquared(double* v);
int getLargest(double* v);

inline bool SameSign(float a, float b) {
    return a*b >= 0.0f;
}

// main
int main(int argc, const char* argv[]) 
{

    float v0[2],v1[2];
    float* m[2];
    m[0] = v0; m[1] = v1;
    
    float ev0[2], ev1[2];
    float* v[2];
    v[0]=ev0; v[1]=ev1;
    float w[2];
    
    float m11,m12,m21,m22;
    
//    printf("\n\nParaview\n\n");
//
//    for( int y=3; y>=-3; y-- )
//    {
//        for( int x=-3; x<=3; x++ )
//        {
//               m11 = m[0][0] = 1.0-2.0*x;
////             m[0][0] = 1.0+2.0*x/3.0;
//               m12 = m[0][1] = y;
//               m21 = m[1][0] = y;
//               m22 = m[1][1] = 1;
//            
//            printf("[%d %d]\t", x, y);
//            printf("[%.0f %.0f %.0f %.0f]\t", m[0][0], m[0][1], m[1][0], m[1][1]);
//
//            vtkMath::JacobiN(m, 2, w, v);
//
//            for( int n=0; n<2; n++) 
//            {
//              float c11 = m11-w[n];
//              float c12 = m12;
//              float c21 = m21;
//              float c22 = m22-w[n];
//              
//              // just to be sure. 
//              assert( SameSign(c11, c12) ==  SameSign(c21, c22));
//
//              //  // now check EVs and flip component if necessary
//              //  float ev1 = v[n][0];
//              //  float ev2 = v[n][1];
//              //  
//              //  if(  (  SameSign(c11, c12) &&  SameSign(ev1, ev2)) ||
//              //       ( !SameSign(c11, c12) && !SameSign(ev1, ev2))  )
//              //  {
//              //      v[n][0] *=-1;
//              //  }
//            }
//            
//            printf("EW: %.4f EV: [%.4f %.4f]\t", w[0], ev0[0], ev0[1]);
//            printf("EW: %.4f EV: [%.4f %.4f]\n", w[1], ev1[0], ev1[1]);
//
//        }
//        printf("\n");
//    }
//    

    printf("\n\nEigen\n\n");

    for( int y=3; y>=-3; y-- )
    {
        for( int x=-3; x<=3; x++ )
        {
            Matrix2f m(2,2);
            m(0,0) = 1.0f;
            //m(0,0) = 1.0+2.0*x/3.0;
            m(0,1) = m(1,0) = 0.0f;
            m(1,1) = 1.0f;
            
            printf("[%d %d]\t", x, y);
            printf("[%.2f %.2f %.2f %.2f]\t", m(0,0), m(0,1), m(1,0), m(1,1));

            Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> es (m, Eigen::ComputeEigenvectors);
            
            Vector2f ew = es.eigenvalues();
            Matrix2f ev = es.eigenvectors();

            for( int n=0; n<2; n++) 
            {
              float c11 = m(0,0)-ew(n);
              float c12 = m(0,1);
              float c21 = m(0,1);
              float c22 = m(1,1)-ew(n);
              
              // just to be sure. 
              assert( SameSign(c11, c12) ==  SameSign(c21, c22));

            //  // now check EVs and flip component if necessary
            //  float ev1 = ev(n,0);
            //  float ev2 = ev(n,1);
            //  
            //  if(  (  SameSign(c11, c12) &&  SameSign(ev1, ev2)) ||
            //       ( !SameSign(c11, c12) && !SameSign(ev1, ev2))  )
            //  {
            //      ev(n,0) *=-1;
            //  }
            }
            
            int max = 0, min=1;
            if( ew(min) > ew(max))
            {
              max=1; min=0;
            }
            
            //std::cout << "[" << ev << "] " << ew << std::endl;  
            printf("EW: %.2f EV: [%.2f %.2f]\t", ew(max), ev(max,0), ev(max,1));
            printf("EW: %.2f EV: [%.2f %.2f]\n", ew(min), ev(min,0), ev(min,1));
            
        }
        printf("\n");
    }

    printf("\n\nLinalg\n\n");

    for( int y=3; y>=-3; y-- )
    {
        for( int x=-3; x<=3; x++ )
        {
            mat2 m;
            m[0][0] = 1.0; // trisec
            //m[0][0] = 1.0+2.0*x/3.0; // double wedge
            m[0][1] = 0.0f;
            m[1][0] = 0.0f;
            m[1][1] = 1.0f;
            
            printf("[%d %d]\t", x, y);
            printf("[%.2f %.2f %.2f %.2f]\t", m[0][0], m[0][1], m[1][0], m[1][1]);

            vec2 lambda;
            mat2eigenvalues(m, lambda);
            
            mat2 ev;
            mat2realEigenvector(m, lambda[0], ev[0]);
            mat2realEigenvector(m, lambda[1], ev[1]);
            
            int max = 0, min=1;
            if( lambda[max] < lambda[min])
            {
              max=1; min=0;
            }
            
            //std::cout << "[" << ev << "] " << ew << std::endl;  
            printf("EW: %.2f EV: [%.2f %.2f]\t", lambda[max], ev[max][0], ev[max][1]);
            printf("EW: %.2f EV: [%.2f %.2f]\n", lambda[min], ev[min][0], ev[min][1]);
            
        }
        printf("\n");
    }
    return 0;
}

