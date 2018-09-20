/*=========================================================================


  Program:   Visualization Toolkit
  Module:    vtkMichelLevy.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkMichelLevy.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkDataSet.h"
#include "vtkExecutive.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkMath.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkTransform.h"
#include "vtkSmartPointer.h"
#include "vtkMath.h"

#include "linalg.h"

vtkStandardNewMacro(vtkMichelLevy);

static float cie31[81][3] = {      //lambda =
  {0.001368, 0.000039, 0.006450},  //380
  {0.002236, 0.000064, 0.010550},  //385
  {0.004243, 0.000120, 0.020050},  //390
  {0.007650, 0.000217, 0.036210},  //395
  {0.014310, 0.000396, 0.067850},  //400
  {0.023190, 0.000640, 0.110200},  //405
  {0.043510, 0.001210, 0.207400},  //410
  {0.077630, 0.002180, 0.371300},  //415
  {0.134380, 0.004000, 0.645600},  //420
  {0.214770, 0.007300, 1.039050},  //425
  {0.283900, 0.011600, 1.385600},  //430
  {0.328500, 0.016840, 1.622960},  //435
  {0.348280, 0.023000, 1.747060},  //440
  {0.348060, 0.029800, 1.782600},  //445
  {0.336200, 0.038000, 1.772110},  //450
  {0.318700, 0.048000, 1.744100},  //455
  {0.290800, 0.060000, 1.669200},  //460
  {0.251100, 0.073900, 1.528100},  //465
  {0.195360, 0.090980, 1.287640},  //470
  {0.142100, 0.112600, 1.041900},  //475
  {0.095640, 0.139020, 0.812950},  //480
  {0.057950, 0.169300, 0.616200},  //485
  {0.032010, 0.208020, 0.465180},  //490
  {0.014700, 0.258600, 0.353300},  //495
  {0.004900, 0.323000, 0.272000},  //500
  {0.002400, 0.407300, 0.212300},  //505
  {0.009300, 0.503000, 0.158200},  //510
  {0.029100, 0.608200, 0.111700},  //515
  {0.063270, 0.710000, 0.078250},  //520
  {0.109600, 0.793200, 0.057250},  //525
  {0.165500, 0.862000, 0.042160},  //530
  {0.225750, 0.914850, 0.029840},  //535
  {0.290400, 0.954000, 0.020300},  //540
  {0.359700, 0.980300, 0.013400},  //545
  {0.433450, 0.994950, 0.008750},  //550
  {0.512050, 1.000000, 0.005750},  //555
  {0.594500, 0.995000, 0.003900},  //560
  {0.678400, 0.978600, 0.002750},  //565
  {0.762100, 0.952000, 0.002100},  //570
  {0.842500, 0.915400, 0.001800},  //575
  {0.916300, 0.870000, 0.001650},  //580
  {0.978600, 0.816300, 0.001400},  //585
  {1.026300, 0.757000, 0.001100},  //590
  {1.056700, 0.694900, 0.001000},  //595
  {1.062200, 0.631000, 0.000800},  //600
  {1.045600, 0.566800, 0.000600},  //605
  {1.002600, 0.503000, 0.000340},  //610
  {0.938400, 0.441200, 0.000240},  //615
  {0.854450, 0.381000, 0.000190},  //620
  {0.751400, 0.321000, 0.000100},  //625
  {0.642400, 0.265000, 0.000050},  //630
  {0.541900, 0.217000, 0.000030},  //635
  {0.447900, 0.175000, 0.000020},  //640
  {0.360800, 0.138200, 0.000010},  //645
  {0.283500, 0.107000, 0.000000},  //650
  {0.218700, 0.081600, 0.000000},  //655
  {0.164900, 0.061000, 0.000000},  //660
  {0.121200, 0.044580, 0.000000},  //665
  {0.087400, 0.032000, 0.000000},  //670
  {0.063600, 0.023200, 0.000000},  //675
  {0.046770, 0.017000, 0.000000},  //680
  {0.032900, 0.011920, 0.000000},  //685
  {0.022700, 0.008210, 0.000000},  //690
  {0.015840, 0.005723, 0.000000},  //695
  {0.011359, 0.004102, 0.000000},  //700
  {0.008111, 0.002929, 0.000000},  //705
  {0.005790, 0.002091, 0.000000},  //710
  {0.004109, 0.001484, 0.000000},  //715
  {0.002899, 0.001047, 0.000000},  //720
  {0.002049, 0.000740, 0.000000},  //725
  {0.001440, 0.000520, 0.000000},  //730
  {0.001000, 0.000361, 0.000000},  //735
  {0.000690, 0.000249, 0.000000},  //740
  {0.000476, 0.000172, 0.000000},  //745
  {0.000332, 0.000120, 0.000000},  //750
  {0.000235, 0.000085, 0.000000},  //755
  {0.000166, 0.000060, 0.000000},  //760
  {0.000117, 0.000042, 0.000000},  //765
  {0.000083, 0.000030, 0.000000},  //770
  {0.000059, 0.000021, 0.000000},  //775
  {0.000042, 0.000015, 0.000000}   //780
};

static double M_rgb[3][3] = {
  { 2.04159, -0.56501, -0.34473},
  {-0.96924,  1.87597,  0.04156},
  { 0.01344, -0.11836,  1.01517},
};

template <class T>
inline void swap( T& a, T& b) {
    T swp=a; a=b; b=swp;
};


//==============================================================================
vtkMichelLevy::vtkMichelLevy() 
: RGB_lin(0)
{
  this->SetNumberOfInputPorts(1);
  
  // by default, process active point tensors
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::TENSORS);
  
  this->AnalyzerAngle = 90.0;
  this->ObjectAngle = 0.0;
  this->Thickness = 1.0;
  this->FringeValue = 1.0;
  
  this->Wavelength = 380.0;
  this->Monochromatic = false;
  this->Lightfield = false;
  this->Intensity = 1.0;
  
  mat3setrows( this->adobe_rgb, M_rgb[0], M_rgb[1], M_rgb[2]);

  createBifringenceTable();
  
}
//==============================================================================
vtkMichelLevy::~vtkMichelLevy()
{
  // clean up
  if( RGB_lin ){
    delete RGB_lin;
    RGB_lin = 0;
  }

}

//==============================================================================
int vtkMichelLevy::RequestData(
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

  vtkDataArray *inTensor;
  double tensor[9];
  vtkIdType numPts, inPtId, ptIncr, i;
  int j;

  
  vtkDebugMacro(<<"Calculating MichelLevy");

  vtkPointData *outPD = output->GetPointData();
  inTensor = this->GetInputArrayToProcess(0, inputVector);
  numPts = input->GetNumberOfPoints();

  if ( !inTensor || numPts < 1 )
    {
    vtkErrorMacro(<<"No data!");
    return 1;
    }
  if (( inTensor->GetNumberOfComponents() != 4 ) &&
      ( inTensor->GetNumberOfComponents() != 9 ))
    {
    vtkErrorMacro("Input array must be a gradient tensor with 4 or 9 components.");
    return 0;
    }

  // allocate mem for output data arrays
  vtkSmartPointer<vtkUnsignedCharArray> colArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colArray->SetName("ColorRGB");
  colArray->SetNumberOfComponents(3);
  colArray->SetNumberOfTuples(numPts);

  vtkSmartPointer<vtkDoubleArray> retardArray = vtkSmartPointer<vtkDoubleArray>::New();
  retardArray->SetName("Retardation");
  retardArray->SetNumberOfComponents(1);
  retardArray->SetNumberOfTuples(numPts);

  vtkSmartPointer<vtkDoubleArray> ew0 = vtkSmartPointer<vtkDoubleArray>::New();
  ew0->SetName("ew0");
  ew0->SetNumberOfComponents(1);
  ew0->SetNumberOfTuples(numPts);
  
  vtkSmartPointer<vtkDoubleArray> ew1 = vtkSmartPointer<vtkDoubleArray>::New();
  ew1->SetName("ew1");
  ew1->SetNumberOfComponents(1);
  ew1->SetNumberOfTuples(numPts);
  
//   vtkSmartPointer<vtkDoubleArray> birefTable = vtkSmartPointer<vtkDoubleArray>::New();
//   birefTable->SetName("biref table");
//   birefTable->SetNumberOfComponents(3); //R-G-B
//   birefTable->SetNumberOfValues(6000/2);
//   
//   for( int i=0; i<6000/2; i++){
//     birefTable->SetTuple( i, RGB_lin[i]);
//   }
//   
  //
  // Traverse all Input points, calculate and store symmetric and antisymmetric tensors
  //

  for (inPtId=0; inPtId < numPts; inPtId++)
  {

      double sig1=0.0, sig2=0.0;
      inTensor->GetTuple(inPtId, tensor);

      if( inTensor->GetNumberOfComponents() == 4 )
      {
          mat2 m;

          for (j=0; j<2; j++)
            for (i=0; i<2; i++)
              m[i][j] = tensor[i+2*j];

          vec2 lambda;
          mat2eigenvalues(m, lambda);
      
          sig1 = std::max<double>(lambda[0], lambda[1]);
          sig2 = std::min<double>(lambda[0], lambda[1]);
      } 
      else 
      {
          mat3 m;

          for (j=0; j<3; j++)
              for (i=0; i<3; i++)
                  m[i][j] = tensor[j+3*i] * tensor[i+3*j];

          vec3 lambda;
          mat3eigenvalues(m, lambda);
          sig1 = std::max<double>( lambda[0], std::max<double>( lambda[1], lambda[2]));
          sig2 = std::min<double>( lambda[0], std::min<double>( lambda[1], lambda[2]));
      }

      ew0->SetValue(inPtId, sig1);
      ew1->SetValue(inPtId, sig2);      

      //vec2 lEV;
      //mat2realEigenvector(m, sigma1, lEV);

      // polarizers priviledged direction 
      //vec2 x_axis = {1.0, 0.0};

      // calculate angle between pol. direction and objects closest priv. direction
      //double alpha = acos( vec2dot(lEV, x_axis));

      // choose closest priv. direction
      //double pi_2 = vtkMath::Pi()/2.0;
      //while( alpha > pi_2 )
      //  alpha -= pi_2;
      //double pi_2m = -1.0*pi_2;
      //while( alpha < pi_2m )
      //  alpha += pi_2;

      //double delta = sigma1 - sigma2;
      double delta = 2.0 * vtkMath::Pi() * this->Thickness * this->FringeValue * ( sig1 - sig2);

      retardArray->SetValue(inPtId, delta);

      vec3 colors = {0.0,0.0,0.0};
      //getRGB(delta, alpha, colors);
      getRGB( delta, colors);

      unsigned char r = (colors[0]>1.0) ? 255 : (colors[0]<0.0) ? 0 : (unsigned char) (colors[0]*255.0 + 0.5);
      unsigned char g = (colors[1]>1.0) ? 255 : (colors[1]<0.0) ? 0 : (unsigned char) (colors[1]*255.0 + 0.5);
      unsigned char b = (colors[2]>1.0) ? 255 : (colors[2]<0.0) ? 0 : (unsigned char) (colors[2]*255.0 + 0.5);

      colArray->SetTuple3(inPtId, r, g, b );
  }
 

  vtkDebugMacro(<<"Calculated " << numPts <<" values");

  //
  // Update output and release memory
  //

  output->ShallowCopy(input);

  vtkPointData *pd = output->GetPointData();
  pd->AddArray(colArray);
  pd->AddArray(retardArray);
  pd->AddArray(ew0);
  pd->AddArray(ew1);
//   pd->AddArray(birefTable);
 
  output->Squeeze();

  return 1;
}

//==============================================================================
void vtkMichelLevy::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void vtkMichelLevy::getRGB( const double& retardation, const double& alpha, vec3& rgb )
{

    // XYZ color coordinates
    vec3 xyz = {0.0, 0.0, 0.0};

    // integrate CIE-colors
    for( int i = 0; i < 81; i++)
    {
      double lambda = 380.0 + 5*i;
      xyz[0] += Trans( retardation, lambda, alpha )*cie31[i][0] ; //* 5.0; // don't required because of normalization
      xyz[1] += Trans( retardation, lambda, alpha )*cie31[i][1] ; //* 5.0;
      xyz[2] += Trans( retardation, lambda, alpha )*cie31[i][2] ; //* 5.0;
    }

    // normalize:
    float norm = xyz[0] + xyz[1] + xyz[2];
    if( norm > 0){
      xyz[0] /= norm;
      xyz[1] /= norm;
      xyz[2] /= norm;
    }
    
    // make 'em RGB...
    mat3vec( adobe_rgb, xyz, rgb);
  
};



void vtkMichelLevy::createBifringenceTable(){

  if( RGB_lin ){
    delete[] RGB_lin;
    RGB_lin = 0;
  }

  int num_values = 6000;
  int step = 2;
  
  // calculate tristimulus values from spectral power distribution
  // for different birefringence values
  RGB_lin = new vec3[num_values/step];
  
  for( int d=0; d< (num_values/step); d++ ){

    vec3 xyz = {0.0, 0.0, 0.0};

    for( int i = 0; i < 81; i++)
    {
      double lambda = 380.0 + 5*i;
      xyz[0] += TransFixed( d*step, lambda )*cie31[i][0] ; //* 5.0; // don't required because of normalization
      xyz[1] += TransFixed( d*step, lambda )*cie31[i][1] ; //* 5.0;
      xyz[2] += TransFixed( d*step, lambda )*cie31[i][2] ; //* 5.0;
    }

    // normalize:
    float norm = xyz[0] + xyz[1] + xyz[2];
    if( norm > 0){
      xyz[0] /= norm;
      xyz[1] /= norm;
      xyz[2] /= norm;
    }
    
    // make 'em RGB...
    mat3vec( adobe_rgb, xyz, RGB_lin[d]);
  }
};

inline double vtkMichelLevy::TransFixed( double Gamma, double lambda) {
  return pow( sin(vtkMath::Pi()*Gamma/lambda ), 2.0);
};

inline double vtkMichelLevy::Trans( double Gamma, double lambda, double alpha) {
  double fixedTerm = TransFixed(Gamma, lambda);
  double phi = this->AnalyzerAngle * vtkMath::Pi() / 180.0;
  double tau = alpha + this->ObjectAngle * vtkMath::Pi() / 180.0;
  return pow(cos(phi),2.0) - pow(sin(tau-phi),2.0) * pow(sin(tau),2.0) * fixedTerm;
};

inline double vtkMichelLevy::T( double delta, double lambda) {
  return pow(sin( delta/(2.0*lambda)),2.0);
};

inline double vtkMichelLevy::TL( double delta, double lambda) {
  return pow(cos( delta/(2.0*lambda)),2.0);
};

void vtkMichelLevy::getRGB( const double& retardation, vec3& rgb)
{
  
    vec3 xyz = {0.0, 0.0, 0.0};

    if( this->Monochromatic ) 
    {
      int i = (this->Wavelength - 380.0 + 2.5) / (int)5; // calculate index from given wave length by rounding and tuncation
      if( this->Lightfield )
      {
        xyz[0] = T( retardation, this->Wavelength )*cie31[i][0] ; 
        xyz[1] = T( retardation, this->Wavelength )*cie31[i][1] ; 
        xyz[2] = T( retardation, this->Wavelength )*cie31[i][2] ;                 
      }
      else
      {
        xyz[0] = TL( retardation, this->Wavelength )*cie31[i][0] ; 
        xyz[1] = TL( retardation, this->Wavelength )*cie31[i][1] ; 
        xyz[2] = TL( retardation, this->Wavelength )*cie31[i][2] ;         
      }   
      vec3scal( xyz, this->Intensity, xyz);
    } 
    else // colored
    {
      if( this->Lightfield )
      {
        // integrate wavelengths of visible light
        for( int i = 0; i < 81; i++)
        {
          double lambda = 380.0 + 5*i;
          xyz[0] += TL( retardation, lambda )*cie31[i][0] ; //* 5.0; // don't required because of normalization
          xyz[1] += TL( retardation, lambda )*cie31[i][1] ; //* 5.0;
          xyz[2] += TL( retardation, lambda )*cie31[i][2] ; //* 5.0;
        }
      }
      else
      {
        // integrate wavelengths of visible light
        for( int i = 0; i < 81; i++)
        {
          double lambda = 380.0 + 5*i;
          xyz[0] += T( retardation, lambda )*cie31[i][0] ; //* 5.0; // don't required because of normalization
          xyz[1] += T( retardation, lambda )*cie31[i][1] ; //* 5.0;
          xyz[2] += T( retardation, lambda )*cie31[i][2] ; //* 5.0;
        }
      }

      vec3scal(xyz, this->Intensity/81.0, xyz);
      }

    // normalize
    if(vec3mag(xyz) > this->Gamma) {
        vec3nrm(xyz, xyz);
        vec3scal(xyz, this->Gamma, xyz );
    }

    // make 'em RGB...
    mat3vec( adobe_rgb, xyz, rgb);
  
}
