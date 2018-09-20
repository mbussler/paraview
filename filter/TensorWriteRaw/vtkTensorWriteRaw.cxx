/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTensorWriteRaw.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkTensorWriteRaw.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkImageData.h"
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
#include "vtkSmartPointer.h"
#include <sstream>
#include <algorithm>

vtkStandardNewMacro(vtkTensorWriteRaw);

inline std::string stripPath( const std::string& str )
{
    std::string s(str);
    s.erase(0, s.find_last_of("\\/") + 1);
    return s;
}

#define clampHigh(a,b) if(a>b) a=b


//==============================================================================
vtkTensorWriteRaw::vtkTensorWriteRaw()
{
  this->DataFileName = NULL;
  
  this->SetNumberOfInputPorts(1);
  
  // by default, process active point tensors
  this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::TENSORS);

}
//==============================================================================
vtkTensorWriteRaw::~vtkTensorWriteRaw()
{
  if( this->DataFileName) {
      delete[] DataFileName;
  }
}

//==============================================================================
int vtkTensorWriteRaw::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkImageData *imageData = vtkImageData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkImageData *output = vtkImageData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDataArray *inData;
  double data[9];
  vtkIdType numPts;
  int j;

  std::string timeStr;
  if( outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())) {
    std::ostringstream os;
    os << outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    timeStr = os.str();
  }

   vtkPointData *outPD = output->GetPointData();
   inData = this->GetInputArrayToProcess(0, inputVector);
   numPts = imageData->GetNumberOfPoints();
 
   if ( !inData || numPts < 1 )
     {
     vtkErrorMacro(<<"No data!");
     return 1;
     }

   int numComponents = inData->GetNumberOfComponents();
   bool bTens6 = numComponents == 6;
   bool bTens9 = numComponents == 9;
   if( numComponents > 6) 
       numComponents = 6;

   //if( !( bTens6 || bTens9 ))
   //  {
   //  vtkErrorMacro("Input array must be a tensor with 6 or 9 components.");
   //  return 0;
   //  }

   const int tensIdx6[6] = { 0, 1, 2, 3, 4, 5 };
   const int tensIdx9[6] = { 0, 4, 8, 1, 5, 2 }; // XX, YY, ZZ, XY, YZ, XZ

   const int* off = bTens9 ? tensIdx9 : tensIdx6;
   
   std::string df(this->DataFileName);
   std::replace( df.begin(), df.end(), '\\', '/');

   std::ostringstream ss;
   ss << std::setw(4) << std::setfill('0') << timeStr;

   df += timeStr.empty() ? "" : "_"+ss.str();

   std::string filename = df + ".raw";

   fstream test(filename.c_str());
   if( test.good() ){
       test.close();
       vtkErrorMacro( << "file "<<filename<<" exists and will not be overwritten!")
       return 0;
   }
   test.close();

   ofstream outfile( filename.c_str(), ios::out | ios::binary );
  
   if( !outfile.is_open()) {
     outfile.close();
     vtkErrorMacro("Could not open file for writing.");
     return 0;
   }
     
   std::cout << "Writing " << inData->GetName() << " to " << filename << std::endl;
     
   
   //
   // Traverse all Input points, calculate and store symmetric and antisymmetric tensors
   //
 
  int* dims = imageData->GetDimensions();
   
   
  std::cout << "Dims: " << " x: " << dims[0] << " y: " << dims[1] << " z: " << dims[2] << std::endl;
 
  std::cout << "Number of points: " << imageData->GetNumberOfPoints() << std::endl;
  std::cout << "Number of cells: " << imageData->GetNumberOfCells() << std::endl;
  
  
  // Retrieve the entries from the image data and print them to the screen
  size_t nBytes = 0;
  for ( int z=0 ; z < dims[2]; z++)
    {
    for ( int y=0 ; y < dims[1]; y++)
      {
      for ( int x=0 ; x < dims[0]; x++)
        {
          vtkIdType ptId = x + y*dims[0] + z*dims[0]*dims[1];
          inData->GetTuple( ptId, data);

          float dataOut[6];
          for (j=0; j<numComponents; j++) {
              dataOut[j] = (float) data[ off[j]];
          }

          // write to file
          size_t size = numComponents*sizeof(float);
          outfile.write((char*)dataOut, size);
          nBytes += size;
        }
      }
    }

  outfile.flush();
  outfile.close();

  std::cout << nBytes << " bytes written to " << std::string(this->DataFileName) << ".dat" << std::endl;

  double spacing[3], origin[3];
  imageData->GetSpacing(spacing);
  imageData->GetOrigin(origin);

  filename = df + ".dat";
  outfile.open( filename.c_str(), ios::out );
  std::string rawFileName = stripPath(df) + ".raw";
  if( outfile.is_open()) {
    outfile << "ObjectFileName: " << rawFileName        <<std::endl;

    outfile << "Format:         FLOAT"                  <<std::endl;
    outfile << "GridType:       EQUIDISTANT"            <<std::endl;
    outfile << "Components:     " << numComponents      <<std::endl;
    outfile << "Dimensions:     3"                      <<std::endl;
    outfile << "TimeSteps:      1"                      <<std::endl;
    outfile << "ByteOrder:      LITTLE_ENDIAN"          <<std::endl;
    outfile << "Resolution:     " << dims[0] << " " 
                                  << dims[1] << " " 
                                  << dims[2]            << std::endl;
    outfile << "SliceThickness: " << spacing[0] << " " 
                                  << spacing[1] << " " 
                                  << spacing[2]         <<std::endl;
    outfile << "Origin:         " << origin[0] << " " 
                                  << origin[1] << " " 
                                  << origin[2]          <<std::endl;
    outfile << "Time:           " << timeStr            <<std::endl;

    outfile.close();
  }


 
  output->ShallowCopy(imageData);

//   vtkPointData *pd = output->GetPointData();
//   pd->AddArray(normArray);
 
  output->Squeeze();

  return 1;
}

//==============================================================================
void vtkTensorWriteRaw::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
