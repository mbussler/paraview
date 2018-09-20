/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkWritePolyData.h

  Copyright (c) Michael Bußler

=========================================================================*/

#include "vtkWritePolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include "vtkPolyData.h"
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
#include <vtkXMLPolyDataWriter.h>

vtkStandardNewMacro(vtkWritePolyData);

inline std::string stripPath( const std::string& str )
{
    std::string s(str);
    s.erase(0, s.find_last_of("\\/") + 1);
    return s;
}

bool checkFileExists( const std::string& filename)
{
    std::string fnStr(filename);
    std::replace( fnStr.begin(), fnStr.end(), '\\', '/'); // replace \ with /
    std::ifstream ftest(fnStr.c_str(), std::ifstream::in);
    if( !ftest.good()) {
        ftest.close();
        //printf("Error finding file '%s'\n", fnStr.c_str());
        return false;
    }
    ftest.close();
    return true;
}

#define clampHigh(a,b) if(a>b) a=b


//==============================================================================
vtkWritePolyData::vtkWritePolyData()
{
  this->OutputFileName = NULL;
  
  this->SetNumberOfInputPorts(1);
  
}

int vtkWritePolyData::RequestUpdateExtent(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector)
{
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    std::string timeStr;
    if( outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())) {
        std::ostringstream os;
        os << outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
        timeStr = os.str();

        std::string df(this->OutputFileName);
        std::replace( df.begin(), df.end(), '\\', '/');

        std::ostringstream ss;
        ss << std::setw(4) << std::setfill('0') << timeStr;

        df += timeStr.empty() ? "" : "_"+ss.str();

        std::string filename = df + ".vtp";

        fstream test(filename.c_str());
        if( test.good() ){
            test.close();
            //vtkErrorMacro( << "file "<<filename<<" exists and will not be overwritten!")
            //return 0;
            //request->Remove(vtkStreamingDemandDrivenPipeline::REQUEST_DATA());
            outInfo->Set(vtkDemandDrivenPipeline::REQUEST_DATA_NOT_GENERATED());
        }
        test.close();
    }

    return 1;}

//==============================================================================
vtkWritePolyData::~vtkWritePolyData()
{
  if( this->OutputFileName) {
      delete[] OutputFileName;
  }
}

//==============================================================================
int vtkWritePolyData::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkIdType numPts;

  std::string timeStr;
  if( outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())) {
    std::ostringstream os;
    os << outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    timeStr = os.str();
  }

   numPts = input->GetNumberOfPoints();
 
   //if ( numPts < 1 )
   //  {
   //  vtkErrorMacro(<<"No data!");
   //  return 1;
   //  }

   std::string df(this->OutputFileName);
   std::replace( df.begin(), df.end(), '\\', '/');

   std::ostringstream ss;
   ss << std::setw(4) << std::setfill('0') << timeStr;

   df += timeStr.empty() ? "" : "_"+ss.str();

   std::string filename = df + ".vtp";

   fstream test(filename.c_str());
   if( test.good() ){
       test.close();
       vtkErrorMacro( << "file "<<filename<<" exists and will not be overwritten!")
       return 0;
   }
   test.close();


   // Write the file
   vtkSmartPointer<vtkXMLPolyDataWriter> writer =  
       vtkSmartPointer<vtkXMLPolyDataWriter>::New();
   writer->SetFileName(filename.c_str());
   writer->SetInputData(input);

   // Optional - set the mode. The default is binary.
   writer->SetDataModeToBinary();
   writer->EncodeAppendedDataOn();
   writer->SetCompressorTypeToZLib();

   writer->Write();
 
  output->ShallowCopy(input);

  output->Squeeze();

  return 1;
}

//==============================================================================
void vtkWritePolyData::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
