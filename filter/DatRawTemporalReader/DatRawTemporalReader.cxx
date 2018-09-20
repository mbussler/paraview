#include "DatRawTemporalReader.h"

#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkImageData.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "datRaw.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <functional>
#include <direct.h>

enum TensorF {
    XX_YY_ZZ_XY_XZ_YZ = 0,
    XX_YY_ZZ_XY_YZ_XZ = 1,
};

const int tensF_0[9] = { 0, 3, 4, 
    3, 1, 5, 
    4, 5, 2 };
const int tensF_1[9] = { 0, 3, 5, 
    3, 1, 4, 
    5, 4, 2 };


inline std::string getPath( const std::string& str )
{
    std::string s(str);
    std::replace( s.begin(), s.end(), '\\', '/'); // replace \ with /
    return s.substr(0, s.find_last_of("/"));
}

vtkStandardNewMacro(DatRawTemporalReader);

DatRawTemporalReader::DatRawTemporalReader()
{
    this->FileName = NULL;
    this->SetNumberOfInputPorts(0);
    this->SetNumberOfOutputPorts(1);
    this->m_NumberOfTimeSteps=0;
    this->bRelativePath = false;
    this->TensorFormat = 0;
}

void DatRawTemporalReader::ReadDatList()
{
    std::ifstream dltFile;
    const int lineSize = 1024;
    char line[lineSize];
    int numEntries = 0;

    dltFile.open(FileName, std::ifstream::in);
    std::string dltFilePath = getPath( this->FileName );

    m_DatList.clear();
    m_TimestepValues.clear();
    m_NumberOfTimeSteps = 0;

    while (dltFile.good()) {

        dltFile.getline(line, 1024);

        // Avoid reading empty lines
        if (strlen(line) > 0) {
            m_DatList.push_back(line);
            numEntries++;
        }

        //-read time step value---
        DatRawFileInfo info;

        // check if file entry is relative to dlt file
        std::ifstream datFile( line );
        if( !datFile.good()) {
            datFile.close();
            datFile.open( std::string (dltFilePath + '/' + line).c_str());
            if( datFile.good() ) {
                this->bRelativePath = true;
                m_DatList.pop_back();
                m_DatList.push_back(dltFilePath + '/' + line);
            }
        }
        datFile.close();


        std::string datFileName = m_DatList.back();
        if (!datRaw_readHeader(datFileName.c_str(), &info, 0)) {
            vtkErrorMacro("Could not read datRaw header from " + datFileName);
            return;
        }

        if( info.time == 0.0 && m_TimestepValues.size() > 0 && m_TimestepValues[0] >= 0.0 )
            m_TimestepValues.push_back(m_TimestepValues.size()); // use file index as time index
        else
            m_TimestepValues.push_back(info.time);

        datRaw_close(&info);
        datRaw_freeInfo(&info);

        //------------------------
    }
    m_NumberOfTimeSteps = numEntries;

    dltFile.close();
}

int DatRawTemporalReader::RequestInformation(vtkInformation*,
    vtkInformationVector**,
    vtkInformationVector* outputVector)
{
    vtkInformation* outInfo = outputVector->GetInformationObject(0);

    if( this->FileName == "") {
        vtkDebugMacro("FileName not set");
        return 0;
    }

    ReadDatList();

    if( m_DatList.size() == 0)
        return 0;

    DatRawFileInfo info;
    std::string datFileName = m_DatList[0];
    if (!datRaw_readHeader(datFileName.c_str(), &info, 0)) {
        vtkErrorMacro("Could not read datRaw header from " + datFileName);
        return 0;
    }

    Extent[0] = 0;
    Extent[1] = info.resolution[0] - 1;
    Extent[2] = 0;
    Extent[3] = info.resolution[1] - 1;
    Extent[4] = 0;
    Extent[5] = info.resolution[2] - 1;

    NumberOfComponents = info.numComponents;

    Spacing[0] = info.sliceDist[0];
    Spacing[1] = info.sliceDist[1];
    Spacing[2] = info.sliceDist[2];

    Origin[0] = info.origin[0];
    Origin[1] = info.origin[1];
    Origin[2] = info.origin[2];

    Time = info.time;

    datRaw_close(&info);
    datRaw_freeInfo(&info);

    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), Extent, 6);
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), 
        &this->m_TimestepValues[0], this->m_TimestepValues.size());
    double timeRange[2] = {m_TimestepValues[0], m_TimestepValues[this->m_TimestepValues.size()-1]};
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),
        timeRange, 2);
    outInfo->Set(vtkDataObject::SPACING(), Spacing, 3);
    outInfo->Set(vtkDataObject::ORIGIN(), Origin, 3);

    return 1;
}

// struct WithinTolerance : public std::binary_function<double,double,bool> {
//   bool operator() (double a, double b) const
//   {  
//     if (std::abs(a-b) < 0.00001)
//       return true;
//     else
//       return false;
//   }
// };

int DatRawTemporalReader::FindClosestTimeStep(double requestedTimeValue)
{
    int ts = 0;
    double mindist = std::abs(m_TimestepValues[0] - requestedTimeValue);

    for (int i = 0; i < m_TimestepValues.size(); i++) {

        double tsv = m_TimestepValues[i];
        double dist = std::abs(tsv - requestedTimeValue);
        if (dist < mindist) {
            mindist = dist;
            ts = i;
        }
    }
    return ts;
}

int DatRawTemporalReader::RequestData(vtkInformation *vtkNotUsed(request),
    vtkInformationVector **vtkNotUsed(inputVector),
    vtkInformationVector *outputVector)
{

    // get the info object
    vtkInformation *outInfo = outputVector->GetInformationObject(0);

    vtkImageData *output = vtkImageData::
        SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    if (output == 0) {
        vtkErrorMacro("Error while creating reader output.");
    }

    int timestepToLoad = 0;
    if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())
        && m_TimestepValues.size() > 0 ) {

            double requestedTimeValue = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

            timestepToLoad = FindClosestTimeStep(requestedTimeValue);
            // std::cout << "timestepToLoad = " << timestepToLoad << std::endl;

            // double ActualTimeStep = vtkstd::
            //   find_if(this->TimestepValues.begin(), this->TimestepValues.end(),
            // 	      vtkstd::bind2nd(WithinTolerance(), requestedTimeValue)) - this->TimestepValues.begin();

            // std::cout << "requestedTimeValue = " << requestedTimeValue << std::endl;
            // std::cout << "ActualTimeStep = " << ActualTimeStep << std::endl;
            // std::cout << "ActualTimeValue = " << TimestepValues[ActualTimeStep] << std::endl;

            // output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), requestedTimeValue);
            // timestepToLoad = ActualTimeStep;

            output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), m_TimestepValues[timestepToLoad]);
    }

    //---------------------------------------------------------------------------
    DatRawFileInfo info;
    std::string datFileName = m_DatList[timestepToLoad];
    if (!datRaw_readHeader(datFileName.c_str(), &info, 0)) {
        vtkErrorMacro("Could not read datRaw header");
        return 0;
    }

    float *dataBuffer = 0;

    std::string path = getPath( datFileName );
    chdir(path.c_str());

    if( datRaw_getNext(&info, (void**)&dataBuffer, DR_FORMAT_FLOAT) <= 0) {
        datRaw_close(&info);
        datRaw_freeInfo(&info);
        return 0;
    }

    datRaw_close(&info);
    datRaw_freeInfo(&info);
    //---------------------------------------------------------------------------

    vtkFloatArray *data = vtkFloatArray::New();

    int numOutputComponents = this->NumberOfComponents;
    if( this->NumberOfComponents == 6 )
        numOutputComponents = 9;

    data->SetNumberOfComponents(numOutputComponents);
    int res[3] = {Extent[1]-Extent[0]+1,
        Extent[3]-Extent[2]+1,
        Extent[5]-Extent[4]+1};
    data->SetNumberOfTuples(res[0]*res[1]*res[2]);
    data->SetName("Data");

    float *tuple = new float[numOutputComponents];


    int offsets[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    const int *off = offsets;

    if (this->NumberOfComponents == 6) {
        if( this->TensorFormat == XX_YY_ZZ_XY_XZ_YZ ) {
            off = tensF_0;
        } else {
            off = tensF_1;
        }
    }

    for (int i = 0; i < res[0]*res[1]*res[2]; i++) 
    {
        for (int j = 0; j < numOutputComponents; j++) 
        {
            tuple[j] = dataBuffer[i*NumberOfComponents+off[j]];
        }
        data->SetTuple(i,tuple);
    }

    delete [] tuple;
    delete [] dataBuffer;

    output->SetExtent(Extent);
    output->SetOrigin(Origin);
    output->SetSpacing(Spacing);
    output->GetPointData()->AddArray(data);

    data->Delete();

    return 1;
}

void DatRawTemporalReader::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);

    os << indent << "File Name: "
        << (this->FileName ? this->FileName : "(none)") << "\n";
}
