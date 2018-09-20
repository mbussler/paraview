/*=========================================================================

  Program:   Visualization Toolkit
  Module:    LoadBonds.h

  Copyright (c) Michael Bußler

=========================================================================*/

#ifndef __ConvertToPolydata_h
#define __ConvertToPolydata_h

#include "vtkPolyDataAlgorithm.h"

#include <vtkSmartPointer.h>
#include <vector>
#include <string>
class vtkPolyData;

#if defined(_MSC_VER) // visual studio is.. different!
    typedef          __int64 int64_t;
    typedef unsigned __int64 uint64_t;
#else
    #include <inttypes.h>
#endif

class Timestep {
public:
    double time;
    int fileId;
};

class LoadBonds : public vtkPolyDataAlgorithm
{
  
  public:
   
    vtkTypeMacro(LoadBonds, vtkPolyDataAlgorithm);

    void PrintSelf(ostream &os, vtkIndent indent);
    
    // Description
    // Calculate Eigenvectors for tensor data
    static LoadBonds *New();
   
    vtkGetMacro(CalculateBrokenBonds, bool);
    vtkSetMacro(CalculateBrokenBonds, bool);
    vtkGetMacro(CalculateBondLength, bool);
    vtkSetMacro(CalculateBondLength, bool);
    vtkGetMacro(CalculateRefBondLength, bool);
    vtkSetMacro(CalculateRefBondLength, bool);
    vtkGetMacro(MaxBondLength, double);
    vtkSetMacro(MaxBondLength, double);
    vtkGetMacro(FilterBondsByLength, double);
    vtkSetMacro(FilterBondsByLength, double);
   
    // Description:
    // Adds names of files to be read. The files are read in the order
    // they are added.
    virtual void AddFileName(const char* fname);

    // Description:
    // Remove all file names.
    virtual void RemoveAllFileNames();

protected:
    LoadBonds();
    ~LoadBonds();

//     int RequestUpdateExtent(vtkInformation *,  vtkInformationVector **, vtkInformationVector *);
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
    int RequestInformation(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector);
    int FillInputPortInformation(int port, vtkInformation *info);
    virtual int RequestUpdateExtent(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector);

  private:
    LoadBonds(const LoadBonds&);  // Not implemented.
    void operator=(const LoadBonds&);  // Not implemented.

    void ParseBondsFile();
    void ProcessBondData( int timestepToLoad, vtkPolyData* pd);
    int FindClosestTimeStep(double requestedTimeValue);
    void tokenize(const std::string line, std::vector<std::string> &tokens);
    vtkSmartPointer<vtkIdList> getConnectedPoints( vtkPolyData * pd, vtkIdType ptId, vtkSmartPointer<vtkIdList> bondIds = NULL);


    bool CalculateBrokenBonds;
    bool CalculateBondLength;
    bool CalculateRefBondLength;
    bool FilterBondsByLength;
    double MaxBondLength;

    int NumberOfTimeSteps;
    std::vector<double> TimestepValues;
    std::vector<Timestep> m_timesteps;

    vtkSmartPointer<vtkPolyData> m_initialMesh;
    std::vector<std::string> m_FileNames;

    int m_requested_time;
};

#endif
