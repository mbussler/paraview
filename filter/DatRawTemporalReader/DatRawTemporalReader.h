#ifndef __DatRawTemporalReader_h
#define __DatRawTemporalReader_h
 
#include "vtkImageAlgorithm.h"
#include <vector>
#include <string>

class DatRawTemporalReader : public vtkImageAlgorithm
{
public:
  vtkTypeMacro(DatRawTemporalReader,vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  static DatRawTemporalReader *New();
 
  // Description:
  // Specify file name of the .dat file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  vtkSetMacro(TensorFormat, int);
  vtkGetMacro(TensorFormat, int);
 
protected:
  DatRawTemporalReader();
  ~DatRawTemporalReader(){}
 
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*);

private:
  DatRawTemporalReader(const DatRawTemporalReader&);  // Not implemented.
  void operator=(const DatRawTemporalReader&);  // Not implemented.

  void ReadDatList();
 
  char* FileName;

  int m_NumberOfTimeSteps;
  std::vector<double> m_TimestepValues;
  std::vector<std::string> m_DatList;

  int FindClosestTimeStep(double requestedTimeValue);

  int Extent[6];
  int NumberOfComponents;
  int TensorFormat;
  double Origin[3];
  double Spacing[3];
  double Time;
  bool bRelativePath;
};
 
#endif//__DatRawTemporalReader_h
