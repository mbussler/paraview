#ifndef __LAMMPSTemporalReader_h
#define __LAMMPSTemporalReader_h

#include <vtkPolyDataAlgorithm.h>
#include <vtkSmartPointer.h>
#include <vector>
#include <string>

class vtkPolyData;
class vtkImplicitFunction;

#if defined(_MSC_VER) // visual studio is.. different!
 typedef          __int64 int64_t;
 typedef unsigned __int64 uint64_t;
#else
 #include <inttypes.h>
#endif

class Atom {
public:
	int id;
	int type;
	float pos[3];
	float nDataFields;
	float* data;
	Atom() : data(0) {};
	~Atom(){ if( data ) delete[] data; }
};

class Timestep {
public:
	int time;
	int nAtoms, nBonds;
	float box[6];
	std::vector<std::string> dataFieldNames;
	int nDataFields;
	int posField[3];
	int bondField[2];
	uint64_t atomStartLine;
	uint64_t bondsStartLine;
};

class LAMMPSTemporalReaderMF : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(LAMMPSTemporalReaderMF,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  static LAMMPSTemporalReaderMF *New();
 
  // Description:
  // Specify i-j-k dimensions on which to sample input points.
  vtkGetVectorMacro(SampleDimensions,int,3);

  // Description:
  // Set the i-j-k dimensions on which to sample the distance function.
  void SetSampleDimensions(int i, int j, int k);

  // Description:
  // Set the i-j-k dimensions on which to sample the distance function.
  void SetSampleDimensions(int dim[3]);

  vtkGetMacro(CalculateBrokenBonds, bool);
  vtkSetMacro(CalculateBrokenBonds, bool);
  vtkGetMacro(CalculateBondLength, bool);
  vtkSetMacro(CalculateBondLength, bool);
  vtkGetMacro(CalculateRefBondLength, bool);
  vtkSetMacro(CalculateRefBondLength, bool);
  vtkGetMacro(LoadBonds, bool);
  vtkSetMacro(LoadBonds, bool);
  vtkGetMacro(ShowLines, bool);
  vtkSetMacro(ShowLines, bool);

  // Description:
  // Adds names of files to be read. The files are read in the order
  // they are added.
  virtual void AddFileName(const char* fname);

  // Description:
  // Remove all file names.
  virtual void RemoveAllFileNames();


protected:
  LAMMPSTemporalReaderMF();
  ~LAMMPSTemporalReaderMF();
 
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int RequestInformation(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
  int FillOutputPortInformation( int port, vtkInformation* info );

private:
  LAMMPSTemporalReaderMF(const LAMMPSTemporalReaderMF&);  // Not implemented.
  void operator=(const LAMMPSTemporalReaderMF&);  // Not implemented.

  void ParseLAMMPSFile();
  void ParseBondsFile();
 
  int NumberOfTimeSteps;
  std::vector<Timestep> m_timesteps;
  std::vector<double> TimestepValues;

  int FindClosestTimeStep(double requestedTimeValue);
  void ReadAtoms( int id, vtkPolyData* pd );
  void tokenize(const std::string line, std::vector<std::string> &tokens);
  void ReadBonds(int timestepToLoad, vtkPolyData * pd);

  vtkSmartPointer<vtkIdList> getConnectedPoints(vtkPolyData * pd, vtkIdType ptId, vtkSmartPointer<vtkIdList> bondIds = NULL);

  int NumberOfComponents;
  bool bRelativePath;
  bool CalculateBrokenBonds;
  bool CalculateBondLength;
  bool CalculateRefBondLength;
  bool LoadBonds;
  bool ShowLines;

  double ModelBounds[6];
  int SampleDimensions[3];

  vtkSmartPointer<vtkPolyData> m_initialMesh;
  std::vector<std::string> m_FileNames;
  std::vector<std::string> m_BondFileNames;
};
 
#endif//__LAMMPSTemporalReader_h
