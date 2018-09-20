#ifndef __LAMMPSReader_h
#define __LAMMPSReader_h

#include <vector>
#include <string>

class vtkPolyData;
class vtkImplicitFunction;

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
	int atomStartLine;
	int BondsStartLine;
};

class LAMMPSReader
{
public:
  LAMMPSReader( const std::string& atomFileName, const std::string& bondsFileName );
  ~LAMMPSReader(){};
 
  void ReadTimeStep( int id, Timestep& ts );
  void ReadTimeStep( int id, vtkPolyData* pd );

  int FindClosestTimeStep(double requestedTimeValue);

  void setClipFunction( vtkImplicitFunction* func );

  int getNumberOfTimesteps() { return NumberOfTimeSteps; };
  std::vector<double> getTimestepValues() {return TimestepValues; };

protected:
	LAMMPSReader() {};

private:
  LAMMPSReader(const LAMMPSReader&);  // Not implemented.
  void operator=(const LAMMPSReader&);  // Not implemented.

  void ReadLAMMPSFile();
  void ReadBondsFile();

  void tokenize(const std::string line, std::vector<std::string> &tokens);

  void ReadTimeStep(      const Timestep& ts, vtkPolyData* pd);
  void ReadTimeStepBonds( const Timestep& ts, vtkPolyData* pd);

  std::string FileName;
  std::string BondsFileName;

  int NumberOfTimeSteps;
  std::vector<Timestep> m_timesteps;
  std::vector<double> TimestepValues;
  vtkImplicitFunction* func;
  //std::tr1::unordered_map pointIds;
  std::vector<int> pointIds;


  int NumberOfComponents;
  bool bRelativePath;
};
 
#endif//__LAMMPSReader_h
