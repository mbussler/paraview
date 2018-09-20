

#include "LAMMPSReader.h"

#include "vtkPolyData.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkSmartPointer.h"
#include "vtkClipPolyData.h"
#include "vtkImplicitFunction.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <functional>
#include <direct.h>

inline std::string getPath( const std::string& str )
{
	std::string s(str);
	std::replace( s.begin(), s.end(), '\\', '/'); // replace \ with /
	return s.substr(0, s.find_last_of("/"));
}

LAMMPSReader::LAMMPSReader( const std::string& atomFileName, const std::string& bondsFileName )
{
	this->FileName = atomFileName;
	this->BondsFileName = bondsFileName;
	this->NumberOfTimeSteps=0;
	this->bRelativePath = false;
	this->func = NULL;

	if( !this->FileName.empty() ){

		ReadLAMMPSFile();

		if( !this->BondsFileName.empty()) {
			ReadBondsFile();
		}
	}
}

void LAMMPSReader::ReadLAMMPSFile()
{
	using namespace std;

	m_timesteps.clear();
	TimestepValues.clear();

	ifstream inFile;
	inFile.open(FileName, ifstream::in);

	Timestep ts;
	ts.atomStartLine = -1;
	ts.time = -1;
	ts.nBonds = -1;
	ts.BondsStartLine = -1;

	string line;
	int lineId = -1;

	while (inFile.good()) {

		getline(inFile, line);
		lineId++;

		vector<string> tokens;
		tokenize(line, tokens);

		if( tokens.size() == 0 ) {
			continue;
		}

		if (tokens[0] == "ITEM:" )
		{
			if( tokens.size() > 1 && tokens[1] == "TIMESTEP") 
			{
				// save current timestep
				if( ts.time >= 0 )
					m_timesteps.push_back(ts);

				// read time
				getline(inFile, line);lineId++;
				
				// init new timestep
				ts = Timestep();
				ts.time = stoi(line);
				ts.nBonds = -1;
				ts.BondsStartLine = -1;
				ts.atomStartLine = -1;
				ts.dataFieldNames.clear();
				ts.nAtoms = -1;
				ts.nDataFields = 0;
			}
			else if ( tokens.size() > 3 && tokens[1] == "NUMBER" && tokens[2] == "OF" && tokens[3] == "ATOMS") 
			{
				// read number of atoms
				getline(inFile, line);lineId++;
				ts.nAtoms = stoi(line);
			}
			else if ( tokens.size() > 3 && tokens[1] == "BOX" && tokens[2] == "BOUNDS")
			{
				// read bounding box
				string::size_type sz;
				for( int i=0; i<3; i++){
					getline(inFile, line);lineId++;
					ts.box[2*i+0] = stof(line, &sz);
					ts.box[2*i+1] = stof(line.substr(sz));
				}
			}
			else if ( tokens.size() > 1 && tokens[1] == "ATOMS" )
			{
				// read atom data
				for( int i=2; i<tokens.size();i++)
				{
					if( tokens[i] == "x" )
						ts.posField[0]=i-2;
					else if( tokens[i] == "y" )
						ts.posField[1]=i-2;
					else if( tokens[i] == "z" )
						ts.posField[2]=i-2;
					else 
						ts.dataFieldNames.push_back(tokens[i]);
				}
				ts.nDataFields = ts.dataFieldNames.size();

				// save line offset
				ts.atomStartLine = lineId+1;

				// skip atom data
				for (int currLineNumber = 0; currLineNumber < ts.nAtoms; ++currLineNumber){
					inFile.ignore(numeric_limits<streamsize>::max(), inFile.widen('\n'));
				}
				lineId+=ts.nAtoms;
			}
		}
	}
	// save last timestep
	if( ts.time >= 0 )
		m_timesteps.push_back(ts);
	
	NumberOfTimeSteps = m_timesteps.size();

	for( int i=0; i<NumberOfTimeSteps; i++)
		TimestepValues.push_back(m_timesteps[i].time);

	inFile.close();
}

void LAMMPSReader::ReadBondsFile()
{
	if( this->BondsFileName.empty() )
		return;

	using namespace std;

	ifstream inFile;
	inFile.open(BondsFileName, ifstream::in);

	string line;
	int lineId = -1;

	Timestep* ts = NULL;

	while (inFile.good()) {

		getline(inFile, line);lineId++;

		vector<string> tokens;
		tokenize(line, tokens);

		if( tokens.size() == 0 ) {
			continue;
		}

		if (tokens[0] == "ITEM:" )
		{
			if( tokens.size() > 1 && tokens[1] == "TIMESTEP")
			{
				// read time
				getline(inFile, line);lineId++;
				int timeFound = stoi(line);
				int timeStepId = FindClosestTimeStep(timeFound);
				if( timeFound == m_timesteps[timeStepId].time ){
					ts = &m_timesteps[timeStepId];
				} else {
					ts = NULL;
				}
			}
			else if ( tokens.size() > 3 && tokens[1] == "NUMBER" && tokens[2] == "OF" && tokens[3] == "ENTRIES") 
			{
				// read number of atoms
				getline(inFile, line);lineId++;
				int numBonds = stoi(line);
				if( ts ) {
					ts->nBonds = numBonds;
				} else {
					// skip entries
					for (int currLineNumber = 0; currLineNumber < numBonds; ++currLineNumber){
						inFile.ignore(numeric_limits<streamsize>::max(), inFile.widen('\n'));
					}
					lineId += numBonds;
				}
			}
			else if ( tokens.size() > 1 && tokens[1] == "ENTRIES" )
			{
				if( ts ) {
					ts->bondField[0]=1;
					ts->bondField[1]=2;
					ts->BondsStartLine = lineId+1;

					// skip entries
					for (int currLineNumber = 0; currLineNumber < ts->nBonds; ++currLineNumber){
						inFile.ignore(numeric_limits<streamsize>::max(), inFile.widen('\n'));
					}
					lineId += ts->nBonds;
				}
			}
		}
	}
	inFile.close();
}


void LAMMPSReader::ReadTimeStep( int id, Timestep& ts )
{
	throw std::exception("Fuction not implemented!");

//	using namespace std;
//
//	if( id < 0 || id > NumberOfTimeSteps )
//		return;
//
//	ts = m_timesteps[id]; // copy meta data
//
//	if( ts.atomStartLine < 0 )
//		return;
//
//	ifstream inFile;
//	inFile.open(this->FileName, ifstream::in);
//
//	int startLine = ts.atomStartLine;
//
//	// skip to start line
//	for (int currLineNumber = 0; currLineNumber < startLine; ++currLineNumber){
//		inFile.ignore(numeric_limits<streamsize>::max(), inFile.widen('\n'));
//	}
//
//	string line;
//	int numLines = ts.nAtoms;
//	while (inFile.good() && numLines ) {
//
//		getline(inFile, line);
//
//		vector<string> tokens;
//		tokenize(line, tokens);
//
//		// Avoid reading empty lines
//		if( tokens.size() == 0)
//			continue;
//
//		if (tokens[0] == "ITEM:" )
//		{
//			if( tokens.size() > 1 && tokens[1] == "TIMESTEP") 
//			{
//				// read time
//				getline(inFile, line);
//				ts.time = stoi(line);
//			}
//			else if ( tokens.size() > 3 && tokens[1] == "NUMBER" && tokens[2] == "OF" && tokens[3] == "ATOMS") 
//			{
//				// read number of atoms
//				getline(inFile, line);
//				ts.nAtoms = stoi(line);
//				ts.atoms.reserve(ts.nAtoms);
//			}
//			else if ( tokens.size() > 3 && tokens[1] == "BOX" && tokens[2] == "BOUNDS")
//			{
//				// read bounding box
//				string::size_type sz;
//				for( int i=0; i<3; i++){
//					getline(inFile, line);
//					ts.box[2*i+0] = stof(line, &sz);
//					ts.box[2*i+1] = stof(line.substr(sz));
//				}
//			}
//			else if ( tokens.size() > 6 && tokens[1] == "ATOMS" )
//			{
//				// read atom data
//				for( int i=7; i<tokens.size();i++)
//					ts.dataFieldNames.push_back(tokens[i]);
//				ts.nDataFields = ts.dataFieldNames.size();
//			}
//		} else if( tokens.size() > 6+ts.nDataFields) {
//
//			//// read atom data
//			//Atom a;
//			//a.id     = stoi(tokens[2]);
//			//a.type   = stoi(tokens[3]);
//			//a.pos[0] = stof(tokens[4]);
//			//a.pos[1] = stof(tokens[5]);
//			//a.pos[2] = stof(tokens[6]);
//			//a.data = new float[ts.nDataFields];
//			//for ( int i=0; i<ts.nDataFields; i++ ){
//			//	a.data[i] = stof(tokens[i+7]);
//			//}
//			//ts.atoms.push_back(a);
//		}
//	}
//
//	inFile.close();
}

void LAMMPSReader::ReadTimeStep( int id, vtkPolyData* pd )
{
	if( id >= 0 && id < NumberOfTimeSteps ){
		Timestep ts = m_timesteps[id];
		ReadTimeStep(ts, pd);

		if( !BondsFileName.empty()){
			ReadTimeStepBonds(ts, pd);
		}
	}
}

void LAMMPSReader::ReadTimeStep( const Timestep& ts, vtkPolyData* pd )
{
	using namespace std;

	if( ts.atomStartLine < 0)
		return;

	ifstream inFile;
	inFile.open(this->FileName, ifstream::in);

	// atom positions
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	// init atom data arrays
	std::vector<vtkSmartPointer<vtkFloatArray> > arrs;
	for( int dataId = 0; dataId< ts.nDataFields; dataId++)
	{
		vtkSmartPointer<vtkFloatArray> data = vtkSmartPointer<vtkFloatArray>::New();
		data->SetName(ts.dataFieldNames[dataId].c_str());
		data->SetNumberOfComponents(1);
		//data->SetNumberOfTuples(points->GetNumberOfPoints());
		arrs.push_back(data);
	}

	// skip to start line of data
	for (int currLineNumber = 0; currLineNumber < ts.atomStartLine; ++currLineNumber){
		inFile.ignore(numeric_limits<streamsize>::max(), inFile.widen('\n'));
	}

	string line;
	int numLines = ts.nAtoms;
	float values[20];
	while (inFile.good() && numLines ) {

		getline(inFile, line);

		// read atom data
		string::size_type sz=0;
		int dataId=0, dataOff=0, sOff=0;
		float x,y,z;
		while( sOff < line.length() )
		{
			try{
				string s = line.substr(sOff);
				if( s == " " )
					break;
				
				float data = stof( s, &sz);
				sOff += sz;

				if( dataOff == ts.posField[0])
					x = data;
				else if ( dataOff == ts.posField[1])
					y = data;
				else if ( dataOff == ts.posField[2])
					z = data;
				else
					//arrs[dataId++]->InsertNextValue(data);
					values[dataId++] = data;
				dataOff++;
			} catch (exception* e) {
				break;
			}
		}

		if( !this->func || this->func->FunctionValue(x,y,z) < 0.0 ) // check if inside implicit function
		{
			points->InsertNextPoint( x,y,z );
			for( int i=0; i<ts.nDataFields; i++)
				arrs[i]->InsertNextValue(values[i]);
			pointIds.push_back((int)values[0]); // store id of point
		}
		numLines--;
	}
	inFile.close();

	pd->SetPoints(points);
	std::vector<vtkSmartPointer<vtkFloatArray> >::iterator a = arrs.begin();
	for( ; a != arrs.end() ; a++ ) {
		pd->GetPointData()->AddArray(*a);
	}
}

void LAMMPSReader::ReadTimeStepBonds( const Timestep& ts, vtkPolyData * pd)
{
	if( BondsFileName.empty() || ts.nBonds < 0 || ts.BondsStartLine < 0)
		return;

	using namespace std;

	ifstream inFile;
	inFile.open(this->BondsFileName, ifstream::in);

	// bonds
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

	// init cell data arrays
	//std::vector<vtkFloatArray*> arrs;
	//for( int dataId = 0; dataId< ts.nDataFields; dataId++){
	//	vtkFloatArray* data = vtkFloatArray::New();
	//	data->SetNumberOfComponents(1);
	//	data->SetNumberOfTuples(points->GetNumberOfPoints());
	//	data->SetName(ts.dataFieldNames[dataId].c_str());
	//	arrs.push_back(data);
	//}

	// skip to start line of data
	for (int currLineNumber = 0; currLineNumber < ts.BondsStartLine; ++currLineNumber){
		inFile.ignore(numeric_limits<streamsize>::max(), inFile.widen('\n'));
	}

	string line;
	int numLines = ts.nBonds;
	while (inFile.good() && numLines>0 ) {

		getline(inFile, line);

		// read bond data
		string::size_type sz=0;
		int dataId=0, dataOff=0, sOff=0;
		int a1=0, a2=0;
		while( sOff < line.length() )
		{
			try{
				string s = line.substr(sOff);
				if( s == " " )
					break;

				int data = stoi( s, &sz);
				sOff += sz;

				if( dataOff == ts.bondField[0])
					a1 = data;
				else if ( dataOff == ts.bondField[1])
					a2 = data;
				//else
				//	arrs[dataId++]->InsertNextValue(data);
				dataOff++;
			} catch (exception* e) {
				break;
			}
		}
		if( a1>0 && a2>0 ){

			if( this->func ) 
			{				
				std::vector<int>::iterator it1, it2;
				it1 = find (pointIds.begin(), pointIds.end(), a1);
				it2 = find (pointIds.begin(), pointIds.end(), a2);
				if( it1 != pointIds.end() && it2 != pointIds.end()) 
				{
					lines->InsertNextCell(2);
					lines->InsertCellPoint(a1-1);
					lines->InsertCellPoint(a2-1);
				}
			} else {
				lines->InsertNextCell(2);
				lines->InsertCellPoint(a1-1);
				lines->InsertCellPoint(a2-1);
			}
		}
		numLines--;

		//this->UpdateProgress ((ts.nAtoms+(ts.nBonds-numLines))/(float)(ts.nAtoms+ts.nBonds));

	}

	inFile.close();
	
	pd->SetLines(lines);

}

int LAMMPSReader::FindClosestTimeStep(double requestedTimeValue)
{
	int ts = 0;
	double mindist = std::abs(TimestepValues[0] - requestedTimeValue);

	for (int i = 0; i < NumberOfTimeSteps; i++) {

		double tsv = TimestepValues[i];
		double dist = std::abs(tsv - requestedTimeValue);
		if (dist < mindist) {
			mindist = dist;
			ts = i;
		}
	}
	return ts;
}

void LAMMPSReader::tokenize(const std::string line, std::vector<std::string> &tokens)
{
	std::stringstream iss(line);
	std::string token;
	while( iss >> token)
		tokens.push_back(token);
}

