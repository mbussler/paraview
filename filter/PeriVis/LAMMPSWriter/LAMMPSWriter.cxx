/*=========================================================================

  Program:   ParaView
  Module:    LAMMPSWriter.cxx

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "LAMMPSWriter.h"

#include "vtkAlgorithm.h"
#include "vtkArrayIteratorIncludes.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkErrorCode.h"
#include "vtkInformation.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkPolyLineToRectilinearGridFilter.h"
#include "vtkTable.h"
#include "vtkSmartPointer.h"

#include <vector>
//#include <vtksys/ios/sstream>

vtkStandardNewMacro(LAMMPSWriter);
//-----------------------------------------------------------------------------
LAMMPSWriter::LAMMPSWriter()
{
  this->FileName = 0;
  this->Precision = 5;

}

//-----------------------------------------------------------------------------
LAMMPSWriter::~LAMMPSWriter()
{
  this->SetFileName(0);
  delete this->Stream;
}

//-----------------------------------------------------------------------------
int LAMMPSWriter::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

//-----------------------------------------------------------------------------
bool LAMMPSWriter::OpenFile()
{
  if ( !this->FileName )
    {
    vtkErrorMacro(<< "No FileName specified! Can't write!");
    this->SetErrorCode(vtkErrorCode::NoFileNameError);
    return false;
    }

  vtkDebugMacro(<<"Opening file for writing...");

  ofstream *fptr = new ofstream(this->FileName, ios::out);

  if (fptr->fail())
    {
    vtkErrorMacro(<< "Unable to open file: "<< this->FileName);
    this->SetErrorCode(vtkErrorCode::CannotOpenFileError);
    delete fptr;
    return false;
    }

  this->Stream = fptr;
  return true;
}

//-----------------------------------------------------------------------------
void LAMMPSWriter::WriteData()
{
  vtkPointSet* ds = vtkPointSet::SafeDownCast(this->GetInput());
  if (ds && ds->GetPointData())
    {
    this->WritePoints(ds);
    }
  else
    {
    vtkErrorMacro(<< "LAMMPSWriter can only write point data.");
    }
}

//-----------------------------------------------------------------------------
void LAMMPSWriter::WritePoints(vtkPointSet* input)
{
  vtkIdType numPoints = input->GetNumberOfPoints();
  if (!this->OpenFile())
    {
    return;
    }

  double b[6];
  input->GetBounds(b);

  // // push the floating point precision/notation type.
  //if (this->UseScientificNotation)
  //{
	 // (*this->Stream) << std::scientific;
  //}
  (*this->Stream) << std::setprecision(this->Precision);
  //(*this->Stream) << std::setprecision(4);

  // Write headers:
  (*this->Stream) << "ITEM: TIMESTEP"             << '\n';
  (*this->Stream) << "0"                          << '\n';
  (*this->Stream) << "ITEM: NUMBER OF ATOMS"      << '\n';
  (*this->Stream) << numPoints                    << '\n';
  (*this->Stream) << "ITEM: BOX BOUNDS ss ss ss"  << '\n';
  (*this->Stream) << b[0] << ' ' << b[1]          << '\n';
  (*this->Stream) << b[2] << ' ' << b[3]          << '\n';
  (*this->Stream) << b[4] << ' ' << b[5]          << '\n';
  (*this->Stream) << "ITEM: ATOMS id bond_type volume dichte x y z"       << '\n';

  double pos[3];
  for( vtkIdType ptId = 0; ptId<numPoints; ptId++)
  {
	  input->GetPoint(ptId, pos);
      (*this->Stream) << (ptId+1) << ' ';
      (*this->Stream) << BondType << ' ';
      (*this->Stream) << Volume << ' ';
	  (*this->Stream) << Dichte << ' ';
	  (*this->Stream) << pos[0]   << ' ';
	  (*this->Stream) << pos[1]   << ' ';
	  (*this->Stream) << pos[2]   << '\n';
  }

  this->Stream->close();
}
