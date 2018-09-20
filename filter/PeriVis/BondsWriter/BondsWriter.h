/*=========================================================================

  Program:   ParaView
  Module:    BondsWriter.h

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME BondsWriter - CSV writer for vtkTable
// Writes a vtkTable as a delimited text file (such as CSV). 
#ifndef __BondsWriter_h
#define __BondsWriter_h

#include "vtkWriter.h"

class vtkStdString;
class vtkPolyData;

class BondsWriter : public vtkWriter
{
public:
  static BondsWriter* New();
  vtkTypeMacro(BondsWriter, vtkWriter);

  // Description:
  // Get/Set the filename for the file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);


  // Description:
  // Get/Set the precision to use for printing numeric values.
  // Default is 5.
  vtkSetClampMacro(Precision, int, 0, VTK_INT_MAX);
  vtkGetMacro(Precision, int);

protected:
  BondsWriter();
  ~BondsWriter();

  bool OpenFile();

  virtual void WriteData();
  virtual void WriteBonds(vtkPolyData* input);

  // see algorithm for more info.
  // This writer takes in vtkTable.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  char* FileName;
  ofstream* Stream;
  int Precision;

private:
  BondsWriter(const BondsWriter&); // Not implemented.
  void operator=(const BondsWriter&); // Not implemented.
};


#endif

