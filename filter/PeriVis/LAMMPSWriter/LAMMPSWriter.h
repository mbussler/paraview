/*=========================================================================

  Program:   ParaView
  Module:    LAMMPSWriter.h

  Copyright (c) Kitware, Inc.
  All rights reserved.
  See Copyright.txt or http://www.paraview.org/HTML/Copyright.html for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME LAMMPSWriter - CSV writer for vtkTable
// Writes a vtkTable as a delimited text file (such as CSV). 
#ifndef __LAMMPSWriter_h
#define __LAMMPSWriter_h

#include "vtkWriter.h"

class vtkStdString;
class vtkPointSet;

class LAMMPSWriter : public vtkWriter
{
public:
  static LAMMPSWriter* New();
  vtkTypeMacro(LAMMPSWriter, vtkWriter);

  // Description:
  // Get/Set the filename for the file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);


  // Description:
  // Get/Set the precision to use for printing numeric values.
  // Default is 5.
  vtkSetClampMacro(Precision, int, 0, VTK_INT_MAX);
  vtkGetMacro(Precision, int);

  vtkSetMacro( BondType, int);
  vtkGetMacro( BondType, int);
  vtkSetMacro( Volume, double) ;
  vtkGetMacro( Volume, double) ;
  vtkSetMacro( Dichte, double) ;
  vtkGetMacro( Dichte, double) ;

protected:
  LAMMPSWriter();
  ~LAMMPSWriter();

  bool OpenFile();

  virtual void WriteData();
  virtual void WritePoints(vtkPointSet* input);

  // see algorithm for more info.
  // This writer takes in vtkTable.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  char* FileName;
  ofstream* Stream;
  int Precision;
  int BondType;
  double Volume;
  double Dichte;

private:
  LAMMPSWriter(const LAMMPSWriter&); // Not implemented.
  void operator=(const LAMMPSWriter&); // Not implemented.
};


#endif

