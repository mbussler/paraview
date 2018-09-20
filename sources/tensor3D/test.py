#!/usr/bin/python
from vtk import *

def printCells(grid):
  for ci in range(grid.GetNumberOfCells()):
    c=grid.GetCell(ci)
    print 'cell '+str(ci)+': [',
    for id in range(c.GetNumberOfPoints()):
      print str(c.GetPointId(id)) + ' ',
    print ']'


exec( file( 'source.py').read() )

printCells(outputGrid)


