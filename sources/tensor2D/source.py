gridX=11
gridY=11

numNodes = gridX * gridY
numCells = 2*(gridX-1)*(gridY-1)

midpoint = ( float(gridX-1) / 2.0, float(gridY-1) / 2.0)
scale    = ( 11.0 / gridX,         11.0 / gridY )

pts = vtk.vtkPoints()
tensor = vtk.vtkFloatArray()
tensor.SetNumberOfComponents(4)
tensor.SetName('tensor')

for ix in range(gridX):
  for iy in range(gridY):
    x = scale[0] * (ix-midpoint[0])
    y = scale[1] * (midpoint[1]-iy)
    z = 0.0
    pts.InsertNextPoint(x,y,z)

    t11 = 1.0-2.0*x
    t12 = y
    t21 = y
    t22 = 1.0

    tensor.InsertNextTuple4( t11, t12, t21, t22 )

try:
  grid = self.GetUnstructuredGridOutput()
except:
  grid = vtk.vtkUnstructuredGrid()

grid.SetPoints(pts)
grid.GetPointData().AddArray(tensor)

# create triangle mesh
cellSize = 3
cellType = 5 # vtk_triangle

grid.Allocate(numCells)

s = gridY
Ids = vtk.vtkIdList()
for ix in range( gridX-1):
  for iy in range( gridY-1):
    i = ix * gridX + iy
    Ids.Reset()
    Ids.InsertNextId(i)
    Ids.InsertNextId(i+1)
    Ids.InsertNextId(i+s)
    grid.InsertNextCell( cellType, Ids)
    Ids.Reset()
    Ids.InsertNextId(i+1)
    Ids.InsertNextId(i+s+1)
    Ids.InsertNextId(i+s)
    grid.InsertNextCell( cellType, Ids)

