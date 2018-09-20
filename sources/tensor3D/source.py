

grid = [11, 11, 11]
size = [22, 22, 22]

numNodes = grid[0] * grid[1] * grid[2]
numCells = (grid[0]-1)*(grid[1]-1)*(grid[2]-1)

midpoint = [(x-1)/2.0 for x in grid]
scale    = ( size[0]/grid[0], size[1]/grid[1], size[2]/grid[2] )

pts = vtk.vtkPoints()
tensor = vtk.vtkFloatArray()
tensor.SetNumberOfComponents(9)
tensor.SetName('tensor')

for iz in range(grid[2]):
  for iy in range(grid[1]):
    for ix in range(grid[0]):
      x = scale[0]*(ix-midpoint[0])
      y = scale[1]*(iy-midpoint[1])
      z = scale[2]*(iz-midpoint[2])
  
      pts.InsertNextPoint(x,y,z)
  
      # diagonals
      t11 = 1.0-2.0*x
      t22 = 1.0
      t33 = 1.0
      # upper left = lower right
      t12 = t21 = y
      t13 = t31 = z
      t23 = t32 = y
  
      tensor.InsertNextTuple9( t11, t12, t13,
                               t21, t22, t23,
                               t31, t32, t33 )

try:
  outputGrid = self.GetUnstructuredGridOutput()
except:
  outputGrid = vtk.vtkUnstructuredGrid()

outputGrid.SetPoints(pts)
outputGrid.GetPointData().AddArray(tensor)

# create triangle mesh
cellSize = 8
cellType = 11 # vtk_voxel

outputGrid.Allocate(numCells)

sz = grid[0]*grid[1]
sy = grid[0]
Ids = vtk.vtkIdList()
for iz in range(grid[2]-1):
  for iy in range(grid[1]-1):
    for ix in range(grid[0]-1):
      i = iz*sz + iy*sy + ix
      Ids.Reset()
      Ids.InsertNextId(i+sy)
      Ids.InsertNextId(i+sy+1)
      Ids.InsertNextId(i)
      Ids.InsertNextId(i+1)
      Ids.InsertNextId(i+sy+sz)
      Ids.InsertNextId(i+sy+sz+1)
      Ids.InsertNextId(i+sz)
      Ids.InsertNextId(i+sz+1)
      outputGrid.InsertNextCell( cellType, Ids)

