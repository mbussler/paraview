from vtk.vtkRenderingFreeType import vtkVectorText
from vtk.vtkFiltersCore import vtkAppendPolyData
from vtk.vtkFiltersGeneral import vtkTransformPolyDataFilter
from vtk.vtkCommonTransforms import vtkTransform

text1 = vtkVectorText()
transform = vtkTransform()
transformF = vtkTransformPolyDataFilter()
transformF.SetInputConnection(text1.GetOutputPort())
transformF.SetTransform(transform)
appender = vtkAppendPolyData()

for y in [0, 10, 20]:
    text1.SetText("Text #%d" % y)

    transform.Translate(0, y, 0)
    transformF.Update()
    clone = transformF.GetOutput().NewInstance()
    clone.UnRegister(None)
    clone.ShallowCopy(transformF.GetOutput())

    appender.AddInputDataObject(clone)

appender.Update()

output = self.GetOutput()
output.ShallowCopy(appender.GetOutput())
