#!/usr/bin/env python

# This example demonstrates the use of vtkVectorText and vtkFollower.
# vtkVectorText is used to create 3D annotation.  vtkFollower is used to
# position the 3D text and to ensure that the text always faces the
# renderer's active camera (i.e., the text is always readable).

import vtk
from vtk import *

# load data

# The source file
file_name = "data/test.vtu"
 
# Read the source file.
reader = vtkXMLUnstructuredGridReader()
reader.SetFileName(file_name)
reader.Update() # Needed because of GetScalarRange
output = reader.GetOutput()
scalar_range = output.GetScalarRange()

# Create the mapper that corresponds the objects of the vtk file
# into graphics elements
mapper = vtkDataSetMapper()
mapper.SetInputConnection(reader.GetOutputPort())
mapper.SetScalarRange(scalar_range)

# Create the Actor
actor = vtkActor()
actor.SetMapper(mapper)
actor.SetScale(100.0,100.0,100.0)

# Create the axes and the associated mapper and actor.
axes = vtk.vtkAxes()
axes.SetOrigin(0, 0, 0)
axesMapper = vtk.vtkPolyDataMapper()
axesMapper.SetInputConnection(axes.GetOutputPort())
axesActor = vtk.vtkActor()
axesActor.SetMapper(axesMapper)

# Create the 3D text and the associated mapper and follower (a type of
# actor).  Position the text so it is displayed over the origin of the
# axes.
atext = vtk.vtkVectorText()
atext.SetText("Origin")
textMapper = vtk.vtkPolyDataMapper()
textMapper.SetInputConnection(atext.GetOutputPort())
textActor = vtk.vtkFollower()
textActor.SetMapper(textMapper)
textActor.SetScale(0.2, 0.2, 0.2)
textActor.AddPosition(0, -0.1, 0)

# Create the Renderer, RenderWindow, and RenderWindowInteractor.
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

# Add the actors to the renderer.
ren.AddActor(axesActor)
ren.AddActor(textActor)
ren.AddActor(actor)

# Zoom in closer.
ren.ResetCamera()
ren.GetActiveCamera().Zoom(1.6)

# Reset the clipping range of the camera; set the camera of the
# follower; render.
ren.ResetCameraClippingRange()
textActor.SetCamera(ren.GetActiveCamera())

iren.Initialize()
renWin.Render()
iren.Start()
