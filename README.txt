
How to build my ParaView-Plugins

1) Build Paraview. Version 4.4 or above should work
2) Plugins are in the filter-directory. Create a build folder and run ccmake or use CMakeGUI:

> mkdir build && cd build && ccmake ..

Important settings are:
ParaView_DIR:               set to where you build ParaView 
CMAKE_INSTALL_PREFIX        set to where the plugins should be installed

3) Prerequisites

For the TetGen plugin (filter/PeriVis/TetGen) the tetlib must be build and installed:
> cd filter/PeriVis/TetGen/Tetgen && mkdir build && cd build && cmake .. && make install
or use CMakeGui and build the INSTALL target in Visual Studio.
This should install the tetlib in filter/PeriVis/TetGen/Tetgen/lib

For the StreamTracerEV plugin, a patch must be applied to VTK.
Therefore, the files vtkAbstractInterpolatedVelocityField.cxx and vtkAbstractInterpolatedVelocityField.h
need to be merged with VTK/Filters/FlowPaths/vtkAbstractInterpolatedVelocityField.* or your VTK installation.
Please also check the *.patch-files.

3) Build and install plugins:

> make install

This should install the given directory defined by CMAKE_INSTALL_PREFIX, with bin/Release appended.

3) Load Plugins

To load the plugins in ParaView, set the PV_PLUGIN_PATH environment variable to where the plugins are installed.


BONUS BONUS BONUS BONUS BONUS BONUS BONUS BONUS BONUS BONUS BONUS BONUS BONUS BONUS BONUS BONUS BONUS 

Check the "screenshot_patch.cpp" to save a screenshot along with every state file in ParaView!