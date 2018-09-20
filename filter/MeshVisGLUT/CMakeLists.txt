
cmake_minimum_required(VERSION 2.6)

project(MeshVis)

option( USE_TETGEN "Include support tetgen generated tetraedral grids" FALSE )

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR/cmake}")

# Libraries
FIND_PACKAGE( VTK REQUIRED COMPONENTS vtkIOCore vtkIOLegacy vtkIOXML vtksys vtkFiltersCore vtkCommonCore vtkCommonDataModel vtkCommonExecutionModel )
FIND_PACKAGE( Boost REQUIRED)
FIND_PACKAGE( CUDA REQUIRED)
FIND_PACKAGE( GLEW REQUIRED)
FIND_PACKAGE( GLUT REQUIRED)
FIND_PACKAGE( OpenGL REQUIRED)

IF( USE_TETGEN )
	add_definitions( -DTETLIBRARY )
	add_definitions( -DSINGLE )
	include_directories( "${CMAKE_CURRENT_SOURCE_DIR}/3rdParty/Tetgen" )
	link_directories( "${CMAKE_CURRENT_SOURCE_DIR}/3rdParty/Tetgen/lib")
	add_subdirectory( 3rdParty/Tetgen )
ENDIF()

include_directories( "${CMAKE_CURRENT_SOURCE_DIR}/3rdParty/glh" )

IF( CUDA_FOUND AND Boost_FOUND AND VTK_FOUND )

	include_directories( ${VTK_INCLUDE_DIRS} )
	include_directories( ${Boost_INCLUDE_DIRS} )
	include_directories( ${GLEW_INCLUDE_DIRS} )
	include_directories( ${GLUT_INCLUDE_DIRS} )
  include_directories( ${OPENGL_INCLUDE_DIR} )
  include_directories( ${CUDA_SDK_ROOT_DIR}/common/inc)
	add_subdirectory( src )
ENDIF()

list( APPEND CMAKE_CXX_FLAGS "-std=c++0x" ${CMAKE_CXX_FLAGS} )