# Usage:
#
#   Pipeline:
#     TODO: Ridge needs pipeline description
#
#   Parameters:
#     TODO: Ridge needs parameter description
#
#   Contribution :
#     TODO: Ridge needs contributor mention
#
#   Author:
#     TODO: Ridge needs author mention

PROJECT(Ridge)

#---CUDA SPECIFIC----------------------------------

FIND_PACKAGE(CUDA REQUIRED)
CUDA_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
CUDA_ADD_LIBRARY(cudaRidge cudaRidge.cu OPTIONS -arch=sm_20 --pre-include cudaRidge.h
)

#--------------------------------------------------

cmake_minimum_required(VERSION 2.8)

IF (ParaView_SOURCE_DIR)
  INCLUDE_DIRECTORIES(
    ${VTK_INCLUDE_DIRS}
  )
ELSE (ParaView_SOURCE_DIR)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF (ParaView_SOURCE_DIR)

LINK_DIRECTORIES(${LOCAL_CUDA_BIN})

ADD_PARAVIEW_PLUGIN(Ridge "0.1"
  SERVER_MANAGER_XML Ridge.xml
  SERVER_MANAGER_SOURCES vtkRidge.cxx
  SOURCES  isosurface.cpp
)

TARGET_LINK_LIBRARIES(Ridge PRIVATE cudaRidge)


INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS Release 
                                 RUNTIME DESTINATION bin/Release 
                                 LIBRARY DESTINATION bin/Release
                                 )

INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS Debug 
                                 RUNTIME DESTINATION bin/Debug 
                                 LIBRARY DESTINATION bin/Debug
                                 )