cmake_minimum_required(VERSION 2.8)

project(TetRidge)

# #--- CudaRidge Lib----------------------------------

# FIND_PACKAGE(CUDA REQUIRED)
# CUDA_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
# CUDA_ADD_LIBRARY(cudaRidge SHARED cudaRidge.cu OPTIONS -arch=sm_20 --pre-include cudaRidge.h )
# #--------------------------------------------------

#linalg
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../common/)

IF (ParaView_SOURCE_DIR)
  INCLUDE_DIRECTORIES(
    ${VTK_INCLUDE_DIRS}
  )
ELSE (ParaView_SOURCE_DIR)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF (ParaView_SOURCE_DIR)

ADD_PARAVIEW_PLUGIN(TetRidge "0.1"
  SERVER_MANAGER_XML TetRidge.xml
  SERVER_MANAGER_SOURCES vtkTetRidge.cxx
)

#TARGET_LINK_LIBRARIES(Ridge PRIVATE cudaRidge )

INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS Release 
                                 RUNTIME DESTINATION bin/Release 
                                 LIBRARY DESTINATION bin/Release
                                 )

INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS Debug 
                                 RUNTIME DESTINATION bin/Debug 
                                 LIBRARY DESTINATION bin/Debug
                                 )
INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS RelWithDebInfo 
                                 RUNTIME DESTINATION bin/RelWithDebInfo 
                                 LIBRARY DESTINATION bin/RelWithDebInfo )
