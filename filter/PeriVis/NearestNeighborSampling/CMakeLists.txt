CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(NearestNeighborSampling)

SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/../../../cmake)

# find and include GSL
FIND_PACKAGE(GSL REQUIRED)
include_directories( ${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})

#linalg
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../common/)

# we are making a paraview plugin, aren't we?
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

set( LIBS 
  ${GSL_LIBRARIES} 
  ${GSLCBLAS_LIBRARIES}
)

ADD_PARAVIEW_PLUGIN( NearestNeighborSampling "1.0"
  SERVER_MANAGER_XML NearestNeighborSampling.xml
  SERVER_MANAGER_SOURCES NearestNeighborSampling.cxx
)
target_link_libraries (${PROJECT_NAME} PRIVATE ${LIBS})

# make distinguishable install directories
INSTALL(
  TARGETS ${PROJECT_NAME} CONFIGURATIONS Release 
  RUNTIME DESTINATION bin/Release 
  LIBRARY DESTINATION bin/Release
)

INSTALL( 
  TARGETS ${PROJECT_NAME} CONFIGURATIONS Debug 
  RUNTIME DESTINATION bin/Debug 
  LIBRARY DESTINATION bin/Debug
)
INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS RelWithDebInfo RUNTIME DESTINATION bin/RelWithDebInfo LIBRARY DESTINATION bin/RelWithDebInfo )