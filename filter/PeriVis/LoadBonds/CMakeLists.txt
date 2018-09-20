CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(LoadBonds)

# find and include Eigen
SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/../../cmake)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

IF (ParaView_SOURCE_DIR)
  INCLUDE_DIRECTORIES(
    ${VTK_INCLUDE_DIRS}
  )  
ELSE (ParaView_SOURCE_DIR)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF (ParaView_SOURCE_DIR)

ADD_PARAVIEW_PLUGIN( LoadBonds "1.0"
  SERVER_MANAGER_XML LoadBonds.xml
  SERVER_MANAGER_SOURCES LoadBonds.cxx
)

TARGET_LINK_LIBRARIES( ${PROJECT_NAME} PRIVATE pqCore )

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