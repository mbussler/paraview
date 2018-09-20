CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(TetClip)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

#linalg
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../common/)

ADD_PARAVIEW_PLUGIN( TetClip "1.0"
  SERVER_MANAGER_XML TetClip.xml
  SERVER_MANAGER_SOURCES TetClip.cxx  
)
  
INSTALL( TARGETS TetClip CONFIGURATIONS Release 
                                 RUNTIME DESTINATION bin/Release 
                                 LIBRARY DESTINATION bin/Release
                                 )

INSTALL( TARGETS TetClip CONFIGURATIONS Debug 
                                 RUNTIME DESTINATION bin/Debug 
                                 LIBRARY DESTINATION bin/Debug
                                 )
INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS RelWithDebInfo RUNTIME DESTINATION bin/RelWithDebInfo LIBRARY DESTINATION bin/RelWithDebInfo )