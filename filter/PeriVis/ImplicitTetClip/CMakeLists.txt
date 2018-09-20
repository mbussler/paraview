CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(ImplicitTetClip)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

#linalg
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../common/)

ADD_PARAVIEW_PLUGIN( ImplicitTetClip "1.0"
  SERVER_MANAGER_XML ImplicitTetClip.xml
  SERVER_MANAGER_SOURCES ImplicitTetClip.cxx  
)
  
INSTALL( TARGETS ImplicitTetClip CONFIGURATIONS Release 
                                 RUNTIME DESTINATION bin/Release 
                                 LIBRARY DESTINATION bin/Release
                                 )

INSTALL( TARGETS ImplicitTetClip CONFIGURATIONS Debug 
                                 RUNTIME DESTINATION bin/Debug 
                                 LIBRARY DESTINATION bin/Debug
                                 )
INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS RelWithDebInfo RUNTIME DESTINATION bin/RelWithDebInfo LIBRARY DESTINATION bin/RelWithDebInfo )