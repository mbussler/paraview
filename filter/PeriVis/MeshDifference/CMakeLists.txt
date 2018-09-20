CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(MeshDifference)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

ADD_PARAVIEW_PLUGIN( ${PROJECT_NAME} "1.0"
  SERVER_MANAGER_XML MeshDifference.xml
  SERVER_MANAGER_SOURCES MeshDifference.cxx  
)
  
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
                                 LIBRARY DESTINATION bin/RelWithDebInfo 
                                 )