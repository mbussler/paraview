CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(ClipPolyData)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

ADD_PARAVIEW_PLUGIN( ClipPolyData "1.0"
  SERVER_MANAGER_XML ClipPolyData.xml
  SERVER_MANAGER_SOURCES ClipPolyData.cxx  
)
  
INSTALL( TARGETS ClipPolyData CONFIGURATIONS Release 
                                 RUNTIME DESTINATION bin/Release 
                                 LIBRARY DESTINATION bin/Release
                                 )

INSTALL( TARGETS ClipPolyData CONFIGURATIONS Debug 
                                 RUNTIME DESTINATION bin/Debug 
                                 LIBRARY DESTINATION bin/Debug
                                 )
INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS RelWithDebInfo RUNTIME DESTINATION bin/RelWithDebInfo LIBRARY DESTINATION bin/RelWithDebInfo )