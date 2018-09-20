CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(ClipLines)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

ADD_PARAVIEW_PLUGIN( ClipLines "1.0"
  SERVER_MANAGER_XML ClipLines.xml
  SERVER_MANAGER_SOURCES ClipLines.cxx  
)
  
INSTALL( TARGETS ClipLines CONFIGURATIONS Release 
                                 RUNTIME DESTINATION bin/Release 
                                 LIBRARY DESTINATION bin/Release
                                 )

INSTALL( TARGETS ClipLines CONFIGURATIONS Debug 
                                 RUNTIME DESTINATION bin/Debug 
                                 LIBRARY DESTINATION bin/Debug
                                 )
INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS RelWithDebInfo RUNTIME DESTINATION bin/RelWithDebInfo LIBRARY DESTINATION bin/RelWithDebInfo )