CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(PolyPointSet)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

ADD_PARAVIEW_PLUGIN( PolyPointSet "1.0"
  SERVER_MANAGER_XML PolyPointSet.xml
  SERVER_MANAGER_SOURCES PolyPointSet.cxx  
)
  
INSTALL( TARGETS PolyPointSet CONFIGURATIONS Release 
                                 RUNTIME DESTINATION bin/Release 
                                 LIBRARY DESTINATION bin/Release
                                 )

INSTALL( TARGETS PolyPointSet CONFIGURATIONS Debug 
                                 RUNTIME DESTINATION bin/Debug 
                                 LIBRARY DESTINATION bin/Debug
                                 )
INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS RelWithDebInfo RUNTIME DESTINATION bin/RelWithDebInfo LIBRARY DESTINATION bin/RelWithDebInfo )