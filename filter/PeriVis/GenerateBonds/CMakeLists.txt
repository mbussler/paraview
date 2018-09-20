CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(GenerateBonds)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

ADD_PARAVIEW_PLUGIN( GenerateBonds "1.0"
  SERVER_MANAGER_XML GenerateBonds.xml
  SERVER_MANAGER_SOURCES GenerateBonds.cxx  
)
  
INSTALL( TARGETS GenerateBonds CONFIGURATIONS Release 
                                 RUNTIME DESTINATION bin/Release 
                                 LIBRARY DESTINATION bin/Release
                                 )

INSTALL( TARGETS GenerateBonds CONFIGURATIONS Debug 
                                 RUNTIME DESTINATION bin/Debug 
                                 LIBRARY DESTINATION bin/Debug
                                 )
INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS RelWithDebInfo RUNTIME DESTINATION bin/RelWithDebInfo LIBRARY DESTINATION bin/RelWithDebInfo )