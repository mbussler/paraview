CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(ImplicitPolydataClip)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

ADD_PARAVIEW_PLUGIN( ImplicitPolydataClip "1.0"
  SERVER_MANAGER_XML ImplicitPolydataClip.xml
  SERVER_MANAGER_SOURCES ImplicitPolydataClip.cxx  
)
  
INSTALL( TARGETS ImplicitPolydataClip CONFIGURATIONS Release 
                                 RUNTIME DESTINATION bin/Release 
                                 LIBRARY DESTINATION bin/Release
                                 )

INSTALL( TARGETS ImplicitPolydataClip CONFIGURATIONS Debug 
                                 RUNTIME DESTINATION bin/Debug 
                                 LIBRARY DESTINATION bin/Debug
                                 )
INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS RelWithDebInfo RUNTIME DESTINATION bin/RelWithDebInfo LIBRARY DESTINATION bin/RelWithDebInfo )