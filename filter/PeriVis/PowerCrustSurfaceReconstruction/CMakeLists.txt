CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(PowerCrustSurfaceReconstruction)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

ADD_PARAVIEW_PLUGIN( PowerCrustSurfaceReconstruction "1.0"
  SERVER_MANAGER_XML PowerCrustSurfaceReconstruction.xml
  SERVER_MANAGER_SOURCES vtkPowerCrustSurfaceReconstruction.cxx  
)
  
INSTALL( TARGETS PowerCrustSurfaceReconstruction CONFIGURATIONS Release 
                                 RUNTIME DESTINATION bin/Release 
                                 LIBRARY DESTINATION bin/Release
                                 )

INSTALL( TARGETS PowerCrustSurfaceReconstruction CONFIGURATIONS Debug 
                                 RUNTIME DESTINATION bin/Debug 
                                 LIBRARY DESTINATION bin/Debug
                                 )
INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS RelWithDebInfo 
                                RUNTIME DESTINATION bin/RelWithDebInfo 
                                LIBRARY DESTINATION bin/RelWithDebInfo 
                                )