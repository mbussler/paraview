cmake_minimum_required(VERSION 2.8)
 
PROJECT(LAMMPSTemporalReaderMF)
 
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

ADD_PARAVIEW_PLUGIN( LAMMPSTemporalReaderMF "1.0"
  SERVER_MANAGER_XML LAMMPSTemporalReaderMF.xml
  SERVER_MANAGER_SOURCES LAMMPSTemporalReaderMF.cxx
  REQUIRED_ON_SERVER)

TARGET_LINK_LIBRARIES( LAMMPSTemporalReaderMF PRIVATE pqCore PeriMultiFitTaylor )
  
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
