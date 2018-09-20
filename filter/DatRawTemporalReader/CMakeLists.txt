cmake_minimum_required(VERSION 2.8)
 
PROJECT(DatRawTemporalReader)
 
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})
FIND_PATH(DATRAW_INCLUDE_PATH datRaw.h HINTS ${CMAKE_CURRENT_SOURCE_DIR})
FIND_FILE(DATRAW_LIBRARY datraw64.lib HINTS ${CMAKE_CURRENT_SOURCE_DIR})

ADD_PARAVIEW_PLUGIN(DatRawTemporalReader "1.0"
  SERVER_MANAGER_XML DatRawTemporalReader.xml
  SERVER_MANAGER_SOURCES DatRawTemporalReader.cxx)

INCLUDE_DIRECTORIES(${DATRAW_INCLUDE_PATH})
TARGET_LINK_LIBRARIES(DatRawTemporalReader PRIVATE ${DATRAW_LIBRARY} )

INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS Release 
                                 RUNTIME DESTINATION bin/Release 
                                 LIBRARY DESTINATION bin/Release
                                 )

INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS Debug 
                                 RUNTIME DESTINATION bin/Debug 
                                 LIBRARY DESTINATION bin/Debug
                                 )
INSTALL( TARGETS ${PROJECT_NAME} CONFIGURATIONS RelWithDebInfo RUNTIME DESTINATION bin/RelWithDebInfo LIBRARY DESTINATION bin/RelWithDebInfo )
                                 