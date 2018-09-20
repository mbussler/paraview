cmake_minimum_required(VERSION 2.8)

project(TetGen)

#linalg
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../common/)

#tetlib
add_definitions( -DTETLIBRARY )
#add_definitions( -DSINGLE )
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}/Tetgen" )
link_directories( "${CMAKE_CURRENT_SOURCE_DIR}/Tetgen/lib")
#add_subdirectory( Tetgen )


IF (ParaView_SOURCE_DIR)
  INCLUDE_DIRECTORIES(
    ${VTK_INCLUDE_DIRS}
  )
ELSE (ParaView_SOURCE_DIR)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF (ParaView_SOURCE_DIR)

ADD_PARAVIEW_PLUGIN(TetGen "0.1"
  SERVER_MANAGER_XML TetGen.xml
  SERVER_MANAGER_SOURCES vtkTetGen.cxx
)

target_link_libraries( TetGen PRIVATE tetlib )

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
                                 LIBRARY DESTINATION bin/RelWithDebInfo )
