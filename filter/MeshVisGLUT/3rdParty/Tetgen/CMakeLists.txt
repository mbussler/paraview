cmake_minimum_required (VERSION 2.6)

set(CMAKE_DEBUG_POSTFIX _d)

set_property(SOURCE predicates.cxx PROPERTY COMPILE_FLAGS "-O0")
set_property(SOURCE tetgen.cxx PROPERTY COMPILE_FLAGS "-Wall -DSELF_CHECK")

project( tetgen )
add_executable( tetgen tetgen.cxx predicates.cxx )
 
project( tetlib )
add_library( tetlib STATIC tetgen.cxx predicates.cxx )
set_target_properties (tetlib PROPERTIES COMPILE_DEFINITIONS "TETLIBRARY" )

install( TARGETS tetlib ARCHIVE DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib )
install( TARGETS tetgen RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin )
