CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(EigenValueTest)

# VTK
#FIND_PACKAGE(VTK REQUIRED COMPONENTS vtkCommonCore)
#include(${VTK_USE_FILE})

# linalg.h
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# find and include eigen
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

ADD_EXECUTABLE( EigenValueTest main.cpp )

target_link_libraries(EigenValueTest ${VTK_LIBRARIES})