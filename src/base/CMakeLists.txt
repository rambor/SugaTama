cmake_minimum_required(VERSION 3.1)

PROJECT(base)

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_STATIC_RUNTIME ON)
#set(Boost_USE_MULTITHREADED ON)
# find_package(OpenMP)

#if(OPENMP_FOUND)
#    MESSAGE("** OPENMP FOUND")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set(Boost_USE_MULTITHREADED ON)
#endif()

SET (QHULL_INCLUDE_DIRS ~/qhull2/src ~/qhull2/src/libqhull_r ~/qhull2/src/libqhullcpp)
include_directories(${QHULL_INCLUDE_DIRS})

SET(BOOST_ROOT ${BOOSTROOT})
find_package(Boost 1.69.0 REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR}/boost)
    MESSAGE("** BOOST FOUND")
endif()

#find_package(LAPACK PATHS /Users/robertrambo/LAPACK/lib)
find_library(QHULL_LIB libqhullstatic_r.a PATHS ~/qhull2/lib)
find_library(QHULLCPP_LIB libqhullcpp.a PATHS ~/qhull2/lib)

find_library(BOOST_FILESYSTEM libboost_filesystem.a PATHS ~/src_libs ~/src_libs/boost_1_69_0/stage/lib)
find_library(BOOST_PO libboost_program_options.a PATHS ~/src_libs ~/src_libs/boost_1_69_0/stage/lib)
find_library(BOOST_RANDOM libboost_random.a ~/src_libs ~/src_libs/boost_1_69_0/stage/lib)
find_library(BOOST_REGEX libboost_regex.a ~/src_libs ~/src_libs/boost_1_69_0/stage/lib )
find_library(BOOST_SYSTEM libboost_system.a ~/src_libs ~/src_libs/boost_1_69_0/stage/lib)
#find_package(OpenMP)

set(SOURCE_FILES
        "${CMAKE_CURRENT_LIST_DIR}/vector3.h"
        "${CMAKE_CURRENT_LIST_DIR}/Bead.h"
        "${CMAKE_CURRENT_LIST_DIR}/Datum.h"
        "${CMAKE_CURRENT_LIST_DIR}/Data.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Objective.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Coords.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/PDBModel.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/ThreadPool.h"
        "${CMAKE_CURRENT_LIST_DIR}/functions.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Component.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/SubUnit.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Example.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Logger.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/IntraSubUnitContact.h"
        "${CMAKE_CURRENT_LIST_DIR}/CEPotential.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/CEConfiguration.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/Histogram.cpp"
        )


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2  ")

add_library(base STATIC ${SOURCE_FILES}) #  builds static lib called libbase.a in sub-directory base
TARGET_LINK_LIBRARIES(base ${BOOST_FILESYSTEM} ${BOOST_PO} ${BOOST_RANDOM} ${BOOST_REGEX} ${BOOST_SYSTEM} ${QHULL_LIB} ${QHULLCPP_LIB} )

target_include_directories (base PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

MESSAGE("** COMPILE FLAGS ${CMAKE_CXX_FLAGS}")
