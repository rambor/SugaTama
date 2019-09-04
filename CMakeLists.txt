cmake_minimum_required(VERSION 3.8)
project(SugaTama)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -g1 -DNDEBUG")

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread -msse4.1 -I${CMAKE_CURRENT_SOURCE_DIR}/src/thirdparty")
MESSAGE("** SYSTEM : ${CMAKE_SYSTEM_NAME}")

set(CMAKE_VERBOSE_MAKEFILE ON )
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

SET (QHULL_INCLUDE_DIRS ~/src_libs/qhull-2015.2/src ~/src_libs/qhull-2015.2/src/libqhull ~/src_libs/qhull-2015.2/src/libqhullcpp)
include_directories(${QHULL_INCLUDE_DIRS})
#
# BOOST SECTION
#
# set Boost directory
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_STATIC_RUNTIME ON)
# if you have boost in a non-standard directory, like you home directory, set the following:
#/Users/robertrambo/boost/boost_1_57_0
#/Users/robertrambo/boost/boost_1_57_0/stage/lib

# Boost root variable set as environment variable in CLION Preferences
MESSAGE("** BOOST ROOT => ${BOOSTROOT}")
SET(BOOST_ROOT ${BOOSTROOT})
SET(BOOST_LIBRARYDIR ${BOOSTROOT}/stage/lib)
SET(BOOST_INCLUDEDIR ${BOOSTROOT})
include_directories(${BOOST_INCLUDEDIR})

set(BOOST_FILESYSTEM ${BOOSTROOT}/stage/lib/libboost_filesystem.a )
set(BOOST_PO ${BOOSTROOT}/stage/lib/libboost_program_options.a )
set(BOOST_RANDOM ${BOOSTROOT}/stage/lib/libboost_random.a )
set(BOOST_REGEX ${BOOSTROOT}/stage/lib/libboost_regex.a )
set(BOOST_SYSTEM ${BOOSTROOT}/stage/lib/libboost_system.a )




MESSAGE("** BOOST INCLUDE => ${BOOST_INCLUDEDIR}")
MESSAGE("** BOOST SPECIAL => ${BOOST_FILESYSTEM}")
MESSAGE("** BOOST SPECIAL => ${BOOST_PO}")
MESSAGE("** BOOST SPECIAL => ${BOOST_RANDOM}")
MESSAGE("** BOOST SPECIAL => ${BOOST_REGEX}")
MESSAGE("** BOOST SPECIAL => ${BOOST_SYSTEM}")
MESSAGE("** BOOST INCLUDE DIRECTORY: ${Boost_INCLUDE_DIR}")

MESSAGE("** CURRENT SOURCE DIRECTORY DIRECTORY: ${CMAKE_CURRENT_SOURCE_DIR}")

find_library(QHULL_LIB libqhullstatic.a PATHS ~/src_libs/qhull-2015.2/lib)
find_library(QHULLCPP_LIB libqhullcpp.a PATHS ~/src_libs/qhull-2015.2/lib)


set(SOURCE_FILES src/main.cpp)
include(src/CMakeLists.txt)
include(src/Annealer/CMakeLists.txt)
include(src/EulerTour/CMakeLists.txt)
include(src/thirdparty/CMakeLists.txt)
include(src/Aligner/CMakeLists.txt)

add_subdirectory (src)
add_subdirectory (src/base)
add_subdirectory (src/Annealer)
add_subdirectory (src/EulerTour)

set(SOURCE_FILES ${SOURCE_FILES} src/Anneal.cpp src/Model.cpp)

add_executable(sugatama ${SOURCE_FILES})
SET_TARGET_PROPERTIES(sugatama PROPERTIES LINKER_LANGUAGE CXX)

TARGET_LINK_LIBRARIES(sugatama base ${BOOST_FILESYSTEM} ${BOOST_PO} ${BOOST_RANDOM} ${BOOST_REGEX} ${BOOST_SYSTEM} ${QHULL_LIB} ${QHULLCPP_LIB})