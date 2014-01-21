cmake_minimum_required(VERSION 2.8)
project( MasterMatcher )
find_package( OpenCV REQUIRED )
if (CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++0x")
endif()

add_subdirectory(OpenSURF/src)
add_subdirectory(src)
add_subdirectory(test)

enable_testing()
add_test( OpenCV-Test RunMatcherTests ${testdriver} )
