cmake_minimum_required(VERSION 2.8)
project( ContourMatcher )
find_package( OpenCV REQUIRED )
find_package( GLUT REQUIRED )
find_package( OpenGL REQUIRED )

set(_glut_libraries)
foreach(_lib ${GLUT_LIBRARIES})
	if(_lib)
		list(APPEND _glut_libraries)
	endif()
endforeach()
set(GLUT_LIBRARIES ${_glut_libraries})


if (CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
endif()

set(CMAKE_CXX_FLAGS "-g -Wall ")

add_subdirectory(src)
add_subdirectory(test)
enable_testing()
add_test( OpenCV-Test RunContourMatcherTests ${ContourMatcherTests} )
