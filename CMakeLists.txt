cmake_minimum_required(VERSION 2.8)
project( MasterMatcher )
find_package( OpenCV REQUIRED )
find_package( GLUT REQUIRED )
find_package( OpenGL REQUIRED )

set(_glut_libraries)
foreach(_lib ${GLUT_LIBRARIES})
	if(_lib)
		list(APPEND _glut_libraries ${_lib}) 
	endif()
endforeach()
set(GLUT_LIBRARIES ${_glut_libraries})

set(CMAKE_CXX_FLAGS "-g -Wall")

if (CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
endif()

add_definitions(-DUSE_GL)

if("${OpenCV_LIBS}" MATCHES ".*nonfree.*")
add_definitions(-DNONFREE_ENABLED)
message(STATUS "OpenCV Non-free module defined.")
message(STATUS "Non-free functions are enabled.")
endif()

link_directories(/home/cg/palantir/levmar/levmar-2.5)

if(EXISTS src/eriolObjs.o)
else()
file(	COPY /home/cg/palantir/atp14/src/eriolObjs.o
    	DESTINATION src/)
endif()

if(EXISTS src/eriolHeader.h)
else()
file(	COPY /home/cg/palantir/atp14/src/eriolHeader.h
    	DESTINATION src/)
endif()


file(MAKE_DIRECTORY build)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)

add_subdirectory(OpenSURF/src)
add_subdirectory(src)
##add_subdirectory(test)

##enable_testing()
##add_test( OpenCV-Test RunContourMatcherTests ${ContourMatcherTests} )
##add_test( OpenCV-Test RunMatcherTests ${testdriver} )
