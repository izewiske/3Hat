file(GLOB_RECURSE OCV-CMSource *.cpp)

#if(EXISTS "${ROOT}/src/eriolObjs.o" AND "${ROOT}/src/eriolHeader.h")
#else()
include_directories(/home/cg/palantir/atp14/src)
#endif()

add_executable(MultiContourDriver MultipleContourDriver.cpp eriolObjs.o)

include_directories(../OpenSURF/src)
target_link_libraries	(MultiContourDriver
						OSURF-BASE
						${OpenCV_LIBS}
						${OPENGL_LIBRARIES}
						${GLUT_LIBRARIES}
						)