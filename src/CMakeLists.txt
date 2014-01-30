file(GLOB_RECURSE OCV-CMSource *.cpp)

add_executable(MultiContourDriver MultipleContourDriver.cpp eriolObjs.o)

include_directories(../OpenSURF/src)
target_link_libraries	(MultiContourDriver
						OSURF-BASE
						${OpenCV_LIBS}
						${OPENGL_LIBRARIES}
						${GLUT_LIBRARIES}
						)