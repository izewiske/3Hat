include_directories(../src)
include_directories(../OpenSURF/src)


add_executable(TestOnContours TestOnContours.cpp  ../src/eriolObjs.o)
target_link_libraries(	TestOnContours
						OSURF-BASE
						${OpenCV_LIBS}
						${OPENGL_LIBRARIES}
						${GLUT_LIBRARIES})

add_executable(SurfTestDriver SURFTests.cpp)

target_link_libraries(	SurfTestDriver
	    					OSURF-BASE
   						${OpenCV_LIBS})
