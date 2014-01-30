include_directories(../src)
include_directories(../OpenSURF/src)


add_executable(TestOnContours TestOnContours.cpp  eriolObjs.o)
target_link_libraries(	TestOnContours
						OSURF-BASE
						${OpenCV_LIBS}
						${OPENGL_LIBRARIES}
						${GLUT_LIBRARIES})

add_executable(testdriver SURFTests.cpp)

target_link_libraries(	testdriver
    					OSURF-BASE
   						${OpenCV_LIBS})