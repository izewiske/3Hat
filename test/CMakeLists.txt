include_directories(../src)


add_executable(TestOnContours TestOnContours.cpp )
add_executable(TestOnContours2 TestOnContours2.cpp )

# temporarily commented out while we work on getting contours from eriol
#
#add_executable(ContourMatcherTests ContourMatcherTests.cpp eriolObjs.o)
#target_link_libraries	(ContourMatcherTests
#						OCV-CM
#						${OpenCV_LIBS}
#						)

target_link_libraries	(TestOnContours
						eriolObjs.o
						${OpenCV_LIBS}
						${OPENGL_LIBRARIES}
						${GLUT_LIBRARIES}
						)

target_link_libraries	(TestOnContours2
						eriolObjs.o
						${OpenCV_LIBS}
						${OPENGL_LIBRARIES}
						${GLUT_LIBRARIES}
						)


