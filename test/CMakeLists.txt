include_directories(../src)


add_executable(TestOnContours TestOnContours.cpp eriolObjs.o showTile.cpp)

# temporarily commented out while we work on getting contours from eriol
#
#add_executable(ContourMatcherTests ContourMatcherTests.cpp eriolObjs.o)
#target_link_libraries	(ContourMatcherTests
#						OCV-CM
#						${OpenCV_LIBS}
#						)

target_link_libraries	(TestOnContours
						${OpenCV_LIBS}
						${OPENGL_LIBRARIES}
						${GLUT_LIBRARIES}
						)