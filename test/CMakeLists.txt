include_directories(../src)


add_executable(TestOnContours TestOnContours.cpp )
add_executable(TestOnContours2 TestOnContours2.cpp )

add_library(Eriol STATIC eriolObjs.o)


SET_TARGET_PROPERTIES(
  Eriol
  PROPERTIES
  LINKER_LANGUAGE C++
  )


# temporarily commented out while we work on getting contours from eriol
#
#add_executable(ContourMatcherTests ContourMatcherTests.cpp eriolObjs.o)
#target_link_libraries	(ContourMatcherTests
#						OCV-CM
#						${OpenCV_LIBS}
#						)

target_link_libraries	(TestOnContours
						Eriol
						${OpenCV_LIBS}
						${OPENGL_LIBRARIES}
						${GLUT_LIBRARIES}
						)

target_link_libraries	(TestOnContours2
						Eriol
						${OpenCV_LIBS}
						${OPENGL_LIBRARIES}
						${GLUT_LIBRARIES}
						)


