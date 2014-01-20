include_directories(../src)
add_executable(ContourMatcherTests ContourMatcherTests.cpp eriolObjs.o)
add_executable(TestOnContours TestOnContours.cpp eriolObjs.o)
target_link_libraries	(ContourMatcherTests
						OCV-CM
						${OpenCV_LIBS}
						)

target_link_libraries	(TestOnContours
						${OpenCV_LIBS}
						)