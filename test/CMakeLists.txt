include_directories(../src)
add_executable(ContourMatcherTests ContourMatcherTests.cpp)
target_link_libraries	(ContourMatcherTests
						OCV-CM
						${OpenCV_LIBS}
						)