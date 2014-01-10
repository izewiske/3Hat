cmake_minimum_required(VERSION 2.8)
project( ContourMatcher )
find_package( OpenCV REQUIRED )

add_subdirectory(src)
add_subdirectory(test)
enable_testing()
add_test( OpenCV-Test RunContourMatcherTests ${ContourMatcherTests} )