cmake_minimum_required(VERSION 2.8)
project( ContourMatcher )
find_package( OpenCV REQUIRED )

file(GLOB_RECURSE OpenCVContourMatcherTestSource "*.cpp")
## file(GLOB_RECURSE OpenCVContourMatcherSource "OpenCV*.cpp") ##
## file(GLOB_RECURSE OpenCVContourMatcherTestSource "*Tests.cpp") ##
file(GLOB_RECURSE OpenCVContourMatcherTestHeaders "*.h")

## add_library(OpenCVContourMatcher SHARED ${OpenCVContourMatcherSource}) ##
add_executable( OpenCVContourMatcherTest ${OpenCVContourMatcherTestSource} )
## target_link_libraries( OpenCVContourMatcher ${OpenCV_LIBS}) ##
target_link_libraries( OpenCVContourMatcherTest ${OpenCV_LIBS} ${OpenCVContourMatcher} )