include_directories(../src)
include_directories(../OpenSURF/src)
add_executable(testdriver SURFTests.cpp)
target_link_libraries(testdriver
    OCV-CM
    OSURF-BASE
    ${OpenCV_LIBS})
#add_executable(ContourMatcherTests ContourMatcherTests.cpp)
#target_link_libraries(ContourMatcherTests
#    OCV-CM
#    OSURF-BASE
#    ${OpenCV_LIBS})
