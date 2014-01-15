include_directories(../src)
include_directories(../OpenSURF/src)
add_executable(testdriver main.cpp)
target_link_libraries(testdriver 3SURF OSURF-BASE ${OpenCV_LIBS})
