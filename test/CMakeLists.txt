include_directories(../src)
add_executable(surf main.cpp)
target_link_libraries(surf OSURF-BASE ${OpenCV_LIBS})
