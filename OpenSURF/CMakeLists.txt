cmake_minimum_required(VERSION 2.8)
project(OpenSURF)
find_package(OpenCV REQUIRED)

add_subdirectory(src)
add_subdirectory(test)
