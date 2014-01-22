include_directories(../src/static)
add_executable(testexample imagedemo1.cpp)
target_link_libraries(testexample
    WAVELET-SHARED
    WAVELET-STATIC
    ${OpenCV_LIBS})
