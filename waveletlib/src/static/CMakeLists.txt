file(GLOB_RECURSE WAVELET-STATIC-Source *.cpp)

add_library(WAVELET-STATIC
    ${WAVELET-STATIC-Source})
