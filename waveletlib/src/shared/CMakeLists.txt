file(GLOB_RECURSE WAVELET-SHARED-Source *.cpp)

add_library(WAVELET-SHARED
    ${WAVELET-SHARED-Source})
