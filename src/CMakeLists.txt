file(GLOB_RECURSE OCV-CMSource *.cpp)

#set(eriol_object_file ${CMAKE_CURRENT_BINARY_DIR}/../eriol/eriolObjs.o)
#set_source_files_properties(
#   ${eriol_object_file}
#   PROPERTIES
#   EXTERNAL_OBJECT true #only link; do not compile
#)

add_library(OCV-CM
#    ${eriol_object_file}
    ${OCV-CMSource} )
