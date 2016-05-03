cmake_minimum_required (VERSION 2.6)
project(exk)

add_executable(exk bwtexk.c)

find_package(ZLIB)
include_directories(${ZLIB_INCLUDE_DIRS})
target_link_libraries(exk ${ZLIB_LIBRARIES})


add_custom_target(build_bwa ALL
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bwa
)
add_library(bwa STATIC IMPORTED)
set_property(TARGET bwa APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(bwa PROPERTIES IMPORTED_LOCATION_NOCONFIG "${CMAKE_CURRENT_SOURCE_DIR}/bwa/libbwa.a")
add_dependencies(bwa build_bwa)

include_directories(bwa)
target_link_libraries(exk bwa ${ZLIB_LIBRARIES})

#add_dependencies(exk bwa)
#target_link_libraries(exk "bwa/")

