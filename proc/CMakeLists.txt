cmake_minimum_required (VERSION 2.6)
project (proc)

include(ExternalProject)

include("contrib/gatbcore.cmake")
include("contrib/zlib.cmake")

add_executable(proc main.cpp)
add_dependencies(proc gatbcore zlib)
target_link_libraries(proc ${gatbcore_LIB} ${zlib_LIB})

