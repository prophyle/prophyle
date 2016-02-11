cmake_minimum_required (VERSION 2.6)
project (proc)


add_definitions(" --std=c++11")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -lprofiler -pg")
#set(CMAKE_CXX_FLAGS " -lprofiler -pg ${CMAKE_CXX_FLAGS} -Wall ")
set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wall -O3 -march=native")

add_executable(proc main.cpp)

include(ExternalProject)

include("cmake-ext/klib.cmake")
add_dependencies(proc klib)

#include("cmake-ext/gatbcore.cmake")
#add_dependencies(proc gatbcore)
#target_link_libraries(proc ${gatbcore_LIB})

include("cmake-ext/zlib.cmake")
add_dependencies(proc zlib)
target_link_libraries(proc ${zlib_LIB})

find_package( Boost 1.46 COMPONENTS program_options REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )
target_link_libraries(proc ${Boost_LIBRARIES})
