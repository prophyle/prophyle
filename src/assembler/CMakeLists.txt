cmake_minimum_required (VERSION 2.6)
project (proc)


#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -lprofiler -pg")
#set(CMAKE_CXX_FLAGS " -lprofiler -pg ${CMAKE_CXX_FLAGS} -Wall ")

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'RELEASE' as none was specified.")
    set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "DEBUG" "RELEASE" "RELWITHDEBINFO")
endif()

set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS}                -std=c++11 -Wall -Wextra -Wno-missing-field-initializers -Wshadow ")
set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE}        -O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3 -march=native -DNDEBUG -g")
set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG}          -O0 -g -pedantic")



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
