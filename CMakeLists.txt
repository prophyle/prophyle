cmake_minimum_required (VERSION 2.6)
project (pseudoassembler)

FIND_PACKAGE( Boost 1.46 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

find_package(ZLIB)
include_directories(${ZLIB_INCLUDE_DIRS})

list( APPEND CMAKE_CXX_FLAGS "-std=gnu++11  -O3 -DNDEBUG ${CMAKE_CXX_FLAGS}")
add_definitions( -Wall -D_FILE_OFFSET_BITS=64 )

add_executable( pseudoassembler main.cpp)
target_link_libraries( pseudoassembler ${ZLIB_LIBRARIES} ${Boost_LIBRARIES})

