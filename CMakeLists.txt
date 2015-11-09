cmake_minimum_required (VERSION 2.6)
project (pseudoassembler)

add_definitions(-Wall -D_FILE_OFFSET_BITS=64 )

find_package(Boost 1.46 COMPONENTS program_options system log log_setup thread REQUIRED )
include_directories(${Boost_INCLUDE_DIR} )
add_definitions(-DBOOST_LOG_DYN_LINK)

find_package(ZLIB)
include_directories(${ZLIB_INCLUDE_DIRS})

list(APPEND CMAKE_CXX_FLAGS "-std=gnu++11 -O3 -DNDEBUG ${CMAKE_CXX_FLAGS}")

add_executable(pseudoassembler main.cpp)
target_link_libraries( pseudoassembler ${ZLIB_LIBRARIES} ${Boost_LIBRARIES})
