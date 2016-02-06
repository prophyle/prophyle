cmake_minimum_required (VERSION 2.6)
project (proc)

include(ExternalProject)

include("contrib/gatbcore.cmake")

add_executable(proc main.cpp)
add_dependencies(proc gatbcore)

