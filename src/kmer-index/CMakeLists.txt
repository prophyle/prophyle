cmake_minimum_required (VERSION 2.6)
project(exk)

add_executable(exk
    bwtexk.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwaseqio.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwashm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwase.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwaseqio.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwtgap.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwtaln.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bamlite.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwape.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/kopen.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/pemerge.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/maxk.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwtsw2_core.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwtsw2_main.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwtsw2_aux.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwt_lite.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwtsw2_chain.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/fastmap.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwtsw2_pair.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/utils.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/kthread.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/kstring.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/ksw.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwt.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bntseq.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwa.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwamem.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwamem_pair.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwamem_extra.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/malloc_wrap.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/QSufSort.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwt_gen.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/rope.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/rle.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/is.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bwa/bwtindex.c
)



find_package(ZLIB)
include_directories(${ZLIB_INCLUDE_DIRS})
target_link_libraries(exk ${ZLIB_LIBRARIES})

find_package (Threads)

find_library(M_LIB m)


add_custom_target(build_bwa ALL
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bwa
)
add_library(bwa STATIC IMPORTED)
set_property(TARGET bwa APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(bwa PROPERTIES IMPORTED_LOCATION_NOCONFIG "${CMAKE_CURRENT_SOURCE_DIR}/bwa/libbwa.a")
add_dependencies(bwa build_bwa)

include_directories(bwa)
#target_link_libraries(exk bwa ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${M_LIB})
target_link_libraries(exk bwa ${ZLIB_LIBRARIES} rt)

#add_dependencies(exk bwa)
#target_link_libraries(exk "bwa/")

