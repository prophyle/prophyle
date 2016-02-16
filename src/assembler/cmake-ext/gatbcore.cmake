set(gatbcore_PREFIX ${CMAKE_BINARY_DIR}/cmake-ext/gatbcore-prefix)
set(gatbcore_INSTALL ${CMAKE_BINARY_DIR}/cmake-ext/gatbcore-install)

ExternalProject_Add(gatbcore
    PREFIX ${gatbcore_PREFIX}
    URL "http://gatb-core.gforge.inria.fr/versions/src/gatb-core-1.1.0-Source.tar.gz"
    #GIT_TAG "tag"
    INSTALL_DIR ${gatbcore_INSTALL}
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${gatbcore_INSTALL}
)

include_directories(${gatbcore_INSTALL}/include)
set(gatbcore_LIB ${gatbcore_INSTALL}/lib/libgatbcore.a)
