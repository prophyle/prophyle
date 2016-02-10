set(klib_PREFIX ${CMAKE_BINARY_DIR}/cmake-ext/klib-prefix)
set(klib_INSTALL ${CMAKE_BINARY_DIR}/cmake-ext/klib-install)

ExternalProject_Add(klib
    PREFIX ${klib_PREFIX}
    GIT_REPOSITORY "https://github.com/attractivechaos/klib"
    #GIT_TAG "tag"
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory
        <SOURCE_DIR> ${klib_INSTALL}
)

include_directories(${klib_INSTALL})
