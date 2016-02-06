set(gatbcore_PREFIX ${CMAKE_BINARY_DIR}/contrib/gatbcore-prefix/src/gatbcore/gatb-core)
set(gatbcore_INSTALL ${CMAKE_BINARY_DIR}/contrib/gatbcore-install)

ExternalProject_Add(gatbcore
	PREFIX ${gatbcore_PREFIX}
	GIT_REPOSITORY "https://github.com/GATB/gatb-core.git"
	#GIT_TAG "fix-signed-one-bit-bitfield"
	UPDATE_COMMAND ""
	BUILD_IN_SOURCE 1
	INSTALL_DIR ${gatbcore_INSTALL}
	CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
	-DCMAKE_INSTALL_PREFIX=${gatbcore_INSTALL}
	#-DAMD64=ON
	#CMAKE_ARGS "${gtest_cmake_args}" "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
	LOG_DOWNLOAD ON
	#LOG_CONFIGURE ON
	LOG_BUILD ON 
	)


include_directories(${gatbcore_INSTALL}/include)
#add_dependencies(google_test)
#set(googletest_LIB ${gatbcore_INSTALL}/lib/libgtest.a ${googletest_INSTALL}/lib/libgtest_main.a)

