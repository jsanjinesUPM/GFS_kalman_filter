CMAKE_MINIMUM_REQUIRED(VERSION 2.9)
PROJECT(GFS_KALMAN_FILTER)

set(MATRIX_PATH ${CMAKE_SOURCE_DIR}/matrix)
set(SUBMODULES_PATH ${CMAKE_SOURCE_DIR}/submodules)
set(TESTS_PATH ${CMAKE_SOURCE_DIR}/tests)

add_library(matrix-interface-lib INTERFACE)
target_include_directories(matrix-interface-lib INTERFACE ${MATRIX_PATH}/interface)

add_library(gfs-matrix-lib-implementation  ${MATRIX_PATH}/gfs_implementation/gfs_matrix_implementation.c)
target_link_libraries(gfs-matrix-lib-implementation PUBLIC matrix-interface-lib)

add_library(unity STATIC ${SUBMODULES_PATH}/Unity/src/unity.c)
target_include_directories(unity PUBLIC ${SUBMODULES_PATH}/Unity/src)

add_executable(matrix_tests ${TESTS_PATH}/matrix_lib_tests.c)
target_link_libraries(matrix_tests
    PUBLIC gfs-matrix-lib-implementation unity)


