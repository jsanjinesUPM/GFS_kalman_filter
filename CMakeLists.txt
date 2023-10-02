CMAKE_MINIMUM_REQUIRED(VERSION 2.9)
PROJECT(GFS_KALMAN_FILTER)

set(MATRIX_PATH ${CMAKE_SOURCE_DIR}/matrix)

add_library(matrix-interface-lib INTERFACE)
target_include_directories(matrix-interface-lib INTERFACE ${MATRIX_PATH}/interface)

add_library(gfs-matrix-lib-implementation  ${MATRIX_PATH}/gfs_implementation/gfs_matrix_implementation.c)
target_link_libraries(gfs-matrix-lib-implementation PUBLIC matrix-interface-lib)