CMAKE_MINIMUM_REQUIRED(VERSION 2.9)
PROJECT(GFS_KALMAN_FILTER)

set(MATRIX_PATH ${CMAKE_SOURCE_DIR}/matrix)
set(FILTER_FUNCTIONS_PATH ${CMAKE_SOURCE_DIR}/filter_functions)
set(ROBOLIB_INC ${CMAKE_SOURCE_DIR}/submodules/librobotcontrol/library/include/rc/math)
set(ROBOLIB_SRC ${CMAKE_SOURCE_DIR}/submodules/librobotcontrol/library/src/math)
set(SUBMODULES_PATH ${CMAKE_SOURCE_DIR}/submodules)
set(TESTS_PATH ${CMAKE_SOURCE_DIR}/tests)

add_library(matrix-interface-lib INTERFACE)
target_include_directories(matrix-interface-lib INTERFACE ${MATRIX_PATH}/interface)
add_library(gfs-matrix-lib-implementation  ${MATRIX_PATH}/gfs_implementation/gfs_matrix_implementation.c)
target_link_libraries(gfs-matrix-lib-implementation PUBLIC matrix-interface-lib)

add_library(unity STATIC ${SUBMODULES_PATH}/Unity/src/unity.c)
target_include_directories(unity PUBLIC ${SUBMODULES_PATH}/Unity/src)

add_library(librobotcontrol STATIC      ${ROBOLIB_SRC}/kalman.c
                                        ${ROBOLIB_SRC}/algebra_common.c
                                        ${ROBOLIB_SRC}/algebra.c
                                        ${ROBOLIB_SRC}/matrix.c
                                        ${ROBOLIB_SRC}/vector.c
                                        ${ROBOLIB_SRC}/other.c
                                        ${ROBOLIB_SRC}/filter.c
                                        ${ROBOLIB_SRC}/polynomial.c
                                        ${ROBOLIB_SRC}/quaternion.c
                                        ${ROBOLIB_SRC}/ring_buffer.c
                                        ${ROBOLIB_SRC}/kalman.c)
                                        
target_include_directories(librobotcontrol PUBLIC   ${ROBOLIB_INC})

add_library(filter_functions INTERFACE)
target_include_directories(filter_functions INTERFACE ${FILTER_FUNCTIONS_PATH}/interface)
add_library(filter_functions_implementation  ${FILTER_FUNCTIONS_PATH}/implementation/filter_functions.c)
target_link_libraries(filter_functions_implementation PUBLIC filter_functions librobotcontrol)

# add_executable(matrix_tests ${TESTS_PATH}/math_lib_tests.c)
# target_link_libraries(matrix_tests
#     PUBLIC gfs-matrix-lib-implementation unity librobotcontrol)

add_executable(olas_tests ${TESTS_PATH}/olas.c)
    target_link_libraries(olas_tests
        PUBLIC gfs-matrix-lib-implementation unity librobotcontrol filter_functions_implementation)

add_executable(filter_tests ${TESTS_PATH}/filter_tests.c)
    target_link_libraries(filter_tests
        PUBLIC gfs-matrix-lib-implementation unity librobotcontrol filter_functions_implementation)
    