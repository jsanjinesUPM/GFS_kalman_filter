#include "unity.h"
#include "filter_functions.h"
#include "kalman.h"
#include "algebra.h"
#include "other.h"
#include <string.h>
#include <stdio.h>

void setUp(void) {
}

void tearDown(void) {
}

void test_read_mesh(void) {
    const char meshfile_path[] = "C:/Users/mateo/GFS_kalman_filter/tests/test_csv/test_mesh.csv";
    TEST_ASSERT_NOT_NULL(meshfile_path);

    int K_rows = mesh_size(meshfile_path);
    TEST_ASSERT(K_rows == 5);

    rc_matrix_t K = RC_MATRIX_INITIALIZER;
    rc_matrix_zeros(&K, K_rows, 2);
    load_mesh(&K, meshfile_path);

    TEST_ASSERT(K.d[0][0] == 0.7969557584733965);
    TEST_ASSERT(K.d[3][1] == 0.7013814045660062);
    rc_matrix_free(&K);
}

void test_get_w(void) {
    int k_rows = 2;
    rc_matrix_t K = RC_MATRIX_INITIALIZER;
    rc_matrix_zeros(&K, k_rows, 2);
    K.d[0][0] = 0.5; K.d[0][1] = 0.7;
    K.d[1][0] = 1.23; K.d[1][1] = 1.56;

    double depth = 20;
    
    rc_vector_t w = RC_VECTOR_INITIALIZER;
    rc_vector_zeros(&w, k_rows);

    get_w(&w, K, depth);
    TEST_ASSERT(round(w.d[0]*1000) == 2903);
    TEST_ASSERT(round(w.d[1]*1000) == 4412);

    rc_matrix_free(&K);
    rc_vector_free(&w);
}

void test_get_H(void) {
    int k_rows = 3;
    rc_matrix_t K = RC_MATRIX_INITIALIZER;
    rc_matrix_t H = RC_MATRIX_INITIALIZER;
    rc_matrix_t X_Y = RC_MATRIX_INITIALIZER;
    rc_matrix_zeros(&H, 3, k_rows*8);
    rc_matrix_zeros(&K, k_rows, 2);
    K.d[0][0] = 0.5; K.d[0][1] = 0.7;
    K.d[1][0] = 2; K.d[1][1] = 3;
    K.d[2][0] = 3; K.d[2][1] = 2;

    rc_matrix_zeros(&X_Y, 3, 2);
    X_Y.d[0][0] = 0; X_Y.d[0][1] = 0;
    X_Y.d[1][0] = 1; X_Y.d[1][1] = 3;
    X_Y.d[2][0] = 3; X_Y.d[2][1] = 2;

    double depth = 20;
    double time = 3.45;
    
    rc_vector_t w = RC_VECTOR_INITIALIZER;
    rc_vector_zeros(&w, k_rows);

    get_w(&w, K, depth);
    TEST_ASSERT(round(w.d[0]*1000) == 2903);
    TEST_ASSERT(round(w.d[1]*1000) == 5944);
    TEST_ASSERT(round(w.d[2]*1000) == 5944);

    get_H(&H, K, X_Y, w, time);
    TEST_ASSERT(round(H.d[0][0]*1000) == 0);
    TEST_ASSERT(round(H.d[0][6]*1000) == -558);
    TEST_ASSERT(round(H.d[2][14]*1000) == 918);


    rc_matrix_free(&K);
    rc_matrix_free(&X_Y);
    rc_vector_free(&w);
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_read_mesh);
    RUN_TEST(test_get_w);
    RUN_TEST(test_get_H);
    return UNITY_END();
}
