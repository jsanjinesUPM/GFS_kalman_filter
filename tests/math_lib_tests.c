#include "unity.h"
#include "kalman.h"
#include <string.h>
#include <stdio.h>

void setUp(void) {
}

void tearDown(void) {
}

void test_identity(void) {
    rc_matrix_t matrix = RC_MATRIX_INITIALIZER;
    rc_matrix_identity(&matrix, 2);
    TEST_ASSERT_NOT_NULL(&matrix);
    TEST_ASSERT_TRUE((matrix.d[0][0] == 1) && (matrix.d[1][1] == 1));
    TEST_ASSERT_TRUE((matrix.d[0][1] == 0) && (matrix.d[1][0] == 0));
    rc_matrix_print(matrix);
}

void test_multiplyMatrices(void) {
    rc_matrix_t a = RC_MATRIX_INITIALIZER;
    rc_matrix_t b = RC_MATRIX_INITIALIZER;
    rc_matrix_t c = RC_MATRIX_INITIALIZER;

    rc_matrix_zeros(&a, 2, 3);
    rc_matrix_zeros(&b, 3, 2);
    
    a.d[0][0] = 1;a.d[0][1] = 1;a.d[0][2] = 1;
    a.d[1][0] = 1;a.d[1][1] = 1;a.d[1][2] = 1;

    b.d[0][0] = 1;b.d[0][1] = 1;
    b.d[1][0] = 1;b.d[1][1] = 1;
    b.d[2][0] = 1;b.d[2][1] = 1;

    printf("A matrix \r\n");
    rc_matrix_print(a);
    printf("B matrix \r\n");
    rc_matrix_print(b);

    rc_matrix_multiply(a, b, &c);
    printf("Result matrix \r\n");
    rc_matrix_print(c);
}

void test_linear_kalman(void){
    
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_identity);
    RUN_TEST(test_multiplyMatrices);
    return UNITY_END();
}
