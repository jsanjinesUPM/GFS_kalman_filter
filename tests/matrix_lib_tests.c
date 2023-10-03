#include "unity.h"
#include "gfs_matrix_interface.h"
#include <string.h>
#include <stdio.h>

void setUp(void) {
}

void tearDown(void) {
}

void test_createMatrix(void) {
    struct gfs_matrix* matrix = createMatrix(2, 3);
    TEST_ASSERT_NOT_NULL(matrix);
    freeMatrix(matrix);
}

void test_multiplyMatrices(void) {
    struct gfs_matrix* matrixA = createMatrix(3, 2);
    TEST_ASSERT_NOT_NULL(matrixA);

    double matrixA_data[3][2] = {
        {1, 2},
        {4, 5},
        {7, 8}
    };

    for (uint32_t row = 0; row < matrixA->rows; row++)
        for (uint32_t column = 0; column < matrixA->cols; column++)  
            matrixA->data[row][column] = matrixA_data[row][column];
    

    printMatrix(matrixA);

    struct gfs_matrix *matrixB = createMatrix(2, 3);
    TEST_ASSERT_NOT_NULL(matrixB);

    double matrixB_data[2][3] = {{1, 2, 3},
                                 {0, 5, 2}};

    for (uint32_t row = 0; row < matrixB->rows; row++)
        for (uint32_t column = 0; column < matrixB->cols; column++)  
            matrixB->data[row][column] = matrixB_data[row][column];

    printMatrix(matrixB);

    struct gfs_matrix *result = multiplyMatrices(matrixA, matrixB);

    TEST_ASSERT_NOT_NULL(result);
    printMatrix(result);
    freeMatrix(matrixA);
    freeMatrix(matrixB);
    freeMatrix(result);
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_createMatrix);
    RUN_TEST(test_multiplyMatrices);
    return UNITY_END();
}
