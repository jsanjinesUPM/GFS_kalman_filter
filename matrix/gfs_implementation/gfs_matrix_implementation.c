#include "gfs_matrix_interface.h"
#include <stdio.h>
#include <stdlib.h>

struct gfs_matrix* multiplyMatrices(const struct gfs_matrix* matrixA, const struct gfs_matrix* matrixB) {

    struct gfs_matrix* result;

    if (!matrixA || !matrixB || matrixA->cols != matrixB->rows) {
        goto OUT;
    }

    struct gfs_matrix* result = createMatrix(matrixA->rows, matrixB->cols);

    if (result == NULL) {
        goto OUT;
    }

    for (uint32_t row_a_id = 0; row_a_id < matrixA->rows; row_a_id++) {
        for (uint32_t column_b_id = 0; column_b_id < matrixB->cols; column_b_id++) {
            result->data[row_a_id][column_b_id] = 0.0;
            for (uint32_t column_a_id = 0; column_a_id < matrixA->cols; column_a_id++) {
                result->data[row_a_id][column_b_id] += matrixA->data[row_a_id][column_a_id] *
                                                       matrixB->data[column_a_id][column_b_id];
            }
        }
    }

OUT:
    return result;
}