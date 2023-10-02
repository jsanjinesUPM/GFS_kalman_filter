#include "gfs_matrix_interface.h"
#include <stdio.h>
#include <stdlib.h>

struct gfs_matrix* multiplyMatrices(const struct gfs_matrix* matrixA, const struct gfs_matrix* matrixB) {

    struct gfs_matrix* result;

    if (!matrixA || !matrixB || matrixA->cols != matrixB->rows) {
        goto OUT;
    }

    result = createMatrix(matrixA->rows, matrixB->cols);

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

struct gfs_matrix* createMatrix(uint32_t rows, uint32_t cols) {
    
    struct gfs_matrix* matrix = (struct gfs_matrix*)malloc(sizeof(struct gfs_matrix));

    if (!matrix) 
        goto OUT; 
    
    matrix->rows = rows;
    matrix->cols = cols;
 
    matrix->data = (double**)malloc(rows * sizeof(double*));

    if (matrix->data == NULL) {
        free(matrix);
        matrix = NULL;
        goto OUT; 
    }

    for (uint32_t row = 0; row < rows; row++) {
        matrix->data[row] = (double*)malloc(cols * sizeof(double));

        if (matrix->data[row] == NULL) {
            
            for (uint32_t j = 0; j < row; j++) 
                free(matrix->data[j]);
            
            free(matrix->data);
            free(matrix);
            matrix = NULL;
            goto OUT;;
        }
    }

OUT:
    return matrix;
}

void freeMatrix(struct gfs_matrix* matrix) {
    if (!matrix)
        return;

    for (uint32_t i = 0; i < matrix->rows; i++) {
        if (matrix->data[i])
            free(matrix->data[i]);
    }

    free(matrix->data);
    free(matrix);
}

void printMatrix(const struct gfs_matrix* matrix) {
    if (!matrix) {
        printf("Matrix not initialized\n");
        return;
    }

    printf("Matrix (%u x %u):\n", matrix->rows, matrix->cols);

    for (uint32_t row = 0; row < matrix->rows; row++) {
        for (uint32_t column = 0; column < matrix->cols; column++) {
            printf("%f ", matrix->data[row][column]);
        }
        printf("\n");
    }
}
