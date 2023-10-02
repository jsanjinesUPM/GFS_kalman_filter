#ifndef GFS_MATRIX_INTERFACE_H
#define GFS_MATRIX_INTERFACE_H

#include <stdint.h>

struct gfs_matrix {
    uint32_t rows;
    uint32_t cols;     
    double** data;
};

struct gfs_matrix* createMatrix(uint32_t rows, uint32_t cols);
void freeMatrix(struct gfs_matrix* matrix);
struct gfs_matrix* multiplyMatrices(const struct gfs_matrix* matrixA, const struct gfs_matrix* matrixB);
void printMatrix(const struct gfs_matrix* matrix);

#endif /* GFS_MATRIX_INTERFACE_H */
