
#include "kalman.h"
#include "algebra.h"
#include "other.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void write_header(FILE *writefile, double time, int x_size, int y_size);
void write_to_file(FILE *writefile, double time, rc_vector_t x, rc_vector_t h, rc_vector_t y);
int mesh_size(const char* meshfile_path);
void load_mesh(rc_matrix_t* K, const char* meshfile_path);
int get_w(rc_vector_t* w, rc_matrix_t K, double depth);
int get_H(rc_matrix_t* H, rc_matrix_t K, rc_matrix_t X_Y, rc_vector_t w, double time);
int add_H_prime(rc_matrix_t* N, rc_matrix_t H, rc_matrix_t K, rc_matrix_t X_Y, rc_vector_t w, double time);
