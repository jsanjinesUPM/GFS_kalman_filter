/*
* @brief for this example, we have taken inspiration from the 
* following example: https://atmos.washington.edu/~breth/classes/AS552/lect/lect26.pdf
* where the goal is to track a moving ball but we will also measure h_pos as library is not
* prepared to take a 1x1 matrix
*
*/

#include "unity.h"
#include "kalman.h"
#include "algebra.h"
#include "other.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


int main(){
    int Nx = 4;
    double time_delta = 1;
    rc_matrix_t D   = RC_MATRIX_INITIALIZER;
    rc_matrix_t F   = RC_MATRIX_INITIALIZER;
    rc_vector_t V   = RC_VECTOR_INITIALIZER;
    rc_vector_t U   = RC_VECTOR_INITIALIZER;
    rc_matrix_identity(&D, Nx);
    rc_matrix_zeros(&F, Nx, Nx);
    rc_vector_zeros(&V, Nx);
    rc_vector_zeros(&U, Nx);

    // D.d[0][2] = 1;
    // D.d[1][3] = 1;

    V.d[3] = -9.8;

    rc_get_state_matrix(Nx, D, V, time_delta, &F, &U);
    printf("\r\nMatrix F:\r\n");
    rc_matrix_print(F);
    printf("\r\nVector U:\r\n");
    rc_vector_print(U);
}