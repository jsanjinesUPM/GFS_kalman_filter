/*
* @brief Test for our function that generates a random number following 
*       a multivariate gaussian distribution
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
    rc_matrix_t A   = RC_MATRIX_INITIALIZER;
    rc_matrix_t L   = RC_MATRIX_INITIALIZER;
    rc_matrix_t Diagonal   = RC_MATRIX_INITIALIZER;
    rc_matrix_t D   = RC_MATRIX_INITIALIZER;
    rc_matrix_t LD   = RC_MATRIX_INITIALIZER;
    rc_matrix_t LT   = RC_MATRIX_INITIALIZER;
    rc_matrix_zeros(&A, 3, 3);
    rc_matrix_zeros(&L, 3, 3);
    rc_matrix_zeros(&Diagonal, 3, 3);
    rc_matrix_zeros(&D, 3, 3);
    double random_normal = 0;
    double random_scaled = 0;

    A.d[0][0] = 4;
    A.d[0][1] = 12;
    A.d[0][2] = -16;

    A.d[1][0] = 12;
    A.d[1][1] = 37;
    A.d[1][2] = -43;

    A.d[2][0] = -16;
    A.d[2][1] = -43;
    A.d[2][2] = 98;

    //Result will be written to a file
    FILE *fpt;
    fpt = fopen("random_gaussian.csv", "w+");
    fprintf(fpt, "random_normal,random_scaled\n"); //CSV file Header

    rc_ldl_decomposition(A, &L, &Diagonal);

    for(int i = 0; i < 100000; i++){
        random_normal = rc_random_normal();
        rc_matrix_duplicate(Diagonal, &D);
        rc_matrix_times_scalar(&D, random_normal);
        rc_matrix_multiply(L, D, &LD);
        rc_matrix_transpose(L, &LT);
        rc_matrix_left_multiply_inplace(LD, &LT);
        random_scaled = LT.d[0][0];

        fprintf(fpt,"%f,%f,\n", random_normal, random_scaled);
    }

    fclose(fpt);

    rc_matrix_free(&A);
    rc_matrix_free(&Diagonal);
    rc_matrix_free(&L);
    rc_matrix_free(&LD);
    rc_matrix_free(&LT);
}