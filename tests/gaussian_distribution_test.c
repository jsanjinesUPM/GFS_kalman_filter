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
    int size = 2;
    rc_matrix_t A   = RC_MATRIX_INITIALIZER;
    rc_matrix_t L   = RC_MATRIX_INITIALIZER;
    rc_matrix_t X   = RC_MATRIX_INITIALIZER;
    rc_matrix_t R   = RC_MATRIX_INITIALIZER;

    rc_matrix_zeros(&A, size, size);
    rc_matrix_zeros(&L, size, size);
    rc_matrix_zeros(&X, size, size);
    rc_matrix_zeros(&R, size, size);
    double random_normal_x = 0;
    double random_normal_y = 0;
    double random_scaled_x = 0;
    double random_scaled_y = 0;

    A.d[0][0] = 4;
    A.d[0][1] = 12;
    //A.d[0][2] = -16;
    //A.d[0][3] = 0;

    A.d[1][0] = 12;
    A.d[1][1] = 37;
    //A.d[1][2] = -43;
    //A.d[1][3] = 3;

    // A.d[2][0] = -16;
    // A.d[2][1] = -43;
    // A.d[2][2] = 98;
    //A.d[2][3] = 0;

    // A.d[3][0] = 0;
    // A.d[3][1] = 3;
    // A.d[3][2] = 0;
    // A.d[3][3] = 9;

    //Result will be written to a file
    FILE *fpt;
    fpt = fopen("random_gaussian.csv", "w+");
    fprintf(fpt, "random_normal_x,random_normal_y,random_scaled_x,random_scaled_y\n"); //CSV file Header

    rc_ll_decomposition(A, &L);
    
    printf("\r\nMatrix L:\r\n");
    rc_matrix_print(L);

    for(int i = 0; i < 1000; i++){
        random_normal_x = rc_random_normal();
        random_normal_y = rc_random_normal();
        R.d[0][0] = random_normal_x;
        R.d[1][1] = random_normal_y;

        rc_matrix_multiply(L,R,&X);
        random_scaled_x = X.d[0][0];
        random_scaled_y = X.d[1][1];

        fprintf(fpt,"%f,%f,%f,%f\n", random_normal_x,  random_normal_y, random_scaled_x, random_scaled_y);
    }

    fclose(fpt);

    rc_matrix_free(&A);
    rc_matrix_free(&L);
}