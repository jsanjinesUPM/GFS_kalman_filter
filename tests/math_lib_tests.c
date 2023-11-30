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

int main(void) {

    //Result will be written to a file
    FILE *fpt;
    fpt = fopen("kalman_results.csv", "w+");
    fprintf(fpt, "time,h_position,v_position,h_velocity,v_velocity,measurement,"
                    "pi_h_pos,pi_v_pos,pi_h_vel,pi_v_vel,real_h_pos,real_v_pos,"
                    "real_h_vel,real_v_vel\n"); //CSV file Header

    int Nx = 4;
    int Ny = 1;
    int Nu = 1;
    double time_delta = 0.1;

    // declare variables
    rc_kalman_t kfilter  = RC_KALMAN_INITIALIZER;
    rc_matrix_t phi = RC_MATRIX_INITIALIZER;
    rc_matrix_t F   = RC_MATRIX_INITIALIZER;
    rc_matrix_t D   = RC_MATRIX_INITIALIZER;
    rc_matrix_t G   = RC_MATRIX_INITIALIZER;
    rc_vector_t u   = RC_VECTOR_INITIALIZER;
    rc_vector_t U   = RC_VECTOR_INITIALIZER;
    rc_vector_t V   = RC_VECTOR_INITIALIZER;
    rc_matrix_t H   = RC_MATRIX_INITIALIZER;
    rc_matrix_t Q   = RC_MATRIX_INITIALIZER;
    rc_matrix_t L   = RC_MATRIX_INITIALIZER; // Used for LL decomposition of Q
    rc_matrix_t M   = RC_MATRIX_INITIALIZER; //Matrix for random normal data
    rc_matrix_t R   = RC_MATRIX_INITIALIZER;
    rc_matrix_t Pi  = RC_MATRIX_INITIALIZER;
    rc_vector_t y   = RC_VECTOR_INITIALIZER;
    rc_vector_t h   = RC_VECTOR_INITIALIZER;
    rc_vector_t x_0   = RC_VECTOR_INITIALIZER;
    rc_vector_t new_x   = RC_VECTOR_INITIALIZER;
    // allocate appropriate memory for system
    rc_matrix_zeros(&F, Nx, Nx);
    rc_matrix_zeros(&G, Nx, Nx);
    rc_matrix_zeros(&D, Nx, Nx);
    rc_matrix_zeros(&H, Ny, Nx);
    rc_matrix_zeros(&Q, Nx, Nx);
    rc_matrix_zeros(&L, Nx, Nx);
    rc_matrix_zeros(&M, Nx, Nx);
    rc_matrix_zeros(&R, Ny, Ny);
    rc_matrix_zeros(&Pi, Nx, Nx);
    rc_vector_zeros(&y, Ny);
    rc_vector_zeros(&h, Ny);
    rc_vector_zeros(&U, Nx);
    rc_vector_zeros(&V, Nx);
    rc_vector_zeros(&u, Nx);
    rc_vector_zeros(&x_0, Nx);
    rc_vector_zeros(&new_x, Nx);

    //Give matrices inicial value
    x_0.d[0]= 0;
    x_0.d[1]= 50;
    x_0.d[2]= 20;
    x_0.d[3]= 50;

    rc_matrix_identity(&F, F.cols);
    D.d[0][2] = 1;
    D.d[1][3] = 1;

    V.d[3] = -9.8;

    R.d[0][0] = 0.00;


    //Manually adding noise to predicted (and real) horizontal pos and horizontal Vel
    Q.d[0][0] = 0.005;
    Q.d[1][3] = 0.005;
    Q.d[1][1] = 0.005;
    Q.d[2][2] = 0.01;
    Q.d[3][1] = 0.005;
    Q.d[3][3] = 0.01;

    rc_ll_decomposition(Q, &L);

    Pi.d[0][0] = 0.1;
    Pi.d[1][1] = 0.1;
    Pi.d[2][2] = 2;
    Pi.d[3][3] = 2;

    //We need to calculate "real" values to generate measurements and compare results
    double real_h_pos = 0.01;
    double real_v_pos = 0.01;
    double real_h_vel = 30;
    double real_v_vel = 60;

    double measured_angle = 0;
    double measured_v = 0;

    //Variables to update H
    double h_pos = 0;
    double v_pos = 0;

    rc_get_state_matrix(Nx, D, V, time_delta, &F, &U);
    if(rc_kalman_new_alloc(&kfilter, F, G, H, Q, R, Pi, u, x_0)==-1) return -1;

    for(int i = 0; i <= 1000; i++){
        if(rc_kalman_predict_ekf(&kfilter, F, G, U, u) == -1) return -1;

        if(i % 100 == 0){
            y.d[0] = atan2(real_v_pos, real_h_pos) + rc_random_normal()*sqrt(R.d[0][0]);
            // Update H
            h_pos = kfilter.x_pre.d[0];
            v_pos = kfilter.x_pre.d[1];
            H.d[0][0] = -v_pos/((h_pos*h_pos)+(v_pos*v_pos));
            H.d[0][1] = h_pos/((h_pos*h_pos)+(v_pos*v_pos));
            // h = h(x_pre) --> h = tan(y/x) for our particular example
            h.d[0] = atan2(v_pos, h_pos);
            if(rc_kalman_prediction_update_ekf(&kfilter, H, y, h) == -1) return -1;

            measured_angle = y.d[0];
            //measured_v = y.d[1];
        }
        fprintf(fpt,"%f, %f, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n",
                                        time_delta*i ,kfilter.x_pre.d[0], kfilter.x_pre.d[1], kfilter.x_pre.d[2],
                                        kfilter.x_pre.d[3], y.d[0], kfilter.P.d[0][0], kfilter.P.d[1][1],
                                        kfilter.P.d[2][2], kfilter.P.d[3][3], real_h_pos, real_v_pos, 
                                        real_h_vel, real_v_vel);

        //Update Real Values
        rc_random_normal(&M);
        rc_matrix_left_multiply_inplace(L,&M);
        real_h_pos = 0 + real_h_vel*time_delta*i + M.d[0][0];
        real_v_pos = 0 + 50*time_delta*i - 9.8*0.5*time_delta*time_delta*i*i + M.d[1][1];
        real_h_vel = 30 + M.d[2][2];
        real_v_vel = 60 - 9.8*time_delta*i + M.d[3][3];

        //Reset Measurements for writing to file purposes
        measured_angle = 0;
        measured_v = 0;
    }

    fclose(fpt);
    //Free memmory
    rc_matrix_free(&F);
    rc_matrix_free(&G);
    rc_matrix_free(&H);
    rc_matrix_free(&Q);
    rc_matrix_free(&R);
    rc_matrix_free(&D);
    rc_matrix_free(&Pi);
    rc_vector_free(&V);
    rc_vector_free(&h);
    rc_vector_free(&y);
    rc_vector_free(&U);
    rc_vector_free(&new_x);
    rc_vector_free(&x_0);
    rc_kalman_free(&kfilter);
    return 0;
}