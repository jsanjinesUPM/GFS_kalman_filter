/*
* @brief for this example, we have taken inspiration from the 
* following example: https://atmos.washington.edu/~breth/classes/AS552/lect/lect26.pdf
* where the goal is to track a moving ball but we will also measure h_pos as library is not
* prepared to take a 1x1 matrix
*
*/

#include "unity.h"
#include "kalman.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {

    //Result will be written to a file
    FILE *fpt;
    fpt = fopen("kalman_results.csv", "w+");
    fprintf(fpt, "time,h_position,v_position,h_velocity,v_velocity,measurement,"
                    "pi_h_pos,pi_v_pos,pi_h_vel,pi_v_vel,real_h_pos,real_v_pos,"
                    "real_h_vel,real_v_vel,measured_angle,measured_v\n"); //CSV file Header

    int Nx = 4;
    int Ny = 2;
    int Nu = 2;
    double time_delta = 0.1;

    // declare variables
    rc_kalman_t kf  = RC_KALMAN_INITIALIZER;
    rc_matrix_t F   = RC_MATRIX_INITIALIZER;
    rc_matrix_t G   = RC_MATRIX_INITIALIZER;
    rc_matrix_t H   = RC_MATRIX_INITIALIZER;
    rc_matrix_t Q   = RC_MATRIX_INITIALIZER;
    rc_matrix_t R   = RC_MATRIX_INITIALIZER;
    rc_matrix_t Pi  = RC_MATRIX_INITIALIZER;
    rc_vector_t u   = RC_VECTOR_INITIALIZER;
    rc_vector_t y   = RC_VECTOR_INITIALIZER;
    rc_vector_t x_0   = RC_VECTOR_INITIALIZER;
    // allocate appropriate memory for system
    rc_matrix_zeros(&F, Nx, Nx);
    rc_matrix_zeros(&G, Nx, Nu);
    rc_matrix_zeros(&H, Ny, Nx);
    rc_matrix_zeros(&Q, Nx, Nx);
    rc_matrix_zeros(&R, Ny, Ny);
    rc_matrix_zeros(&Pi, Nx, Nx);
    rc_vector_zeros(&u, Nu);
    rc_vector_zeros(&y, Ny);
    rc_vector_zeros(&x_0, Nx);

    //Give matrices inicial value
    x_0.d[0]= 0;
    x_0.d[1]= 0;
    x_0.d[2]= 30;
    x_0.d[3]= 50;

    u.d[0] = 9.8;

    rc_matrix_identity(&F, F.cols);
    F.d[0][2] = time_delta;
    F.d[1][3] = time_delta;
                    //      [-y/x^2+y^2  x/x^2+y^2   0   0]
    H.d[1][1] = 1;  // H =  [     0          1       0   0]

    G.d[1][0] = -time_delta*time_delta/2;   //      [0      0]
    G.d[3][0] = -time_delta;                // G =  [-t^2/2 0]
                                            //      [0      0]
    // R.d[0][0] = 0.05;                       //      [-t     0]
    // R.d[0][1] = 0.5;
    // R.d[1][1] = 0.5;
    // R.d[1][1] = 5;
    R.d[0][0] = 0.0001;
    R.d[0][1] = 0.0001;
    R.d[1][1] = 0.0001;
    R.d[1][1] = 0.0001;

    rc_matrix_t tmp1 = RC_MATRIX_INITIALIZER;
    rc_matrix_transpose(G, &tmp1);
    rc_matrix_multiply(G, tmp1, &Q);
    rc_matrix_times_scalar(&Q,100);
    rc_matrix_free(&tmp1);

    Pi.d[0][0] = 0.1;
    Pi.d[0][0] = 0.1;
    Pi.d[0][0] = 2;
    Pi.d[0][0] = 2;

    //We need to calculate "real" values to generate measurements and compare results
    double real_h_pos = 0;
    double real_v_pos = 0;
    double real_h_vel = 20;
    double real_v_vel = 60;

    double measured_angle = 0;
    double measured_v = 0;

    //Variables to update H
    double h_pos = 0;
    double v_pos = 0;

    if(rc_kalman_new_alloc(&kf, F, G, H, Q, R, Pi, x_0)==-1) return -1;

    for(int i = 0; i <= 1000; i++){
        if(rc_kalman_predict_ekf(&kf, F, G, u) == -1) return -1;
        if(i % 100 == 0){
            y.d[0] = atan2(real_v_pos, real_h_pos) + ((rand()%1000)/1000 - 0.5)*sqrt(R.d[0][0]);
            y.d[1] = real_v_pos + ((rand()%1000)/1000 - 0.5)*sqrt(R.d[1][1]);
            if(rc_kalman_prediction_update_ekf(&kf, H, y) == -1) return -1;

            measured_angle = y.d[0];
            measured_v = y.d[1];
        }
        fprintf(fpt,"%f, %f, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n",
                                        time_delta*i ,kf.x_pre.d[0], kf.x_pre.d[1], kf.x_pre.d[2],
                                        kf.x_pre.d[3], y.d[0], kf.P.d[0][0], kf.P.d[1][1],
                                        kf.P.d[2][2], kf.P.d[3][3], real_h_pos, real_v_pos, 
                                        real_h_vel, real_v_vel, measured_angle, measured_v);
        // Update H
        h_pos = kf.x_pre.d[0];
        v_pos = kf.x_pre.d[1];
        H.d[0][0] = -v_pos/((h_pos*h_pos)+(v_pos*v_pos));
        H.d[1][0] = h_pos/((h_pos*h_pos)+(v_pos*v_pos));

        //Update Real Values
        real_h_pos = 0 + real_h_vel*time_delta*i;
        real_v_pos = 0 + 50*time_delta*i - 9.8*0.5*time_delta*time_delta*i*i ;
        real_h_vel = 30;
        real_v_vel = 50 - 9.8*time_delta*i;

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
    rc_matrix_free(&Pi);
    rc_vector_free(&u);
    rc_vector_free(&y);
    rc_kalman_free(&kf);
    return 0;
}
