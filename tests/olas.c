/*
* @brief we aim to test our previously implemented kalman filter
* to predict sea state
*
*/
#include "unity.h"
#include "kalman.h"
#include "algebra.h"
#include "other.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define MAXCHAR 1024
//#define eta_analyzer 1


void write_header(FILE *writefile, double time, int x_size, int y_size){
    //Writing header to CSV file
    fprintf(writefile, "time,");
    for(int i = 0; i < x_size; i++){
        fprintf(writefile, "a%d,", i);
    }
    for(int i = 0; i < y_size; i++){
        fprintf(writefile, "h%d,", i);
        fprintf(writefile, "y%d,", i);
    }
    fprintf(writefile,"\n");
    return;
}

void write_to_file(FILE *writefile, double time, rc_vector_t x, rc_vector_t h, rc_vector_t y){
    //Writing to file
    fprintf(writefile, "%f,", time);
    for(int i = 0; i < x.len; i++){
        fprintf(writefile, "%lf,", x.d[i]);
    }
    for(int i = 0; i < h.len; i++){
        fprintf(writefile, "%lf,", h.d[i]);
        fprintf(writefile, "%lf,", y.d[i]);
    }
    fprintf(writefile,"\n");
    return;
}

int mesh_size(const char* meshfile_path){
    //Data will be read from a file
    FILE *meshfile;
    // readfile = fopen("../amp_0-2_t_2_3-5.csv", "r+");
    meshfile = fopen(meshfile_path, "r+");
    char row[MAXCHAR];
    int lines = 0;
    while(1){
        if(fgets(row, MAXCHAR, meshfile) == NULL) break;
        lines++;
    }
    fclose(meshfile);
    return lines;
}

void load_mesh(rc_matrix_t* K, const char* meshfile_path){
    //Data will be read from a file
    FILE *meshfile;
    // readfile = fopen("../amp_0-2_t_2_3-5.csv", "r+");
    meshfile = fopen(meshfile_path, "r+");
    char row[MAXCHAR];
    char *token;
    int column = 0; //Will be used to iterate though CSV columns
    int row_number = 0;
    while(1){
        if(fgets(row, MAXCHAR, meshfile) == NULL) break;
        token = strtok(row, ",");
        column = 0;
        while(token != NULL)
        {   
            if(column == K->cols) break;
            K->d[row_number][column] = strtod(token,NULL);
            token = strtok(NULL, ",");
            column++;
        }
        row_number++;
    }
    fclose(meshfile);
}

int get_w(rc_vector_t* w, rc_matrix_t K, double depth){
    if(w->len != K.rows){
        return -1;
    }
    double k_norm = 0;
    double g = 9.8;
    for(int i = 0; i < w->len; i++){
        k_norm = sqrt(K.d[i][0]*K.d[i][0] + K.d[i][1]*K.d[i][1]);
        w->d[i] = sqrt(g*k_norm*tanh(k_norm*depth));
    }
    return 0;
}

int get_H(rc_matrix_t* H, rc_matrix_t K, rc_matrix_t X_Y, rc_vector_t w, double time){
    if(H->cols != K.rows*8){
        return -1;
    }
    for(int i = 0; i < H->rows; i++){
        for(int j = 0; j < K.rows; j++){
            H->d[i][j*8+0] = sin(K.d[j][0]*X_Y.d[i][0])*sin(K.d[j][1]*X_Y.d[i][1])*sin(w.d[j]*time);
            H->d[i][j*8+1] = sin(K.d[j][0]*X_Y.d[i][0])*sin(K.d[j][1]*X_Y.d[i][1])*cos(w.d[j]*time);
            H->d[i][j*8+2] = sin(K.d[j][0]*X_Y.d[i][0])*cos(K.d[j][1]*X_Y.d[i][1])*sin(w.d[j]*time);
            H->d[i][j*8+3] = sin(K.d[j][0]*X_Y.d[i][0])*cos(K.d[j][1]*X_Y.d[i][1])*cos(w.d[j]*time);
            H->d[i][j*8+4] = cos(K.d[j][0]*X_Y.d[i][0])*sin(K.d[j][1]*X_Y.d[i][1])*sin(w.d[j]*time);
            H->d[i][j*8+5] = cos(K.d[j][0]*X_Y.d[i][0])*sin(K.d[j][1]*X_Y.d[i][1])*cos(w.d[j]*time);
            H->d[i][j*8+6] = cos(K.d[j][0]*X_Y.d[i][0])*cos(K.d[j][1]*X_Y.d[i][1])*sin(w.d[j]*time);
            H->d[i][j*8+7] = cos(K.d[j][0]*X_Y.d[i][0])*cos(K.d[j][1]*X_Y.d[i][1])*cos(w.d[j]*time);
        }
    }
    return 0;
}

int add_H_prime(rc_matrix_t* N, rc_matrix_t H, rc_matrix_t K, rc_matrix_t X_Y, rc_vector_t w, double time){
    if(H.cols != K.rows*8){
        return -1;
    }
    for(int i = 0; i < H.rows; i++){
        for(int j = 0; j < K.rows; j++){
            N->d[i][j*8+0] = w.d[j]*sin(K.d[j][0]*X_Y.d[i][0])*sin(K.d[j][1]*X_Y.d[i][1])*cos(w.d[j]*time);
            N->d[i][j*8+1] = -w.d[j]*sin(K.d[j][0]*X_Y.d[i][0])*sin(K.d[j][1]*X_Y.d[i][1])*sin(w.d[j]*time);
            N->d[i][j*8+2] = w.d[j]*sin(K.d[j][0]*X_Y.d[i][0])*cos(K.d[j][1]*X_Y.d[i][1])*cos(w.d[j]*time);
            N->d[i][j*8+3] = -w.d[j]*sin(K.d[j][0]*X_Y.d[i][0])*cos(K.d[j][1]*X_Y.d[i][1])*sin(w.d[j]*time);
            N->d[i][j*8+4] = w.d[j]*cos(K.d[j][0]*X_Y.d[i][0])*sin(K.d[j][1]*X_Y.d[i][1])*cos(w.d[j]*time);
            N->d[i][j*8+5] = -w.d[j]*cos(K.d[j][0]*X_Y.d[i][0])*sin(K.d[j][1]*X_Y.d[i][1])*sin(w.d[j]*time);
            N->d[i][j*8+6] = w.d[j]*cos(K.d[j][0]*X_Y.d[i][0])*cos(K.d[j][1]*X_Y.d[i][1])*cos(w.d[j]*time);
            N->d[i][j*8+7] = -w.d[j]*cos(K.d[j][0]*X_Y.d[i][0])*cos(K.d[j][1]*X_Y.d[i][1])*sin(w.d[j]*time);
        }
    }
    return 0;
}

int main(void) {

    const char meshfile_path[] = "../mesh.csv";

    double depth = 20;
    double time = 0;

    int K_rows = mesh_size(meshfile_path);
    int Num_x_y = 9;
    int y_size = Num_x_y;
    int x_size = K_rows*8;

    //Result will be written to a file
    FILE *writefile;
    writefile = fopen("wave_results.csv", "w+");
    write_header(writefile, time, x_size, y_size);
    
    //Data will be read from a file
    FILE *readfile;
    // readfile = fopen("../amp_0-2_t_2_3-5.csv", "r+");
    readfile = fopen("../wave_data.csv", "r+");
    char row[MAXCHAR];
    char *token;
    int i = 0; //Will be used to iterate though CSV columns


    rc_kalman_t kfilter  = RC_KALMAN_INITIALIZER;

    rc_matrix_t K = RC_MATRIX_INITIALIZER;
    rc_matrix_t X_Y = RC_MATRIX_INITIALIZER;
    rc_matrix_t H = RC_MATRIX_INITIALIZER;
    rc_matrix_t F = RC_MATRIX_INITIALIZER;
    rc_matrix_t G = RC_MATRIX_INITIALIZER;
    rc_matrix_t R = RC_MATRIX_INITIALIZER;
    rc_matrix_t Q = RC_MATRIX_INITIALIZER;
    rc_matrix_t Pi = RC_MATRIX_INITIALIZER;

    rc_vector_t w = RC_VECTOR_INITIALIZER;
    rc_vector_t X = RC_VECTOR_INITIALIZER;
    rc_vector_t y = RC_VECTOR_INITIALIZER;
    rc_vector_t u = RC_VECTOR_INITIALIZER;
    rc_vector_t U = RC_VECTOR_INITIALIZER;
    rc_vector_t x_0 = RC_VECTOR_INITIALIZER;
    rc_vector_t h = RC_VECTOR_INITIALIZER;

    rc_matrix_zeros(&K, K_rows, 2);
    rc_matrix_zeros(&X_Y, y_size, 2);
    rc_matrix_zeros(&H, y_size, x_size);
    rc_matrix_zeros(&F, x_size, x_size);
    rc_matrix_zeros(&Q, x_size, x_size);
    rc_matrix_zeros(&Pi, x_size, x_size);
    rc_matrix_zeros(&R, y_size, y_size);
    rc_matrix_zeros(&G, x_size, x_size);

    rc_vector_zeros(&X, x_size);
    rc_vector_zeros(&w, K_rows);
    rc_vector_zeros(&y, y_size);
    rc_vector_zeros(&h, y_size);
    rc_vector_zeros(&U, x_size);
    rc_vector_zeros(&u, x_size);
    rc_vector_zeros(&x_0, x_size);

    rc_matrix_identity(&F, F.cols);
    rc_matrix_identity(&Pi, Pi.cols);
    

    rc_matrix_identity(&R, R.cols);
    rc_matrix_times_scalar(&R, pow(0.02, 2));

    rc_matrix_identity(&Q, Q.cols);
    rc_matrix_times_scalar(&Q, 0.025/10);

    load_mesh(&K, meshfile_path);

    X_Y.d[0][0] = 1;
    X_Y.d[0][1] = 1;

    X_Y.d[1][0] = 2;
    X_Y.d[1][1] = 2;

    X_Y.d[2][0] = 2;
    X_Y.d[2][1] = 1;

    X_Y.d[3][0] = 1;
    X_Y.d[3][1] = 2;

    // X_Y.d[4][0] = 1;
    // X_Y.d[4][1] = 3;

    // X_Y.d[5][0] = 2;
    // X_Y.d[5][1] = 3;

    // X_Y.d[6][0] = 3;
    // X_Y.d[6][1] = 3;

    // X_Y.d[7][0] = 3;
    // X_Y.d[7][1] = 2;

    // X_Y.d[8][0] = 3;
    // X_Y.d[8][1] = 1;


    get_w(&w,K, depth);
    get_H(&H, K, X_Y, w, time);


    if(rc_kalman_new_alloc(&kfilter, F, G, H, Q, R, Pi, u, x_0)==-1) return -1;

    while (1){
        clock_t begin_predict = clock();
        if(rc_kalman_predict_simple(&kfilter, F) == -1) return -1;
        #ifdef eta_analyzer
            get_H(&H, K, X_Y, w, time);
            rc_matrix_times_col_vec(H, kfilter.x_pre, &h);
        #endif
        clock_t end_predict = clock();
        #ifdef eta_analyzer
        if((int)(time*1000)%25 == 0){
        #endif
            if(fgets(row, MAXCHAR, readfile) == NULL) break;
            token = strtok(row, ",");
            i = 0;
            while(token != NULL)
            {   
                if(i == y.len) break;
                y.d[i] = strtod(token,NULL);
                y.d[i] += rc_random_normal()*sqrt(R.d[i][i]);
                token = strtok(NULL, ",");
                i++;
            }
            i = 0;
            while(token != NULL)
            {   
                if(i == y.len) break;
                X_Y.d[i][0] = strtod(token,NULL);
                token = strtok(NULL, ",");
                X_Y.d[i][1] = strtod(token,NULL);
                token = strtok(NULL, ",");
                i++;
            }
            get_H(&H, K, X_Y, w, time);
            rc_matrix_times_col_vec(H, kfilter.x_pre, &h);
            clock_t begin_update = clock();
            if(rc_kalman_prediction_update_ekf(&kfilter, H, y, h) == -1) return -1;
            clock_t end_update = clock();
            double time_predict = (double)(end_predict - begin_predict) / CLOCKS_PER_SEC;
            double time_update = (double)(end_update - begin_update)/ CLOCKS_PER_SEC;
        // printf("Time per iteration predict: %lf \n", time_predict);
        // printf("Time per iteration update: %lf \n \n", time_update);
        #ifdef eta_analyzer
        }
        #endif
        //Writing to file
        write_to_file(writefile, time, kfilter.x_pre, h, y);
        #ifdef eta_analyzer
        time += 0.001;
        #else
        time += 0.025;
        #endif
    }
    #ifdef eta_analyzer
    for(int extra = 0; extra < 10000; extra ++){
        if(rc_kalman_predict_simple(&kfilter, F) == -1) return -1;
        get_H(&H, K, X_Y, w, time);
        rc_matrix_times_col_vec(H, kfilter.x_pre, &h);
        write_to_file(writefile, time, kfilter.x_pre, h, y);
        time += 0.001;
    }
    #endif

    fclose(writefile);
    fclose(readfile);

    rc_matrix_free(&K);
    rc_matrix_free(&X_Y);
    rc_matrix_free(&H);
    rc_matrix_free(&F);
    rc_matrix_free(&Q);
    rc_matrix_free(&Pi);
    rc_matrix_free(&R);
    rc_matrix_free(&G);

    rc_vector_zeros(&X, x_size);
    rc_vector_zeros(&w, K_rows);
    rc_vector_zeros(&y, y_size);
    rc_vector_zeros(&h, y_size);
    rc_vector_zeros(&U, x_size);
    rc_vector_zeros(&u, x_size);
    rc_vector_zeros(&x_0, x_size);

    rc_vector_free(&X);
    rc_vector_free(&w);
    rc_vector_free(&y);
    rc_vector_free(&h);
    rc_vector_free(&U);
    rc_vector_free(&u);
    rc_vector_free(&x_0);
    return 0;
}
