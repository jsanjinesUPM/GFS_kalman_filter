#include "filter_functions.h"
#include "kalman.h"
#include "algebra.h"
#include "other.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAXCHAR 1024

void write_header(FILE *writefile, double time, int x_size, int y_size){
    //Writing header to CSV file
    fprintf(writefile, "time,");
    for(int i = 0; i < x_size-3; i++){
        fprintf(writefile, "a%d,", i);
    }
    fprintf(writefile, "z,");
    fprintf(writefile, "vz,");
    fprintf(writefile, "az,");
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
    if(H->cols != K.rows*8 + 3){
        return -1;
    }
    for(int i = 0; i < H->rows-1; i++){
        for(int j = 0; j < K.rows; j++){
            H->d[i][j*8+0] = -sin(K.d[j][0]*X_Y.d[i][0])*sin(K.d[j][1]*X_Y.d[i][1])*sin(w.d[j]*time);
            H->d[i][j*8+1] = -sin(K.d[j][0]*X_Y.d[i][0])*sin(K.d[j][1]*X_Y.d[i][1])*cos(w.d[j]*time);
            H->d[i][j*8+2] = -sin(K.d[j][0]*X_Y.d[i][0])*cos(K.d[j][1]*X_Y.d[i][1])*sin(w.d[j]*time);
            H->d[i][j*8+3] = -sin(K.d[j][0]*X_Y.d[i][0])*cos(K.d[j][1]*X_Y.d[i][1])*cos(w.d[j]*time);
            H->d[i][j*8+4] = -cos(K.d[j][0]*X_Y.d[i][0])*sin(K.d[j][1]*X_Y.d[i][1])*sin(w.d[j]*time);
            H->d[i][j*8+5] = -cos(K.d[j][0]*X_Y.d[i][0])*sin(K.d[j][1]*X_Y.d[i][1])*cos(w.d[j]*time);
            H->d[i][j*8+6] = -cos(K.d[j][0]*X_Y.d[i][0])*cos(K.d[j][1]*X_Y.d[i][1])*sin(w.d[j]*time);
            H->d[i][j*8+7] = -cos(K.d[j][0]*X_Y.d[i][0])*cos(K.d[j][1]*X_Y.d[i][1])*cos(w.d[j]*time);
        }
        H->d[i][(K.rows)*8 + 0] = 1; //Position Column
        H->d[i][(K.rows)*8 + 1] = 0; //Velocity Column
        H->d[i][(K.rows)*8 + 2] = 0; //Acceleration Column
    }
    //Constructing Acceleration Measurement Row 
    for(int i = 0; i < H->cols-3; i++){
        H->d[H->rows-1][i] = 0;
    }
    H->d[H->rows-1][H->cols-3] = 0;
    H->d[H->rows-1][H->cols-2] = 0;
    H->d[H->rows-1][H->cols-1] = 1;
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