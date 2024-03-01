import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df_kalman = pd.read_csv('c:/Users/mateo/GFS_kalman_filter/build/wave_results.csv', skip_blank_lines=True, header=0)
df_amplitudes = pd.read_csv('c:/Users/mateo/GFS_kalman_filter/amplitud_data.csv', skip_blank_lines=True, header=0)

manual_filter_k_num = 2
k_num = max(int(df_amplitudes.shape[1]/8), manual_filter_k_num) #number of columns in dataframe/8
for k in range(0,k_num):
    fig, ax = plt.subplots(2,4)
    for i in range(0,2):
        for j in range(0,4):
            try:
                ax[i,j].plot(df_kalman['time'], df_kalman.iloc[:,8*k+i*4+j+1], color='green')
            except:
                print("kalman values went wrong")
            try:
                ax[i,j].plot(df_kalman['time'], df_amplitudes.iloc[:,8*k+i*4+j], color='red')
            except:
                print("Wave data went wrong")
            ax[i,j].set_xlabel('t[s]')
            ax[i,j].set_ylabel('A[m]')
            ax[i,j].set_title('Amplitud '+str(i*4+j)+' frente al tiempo')
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.5, hspace=0.4)
    plt.show()
