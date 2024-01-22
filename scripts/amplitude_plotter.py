import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df_kalman = pd.read_csv('c:/Users/mateo/GFS_kalman_filter/build/wave_results.csv', skip_blank_lines=True, header=0)
df_amplitudes = pd.read_csv('c:/Users/mateo/GFS_kalman_filter/amplitud_data.csv', skip_blank_lines=True, header=0)

fig, ax = plt.subplots(2,4)
#Plot for X position, Plotting Estimation, Real Value and Stimation Deviation
for i in range(0,2):
    for j in range(0,4):
        ax[i,j].plot(df_kalman['time'], df_kalman.iloc[:,i*4+j+1], color='green')
        ax[i,j].plot(df_kalman['time'], df_amplitudes.iloc[:,i*4+j], color='red')
        ax[i,j].set_xlabel('t[s]')
        ax[i,j].set_ylabel('A[m]')
        ax[i,j].set_title('Amplitud '+str(i*4+j)+' frente al tiempo')
fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.5, hspace=0.4)
plt.show()
