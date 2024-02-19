import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df_kalman = pd.read_csv('c:/Users/mateo/GFS_kalman_filter/build/wave_results.csv', skip_blank_lines=True, header=0)

h0 = df_kalman.columns.get_loc('h0')
y0 = df_kalman.columns.get_loc('y0')

# fig, ax = plt.subplots(2,1)
# ax[0].plot(df_kalman['time'], df_kalman.iloc[:,h0]-df_kalman.iloc[:,y0], color='green')
# ax[0].set_xlabel('t[s]')
# ax[0].set_ylabel('h[m]')
# ax[0].set_title('Kalman Wave Height - Real Wave Height')

plt.plot(df_kalman['time'], df_kalman.iloc[:,h0], color='green', label='Kalman')
plt.plot(df_kalman['time'], df_kalman.iloc[:,y0], color='red', label='Real')
plt.xlabel('t[s]')
plt.ylabel('h[m]')
plt.title('Kalman Wave Height vs Real Wave Height')
# ax[1].set_xlabel('t[s]')
# ax[1].set_ylabel('h[m]')
# ax[1].set_title('Kalman Wave Height vs Real Wave Height')
#ax[1].subplots_adjust(left=0.075, bottom=0.1, right=0.95, top=0.9, wspace=0.5, hspace=0.4)

plt.show()
