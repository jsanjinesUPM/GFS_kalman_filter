import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import time

def get_vx(t,x,y,k_1,k_2, k_mod, w, a, g, h, z):
    S = k_1*g*math.cosh(k_mod*(h+z))/(w*math.cosh(k_mod*h))
    v  = -a[0]*math.cos(k_1*x)*math.sin(k_2*y)*math.cos(w*t)
    v += a[1]*math.cos(k_1*x)*math.sin(k_2*y)*math.sin(w*t)
    v += -a[2]*math.cos(k_1*x)*math.cos(k_2*y)*math.cos(w*t)
    v += a[3]*math.cos(k_1*x)*math.cos(k_2*y)*math.sin(w*t)
    v += a[4]*math.sin(k_1*x)*math.sin(k_2*y)*math.cos(w*t)
    v += -a[5]*math.sin(k_1*x)*math.sin(k_2*y)*math.sin(w*t)
    v += a[6]*math.sin(k_1*x)*math.cos(k_2*y)*math.cos(w*t)
    v += -a[7]*math.sin(k_1*x)*math.cos(k_2*y)*math.sin(w*t)
    return S*v

def get_vy(t,x,y,k_1,k_2, k_mod, w, a, g, h, z):
    S = k_2*g*math.cosh(k_mod*(h+z))/(w*math.cosh(k_mod*h))
    v  = -a[0]*math.sin(k_1*x)*math.cos(k_2*y)*math.cos(w*t)
    v += a[1]*math.sin(k_1*x)*math.cos(k_2*y)*math.sin(w*t)
    v += a[2]*math.sin(k_1*x)*math.sin(k_2*y)*math.cos(w*t)
    v += -a[3]*math.sin(k_1*x)*math.sin(k_2*y)*math.sin(w*t)
    v += -a[4]*math.cos(k_1*x)*math.cos(k_2*y)*math.cos(w*t)
    v += a[5]*math.cos(k_1*x)*math.cos(k_2*y)*math.sin(w*t)
    v += a[6]*math.cos(k_1*x)*math.sin(k_2*y)*math.cos(w*t)
    v += -a[7]*math.cos(k_1*x)*math.sin(k_2*y)*math.sin(w*t)
    return S*v

def get_vz(t,x,y,k_1,k_2, k_mod, w, a, g, h, z):
    S = k_mod*g*math.sinh(k_mod*(h+z))/(w*math.cosh(k_mod*h))
    v  = -a[0]*math.sin(k_1*x)*math.sin(k_2*y)*math.cos(w*t)
    v += a[1]*math.sin(k_1*x)*math.sin(k_2*y)*math.sin(w*t)
    v += -a[2]*math.sin(k_1*x)*math.cos(k_2*y)*math.cos(w*t)
    v += a[3]*math.sin(k_1*x)*math.cos(k_2*y)*math.sin(w*t)
    v += -a[4]*math.cos(k_1*x)*math.sin(k_2*y)*math.cos(w*t)
    v += a[5]*math.cos(k_1*x)*math.sin(k_2*y)*math.sin(w*t)
    v += -a[6]*math.cos(k_1*x)*math.cos(k_2*y)*math.cos(w*t)
    v += a[7]*math.cos(k_1*x)*math.cos(k_2*y)*math.sin(w*t)
    return S*v


df_real_vel = pd.read_csv('c:/Users/mateo/GFS_kalman_filter/velocity_data.csv', skip_blank_lines=True, header=0)
df_kalman_amplitudes = pd.read_csv('c:/Users/mateo/GFS_kalman_filter/build/wave_results.csv', skip_blank_lines=True, header=0)
sensors = pd.read_csv('c:/Users/mateo/GFS_kalman_filter/sensors.csv', skip_blank_lines=True, header=None)
times = df_kalman_amplitudes['time'].tolist()
a0 = df_kalman_amplitudes.columns.get_loc('a0')
vx_sensor1 = df_real_vel.columns.get_loc('vxSensor1')

mesh = pd.read_csv('c:/Users/mateo/GFS_kalman_filter/mesh.csv', skip_blank_lines=True)
k= mesh.to_numpy()

g = 9.8
h = 20
z = -0.6

num_sensors = int(len(sensors.columns)/2)
num_k = np.shape(k)[0]
num_rows = len(df_kalman_amplitudes)

vx = np.zeros((num_rows,num_sensors))
vy = np.zeros((num_rows,num_sensors))
vz = np.zeros((num_rows,num_sensors))

for i in range(0, num_rows):
    for j in range(0, num_sensors):
        for k_pos in range(0, num_k):
            kx = k[k_pos][0]
            ky = k[k_pos][1]
            k_mod = math.sqrt(kx**2+ky**2)
            w = math.sqrt(g*k_mod*math.tanh(k_mod*h))
            start_time = time.time()

            x = sensors.iloc[i,j*2+0]
            y = sensors.iloc[i,j*2+1]
            a = df_kalman_amplitudes.iloc[i,k_pos*8+a0:k_pos*8+8+a0].tolist()

            vx[i][j] += get_vx(times[i],x,y,kx,ky,k_mod, w, a, g, h, z)
            vy[i][j] += get_vy(times[i],x,y,kx,ky,k_mod, w, a, g, h, z)
            vz[i][j] += get_vz(times[i],x,y,kx,ky,k_mod, w, a, g, h, z)

fig, ax = plt.subplots(2,1)
ax[0].plot(times, vx[:,0], color='green', label='Kalman')
ax[0].plot(times, df_real_vel.iloc[:,vx_sensor1], color='red', label='Real')
ax[0].legend(loc='best')
ax[0].set_xlabel('t[s]')
ax[0].set_ylabel('v[m/s]')
ax[0].set_title('Kalman Wave Speed vs Real Wave Speed')


ax[1].plot(times, abs(df_real_vel.iloc[:,vx_sensor1]-vx[:,0]), color='Blue')
ax[1].set_xlabel('t[s]')
ax[1].set_ylabel('Error[m/s]')
ax[1].set_title('Speed Absolute Error')
fig.subplots_adjust(left=0.075, bottom=0.1, right=0.95, top=0.9, wspace=0.5, hspace=0.4)
plt.show()