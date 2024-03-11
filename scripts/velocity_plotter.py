import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import time
from waves import g,h,z
from waves import k as k_sim

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


df_kalman_amplitudes = pd.read_csv('c:/Users/mateo/GFS_kalman_filter/build/wave_results.csv', skip_blank_lines=True, header=0)
df_sim_amplitudes = pd.read_csv('c:/Users/mateo/GFS_kalman_filter/amplitud_data.csv', skip_blank_lines=True, header=0)

points = pd.read_csv('c:/Users/mateo/GFS_kalman_filter/measuring_point.csv', skip_blank_lines=True, header=None)
times = df_kalman_amplitudes['time'].tolist()
a0 = df_kalman_amplitudes.columns.get_loc('a0')
sim_a0 = df_sim_amplitudes.columns.get_loc('a0')

mesh = pd.read_csv('c:/Users/mateo/GFS_kalman_filter/mesh.csv', skip_blank_lines=True, header=None)
mesh = pd.read_csv('c:/Users/mateo/GFS_kalman_filter/mesh.csv', skip_blank_lines=True, header=None)
k_mesh= mesh.to_numpy()

num_points = int(len(points.columns)/2)
num_k = np.shape(k_mesh)[0]
num_sim_k = np.shape(k_sim)[0]
num_rows = len(df_kalman_amplitudes)

vx = np.zeros((num_rows,num_points))
vy = np.zeros((num_rows,num_points))
vz = np.zeros((num_rows,num_points))
vx_sim = np.zeros((num_rows,num_points))
vy_sim = np.zeros((num_rows,num_points))
vz_sim = np.zeros((num_rows,num_points))


for i in range(0, num_rows):
    for j in range(0, num_points):
        for k_pos in range(0, num_k):
            kx = k_mesh[k_pos][0]
            ky = k_mesh[k_pos][1]
            k_mod = math.sqrt(kx**2+ky**2)
            w = math.sqrt(g*k_mod*math.tanh(k_mod*h))

            x = points.iloc[i,j*2+0]
            y = points.iloc[i,j*2+1]

            a = df_kalman_amplitudes.iloc[i,k_pos*8+a0:k_pos*8+8+a0].tolist()

            vx[i][j] += get_vx(times[i],x,y,kx,ky,k_mod, w, a, g, h, z)
            vy[i][j] += get_vy(times[i],x,y,kx,ky,k_mod, w, a, g, h, z)
            vz[i][j] += get_vz(times[i],x,y,kx,ky,k_mod, w, a, g, h, z)
        for k_pos in range(0, num_sim_k):
            kx = k_sim[k_pos][0]
            ky = k_sim[k_pos][1]
            k_mod = math.sqrt(kx**2+ky**2)
            w = math.sqrt(g*k_mod*math.tanh(k_mod*h))

            x = points.iloc[i,j*2+0]
            y = points.iloc[i,j*2+1]

            a = df_sim_amplitudes.iloc[i,k_pos*8+sim_a0:k_pos*8+8+sim_a0].tolist()

            vx_sim[i][j] += get_vx(times[i],x,y,kx,ky,k_mod, w, a, g, h, z)
            vy_sim[i][j] += get_vy(times[i],x,y,kx,ky,k_mod, w, a, g, h, z)
            vz_sim[i][j] += get_vz(times[i],x,y,kx,ky,k_mod, w, a, g, h, z)

#Lets calculate velocities parallel and perpendicular to the boats' heading
heading = math.atan((points.iloc[1,1]-points.iloc[0,1])/(points.iloc[1,0]-points.iloc[0,0]))
v_parallel = math.cos(heading)*vx + math.sin(heading)*vy
v_parallel_sim = math.cos(heading)*vx_sim + math.sin(heading)*vy_sim

v_perpendicular = -math.sin(heading)*vx + math.cos(heading)*vy
v_perpendicular = -math.sin(heading)*vx + math.cos(heading)*vy


fig, ax = plt.subplots(2,1)
# ax[0].plot(times, v_parallel[:,0], color='green', label='Kalman')
# ax[0].plot(times, v_parallel_sim[:,0], color='red', label='Real')
ax[0].plot(times, vz[:,0], color='green', label='Kalman')
ax[0].plot(times, vz_sim[:,0], color='red', label='Real')
ax[0].legend(loc='best')
ax[0].set_xlabel('t[s]')
ax[0].set_ylabel('v[m/s]')
ax[0].set_title('Kalman Wave Speed vs Real Wave Speed')


ax[1].plot(times, abs(v_parallel_sim[:,0]-v_parallel[:,0]), color='Blue')
ax[1].set_xlabel('t[s]')
ax[1].set_ylabel('Error[m/s]')
ax[1].set_title('Speed Absolute Error')
fig.subplots_adjust(left=0.075, bottom=0.1, right=0.95, top=0.9, wspace=0.5, hspace=0.4)
plt.show()