import numpy as np
import matplotlib.pyplot as plt
import math
import csv


file =  open('wave_data.csv', 'w', newline='')
writer = csv.writer(file)

amplitude_file =  open('amplitud_data.csv', 'w', newline='')
amplitude_writer = csv.writer(amplitude_file)
amplitude_file.write('a0,a1,a2,a3,a4,a5,a6,a7,a10,a11,a12,a13,a14,a15,a16,a17\n')


def get_eta(t,x,y,k_1,k_2,w, a):
    n  = a[0]*math.sin(k_1*x)*math.sin(k_2*y)*math.sin(w*t)
    n += a[1]*math.sin(k_1*x)*math.sin(k_2*y)*math.cos(w*t)
    n += a[2]*math.sin(k_1*x)*math.cos(k_2*y)*math.sin(w*t)
    n += a[3]*math.sin(k_1*x)*math.cos(k_2*y)*math.cos(w*t)
    n += a[4]*math.cos(k_1*x)*math.sin(k_2*y)*math.sin(w*t)
    n += a[5]*math.cos(k_1*x)*math.sin(k_2*y)*math.cos(w*t)
    n += a[6]*math.cos(k_1*x)*math.cos(k_2*y)*math.sin(w*t)
    n += a[7]*math.cos(k_1*x)*math.cos(k_2*y)*math.cos(w*t)
    return n

measurement_period = 0.025
tRenovation = 10
time = np.arange(0,10,0.025)
g = 9.8
h = 20
k = np.array([[0.707,0.707]])
x_y = np.array([[1,1],[2,2],[2,1],[1,2]])

###### Coeficients ######
a = np.array([0,0,0,0.5,0,0.5,0,0,  0,0,0,0,0,0,0,0,], dtype=float)
#a = np.array([[0,0,0,0,0,0,0,0],], dtype=float)
########################

num_rows = np.shape(time)[0]
num_cols = np.shape(x_y)[0]
num_k = np.shape(k)[0]
y = np.zeros((num_rows,num_cols))

##### Noise Matrix ########
P = np.identity(num_k*8)
Q = P*measurement_period/tRenovation

for i in range(0, num_rows):
    row = ''
    for j in range(0, num_cols):
        M = np.random.normal(size=np.shape(Q))
        M = M*Q
        for k_pos in range(0, num_k):
            k_mod = math.sqrt(k[k_pos][0]**2+k[k_pos][1]**2)
            w = math.sqrt(g*k_mod*math.tanh(k_mod*h))
            for m in range(0,np.shape(Q)[0]):
                a[m] = a[m] + M[m][m]
            y[i][j] += get_eta(time[i],x_y[j][0],x_y[j][1],k[k_pos][0],k[k_pos][1],w, a[k_pos*8:k_pos*8+8])
            row = row + str(y[i][j])+','
    file.write(row+"\n")
    for k_pos in range(0, num_k):
        amplitude_file.write(str(a[k_pos*8+0])+','+str(a[k_pos*8+1])+','+str(a[k_pos*8+2])+','+str(a[k_pos*8+3])+','
                            +str(a[k_pos*8+4])+','+str(a[k_pos*8+5])+','+str(a[k_pos*8+6])+','+str(a[k_pos*8+7])+',')
    amplitude_file.write("\n")
file.close()
amplitude_file.close()
