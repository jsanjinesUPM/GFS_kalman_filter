import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import random


file =  open('wave_data.csv', 'w', newline='')
writer = csv.writer(file)

velocity_file =  open('velocity_data.csv', 'w', newline='')
writer = csv.writer(velocity_file)
#These are X,Y,Z coordinates, will later need to be transformed according to boat velocity and direction
velocity_file.write('vxSensor1,vySensor1,vzSensor1,vxSensor2,vySensor2,vzSensor2,'
                    'vxSensor3,vySensor3,vzSensor3,vxSensor4,vySensor4,vzSensor4,\n')

amplitude_file =  open('amplitud_data.csv', 'w', newline='')
amplitude_writer = csv.writer(amplitude_file)
amplitude_file.write('a0,a1,a2,a3,a4,a5,a6,a7,a10,a11,a12,a13,a14,a15,a16,a17\n')

sensors_file = open('sensors.csv', 'w', newline='')


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

measurement_period = 0.025
tRenovation = 10
time = np.arange(0,20,0.025)
g = 9.8
h = 20
z = -0.6
####Generating Random K####
angle = random.randint(0,900)*math.pi/(10*180)
module = random.randint(400,2500)*2/1000
k = np.array([[0.7969557584733965, 0.06972459419812653]])
#k = np.array([[module*math.cos(angle),module*math.sin(angle)]])
print("module: "+ str(module))
print("angle: "+ str(angle*180/math.pi))
###########################
x_y_0 = np.array([[1,1],[2,2],[2,1],[1,2],[1,3],[2,3],[3,3],[3,1],[3,2]])
x_y = np.zeros(np.shape(x_y_0))
x_boat_speed = 0 #in m/s X SPEED IS NOW TIME DEPENDANT AND IS DEFINED IN THE MAIN FOR LOOP
y_boat_speed = 0 #in m/s

###### Amplitudes ######
a = np.array([0,0,0,2,0,2,0,0,], dtype=float)
########################

num_rows = np.shape(time)[0]
num_cols = np.shape(x_y)[0]
num_k = np.shape(k)[0]
y = np.zeros((num_rows,num_cols))
vx = np.zeros((num_rows,num_cols))
vy = np.zeros((num_rows,num_cols))
vz = np.zeros((num_rows,num_cols))

##### Noise Matrix ########
P = np.identity(num_k*8)
Q = P*measurement_period/tRenovation

for i in range(0, num_rows):
    row = ''
    sensors_row = ''
    velocity_row = ''
    for j in range(0, num_cols):
        M = np.random.normal(size=np.shape(Q))
        M = M*Q
        for k_pos in range(0, num_k):
            k_mod = math.sqrt(k[k_pos][0]**2+k[k_pos][1]**2)
            w = math.sqrt(g*k_mod*math.tanh(k_mod*h))
            for m in range(0,np.shape(Q)[0]):
                a[m] = a[m] + M[m][m]
            x_y[j][0] = x_y_0[j][0] + x_boat_speed * time[i]
            x_y[j][1] = x_y_0[j][1] + y_boat_speed * time[i]
            y[i][j] += get_eta(time[i],x_y[j][0],x_y[j][1],k[k_pos][0],k[k_pos][1],w, a[k_pos*8:k_pos*8+8])
            vx[i][j] += get_vx(time[i],x_y[j][0],x_y[j][1],k[k_pos][0],k[k_pos][1],k_mod, w, a[k_pos*8:k_pos*8+8], g, h, z)
            vy[i][j] += get_vy(time[i],x_y[j][0],x_y[j][1],k[k_pos][0],k[k_pos][1],k_mod, w, a[k_pos*8:k_pos*8+8], g, h, z)
            vz[i][j] += get_vz(time[i],x_y[j][0],x_y[j][1],k[k_pos][0],k[k_pos][1],k_mod, w, a[k_pos*8:k_pos*8+8], g, h, z)
            row = row + str(y[i][j])+','
            sensors_row = sensors_row + str(x_y[j][0])+',' + str(x_y[j][1])+','
            velocity_row = velocity_row + str(vx[i][j])+','+ str(vy[i][j]) + ',' + str(vz[i][j]) + ','
    file.write(row+sensors_row+"\n")
    sensors_file.write(sensors_row+"\n")
    velocity_file.write(velocity_row+"\n")
    for k_pos in range(0, num_k):
        amplitude_file.write(str(a[k_pos*8+0])+','+str(a[k_pos*8+1])+','+str(a[k_pos*8+2])+','+str(a[k_pos*8+3])+','
                            +str(a[k_pos*8+4])+','+str(a[k_pos*8+5])+','+str(a[k_pos*8+6])+','+str(a[k_pos*8+7])+',')
    amplitude_file.write("\n")
file.close()
sensors_file.close()
amplitude_file.close()
velocity_file.close()
