import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import random

g = 9.8
h = 20
z = -0.6
k = np.array([[0.7969557584733965, 0.06972459419812653]]) #Careful when using random K this will be overwritten

if __name__ == '__main__':
    file =  open('wave_data.csv', 'w', newline='')
    writer = csv.writer(file)

    amplitude_file =  open('amplitud_data.csv', 'w', newline='')
    amplitude_writer = csv.writer(amplitude_file)
    amplitude_file.write('a0,a1,a2,a3,a4,a5,a6,a7,a10,a11,a12,a13,a14,a15,a16,a17\n')

    sensors_file = open('sensors.csv', 'w', newline='')
    points_file = open('measuring_point.csv', 'w', newline='')


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

    ####Generating Random K####
    angle = random.randint(0,900)*math.pi/(10*180)
    module = random.randint(400,2500)*2/1000
    #k = np.array([[module*math.cos(angle),module*math.sin(angle)]])
    print("module: "+ str(module))
    print("angle: "+ str(angle*180/math.pi))
    ###########################
    x_y_0 = np.array([[1,1],[2,2],[2,1],[1,2],])#[1,3],[2,3],[3,3],[3,2],[3,1]])
    x_y = np.zeros(np.shape(x_y_0), float)
    measure_point = np.array([[0.5,0.5]], float)
    x_boat_speed = 8 #in m/s X iS NOW TIME DEPENDANT
    y_boat_speed = 8 #in m/s

    ###### Amplitudes ######
    a = np.array([0,0,0,2,0,2,0,0,], dtype=float)
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
        sensors_row = ''
        points_row = ''
        for j in range(0, num_cols):
            M = np.random.normal(size=np.shape(Q))
            M = M*Q
            for k_pos in range(0, num_k):
                kx = k[k_pos][0]
                ky = k[k_pos][1]
                k_mod = k_mod = math.sqrt(kx**2+ky**2)
                w = math.sqrt(g*k_mod*math.tanh(k_mod*h))
                for m in range(0,np.shape(Q)[0]):
                    a[m] = a[m] + M[m][m]
                x_y[j][0] = x_y_0[j][0] + x_boat_speed * time[i]
                x_y[j][1] = x_y_0[j][1] + y_boat_speed * time[i]
                if(j == 0):
                    measure_point[j][0] = measure_point[j][0] + x_boat_speed * time[i] #+ math.sin(2*math.pi*0.3*time[i])*math.sin(math.atan(y_boat_speed/x_boat_speed))
                    measure_point[j][1] = measure_point[j][1] + y_boat_speed * time[i] #+ math.sin(2*math.pi*0.3*time[i])*math.cos(math.atan(y_boat_speed/x_boat_speed))
                    points_row = points_row + str(measure_point[j][0])+',' + str(measure_point[j][1])+','
                y[i][j] += get_eta(time[i],x_y[j][0],x_y[j][1],kx,ky,w, a[k_pos*8:k_pos*8+8])
                row = row + str(y[i][j])+','
                sensors_row = sensors_row + str(x_y[j][0])+',' + str(x_y[j][1])+','
        file.write(row+sensors_row+"\n")
        sensors_file.write(sensors_row+"\n")
        points_file.write(points_row+"\n")
        for k_pos in range(0, num_k):
            amplitude_file.write(str(a[k_pos*8+0])+','+str(a[k_pos*8+1])+','+str(a[k_pos*8+2])+','+str(a[k_pos*8+3])+','
                                +str(a[k_pos*8+4])+','+str(a[k_pos*8+5])+','+str(a[k_pos*8+6])+','+str(a[k_pos*8+7])+',')
        amplitude_file.write("\n")
    file.close()
    sensors_file.close()
    amplitude_file.close()
    points_file.close()
