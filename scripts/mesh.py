import numpy as np
import matplotlib.pyplot as plt
import math
import csv

start_k = 0.8
stop_k = 5
number_modules = 1
start_angle = 5
stop_angle = 80
number_angles = 4

file =  open('mesh.csv', 'w', newline='')
writer = csv.writer(file)


modules = np.linspace(start_k, stop_k, number_modules)
angles = np.linspace(start_angle, stop_angle, number_angles)
angles = angles * math.pi/180
i = 0
for m in modules:
    for a in angles:
        file.write(str(m*math.cos(a))+','+str(m*math.sin(a))+"\n")
        i += 1
file.close()