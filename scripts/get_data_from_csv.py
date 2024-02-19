import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

real_distance = 0.7855
path = 'C:/Users/mateo/GFS_kalman_filter/real_measured_data/OLAS_0_2_PERIODO_2_Y_3-5.csv'
measurements = pd.read_csv(path, skip_blank_lines=True)
measurements['True_Distance'] = (measurements['Distance']*100 + measurements['Distance1'])/10000 - real_distance
measurements['True_Distance'].to_csv('amp_0-2_t_2_3-5.csv', index=False)
