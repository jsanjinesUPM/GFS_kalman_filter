import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv('c:/Users/mateo/GFS_kalman_filter/build/random_gaussian.csv', skip_blank_lines=True, header=0, index_col=False)


fig, ax = plt.subplots(2,1, sharex=True, sharey=True)
ax[1].scatter(df['random_scaled_x'], df['random_scaled_y'], label="Random Scaled")
ax[1].legend()
ax[0].scatter(df['random_normal_x'], df['random_normal_y'], label="Random Normal")
ax[0].legend()
plt.show()