import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv('c:/Users/mateo/GFS_kalman_filter/build/random_gaussian.csv', skip_blank_lines=True, header=0, index_col=False)


fig, ax = plt.subplots(2,1)
ax[0].hist(df['random_normal'], bins=100, label="Random Normal")
ax[0].legend()
ax[1].hist(df['random_scaled'], bins=100, label="Random Scaled")
ax[1].legend()
plt.show()