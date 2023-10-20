import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv('../build/kalman_results.csv', skip_blank_lines=True, header=0)

fig, ax = plt.subplots(2,2)
#Plot for X position, Plotting Estimation, Real Value and Stimation Deviation
ax[0,0].plot(df['time'], df['h_position'], color='green', label="Estimated X")
ax[0,0].plot(df['time'], df['real_h_pos'], color='red', label="Real X")
ax[0,0].plot(df['time'], df['h_position'] + df['pi_h_pos'], color='black', linestyle="--", linewidth="0.5", label="X Estimate Deviation")
ax[0,0].plot(df['time'], df['h_position'] - df['pi_h_pos'], color='black', linestyle="--", linewidth="0.5")
ax[0,0].legend()
#Plot for Y position, Plotting Estimation, Real Value, Stimation Deviation and Measurements
ax[1,0].plot(df['time'], df['v_position'], color='green', label="Estimated Y")
ax[1,0].plot(df['time'], df['real_v_pos'], color='red', label="Real Y")
ax[1,0].plot(df['time'], df['v_position'] + df['pi_v_pos'], color='black', linestyle="--", linewidth="0.5", label="X Estimate Deviation")
ax[1,0].plot(df['time'], df['v_position'] - df['pi_v_pos'], color='black', linestyle="--", linewidth="0.5")
df['measured_v'] = df['measured_v'].replace({0:np.nan})
ax[1,0].scatter(df['time'], df['measured_v'], color='black', label="Measurements")
ax[1,0].legend()
#Plot for X Velocity position, Plotting Estimation, Real Value and Stimation Deviation
ax[0,1].plot(df['time'], df['h_velocity'], color='green', label="Estimated X Velocity")
ax[0,1].plot(df['time'], df['real_h_vel'], color='red', label="Real X Velocity")
ax[0,1].plot(df['time'], df['h_velocity'] + df['pi_h_vel'], color='black', linestyle="--", linewidth="0.5", label="X Estimate Deviation")
ax[0,1].plot(df['time'], df['h_velocity'] - df['pi_h_vel'], color='black', linestyle="--", linewidth="0.5")
ax[0,1].legend()
#Plot for Y Velocity position, Plotting Estimation, Real Value and Stimation Deviation
ax[1,1].plot(df['time'], df['v_velocity'], color='green', label="Estimated Y Velocity")
ax[1,1].plot(df['time'], df['real_v_vel'], color='red', label="Real Y Velocity")
ax[1,1].plot(df['time'], df['v_velocity'] + df['pi_v_vel'], color='black', linestyle="--", linewidth="0.5", label="X Estimate Deviation")
ax[1,1].plot(df['time'], df['v_velocity'] - df['pi_v_vel'], color='black', linestyle="--", linewidth="0.5")
ax[1,1].legend()
plt.show()