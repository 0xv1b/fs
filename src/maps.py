import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as  plt
import h5py
import glob
from scipy.optimize import curve_fit, leastsq
from scipy.interpolate import CubicSpline
import re
#from obspy.signal.detrend import polynomial
from .helpers import common, detrend, constants
from pprint import pp

from matplotlib import rcParams

rcParams['axes.linewidth'] = 1.5

xscan180_keys = [1.14,1.142,1.144,1.146,1.148,1.15,1.152,1.155,1.16,1.1649999999999998,1.1699999999999997,1.17,1.1749999999999998]
xscan_keys = [1.148,1.1529999999999998,1.1579999999999997,1.1629999999999996,1.1679999999999995,1.1729999999999994,1.1779999999999993]

phase_map = 113
current_map = common.func_current(constants.current_magnitude, constants.current_phase, phase_map)


x_grid, y_grid = np.meshgrid(np.linspace(1.13, 1.19, 61), np.linspace(13.35, 13.395, 36))
#fig, axes = plt.subplots(2, 1, figsize=(8, 8))
#ax.text(1.129, 13.395, 'a)', fontsize=12, color='white',bbox=dict(facecolor='darkgray', edgecolor='none'))
plt.figure(1, figsize=(6,6))
plt.xlabel(r'x-axis ($\mu$m)',fontsize=18)
plt.ylabel(r'y-axis ($\mu$m)',fontsize=18)
#for ax in axes:
 #   ax.tick_params(axis='x', direction='in', labelsize=12)
  #  ax.tick_params(axis='y', direction='in', labelsize=12)
 #   ax.set_xlabel('x position (mm)', fontsize=16)
  #  ax.set_ylabel('y position (mm)', fontsize=16)
  #  ax.set_xlim(x_grid.min(), x_grid.max())
 #   ax.set_ylim(y_grid.min(), y_grid.max())
  #  ax.set_aspect('equal', adjustable='box')

normed_current_map = current_map.T - (( np.max(current_map) + np.min(current_map))/2)
photo_plot = plt.pcolormesh(x_grid, y_grid, (current_map.T - (( np.max(current_map) + np.min(current_map))/2))/np.max(normed_current_map), cmap='seismic', shading='auto', rasterized=True)
#photo_plot = axes[0].pcolormesh(x_grid, y_grid, current_map.T, cmap='seismic', shading='auto')

x = 1.14612
y = 13.36196

plt.plot(x,y, marker='x', markersize=16,markeredgewidth=4, color='black')


# for key in xscan180_keys:
#   fig = plt.gcf()
#   ax = fig.gca()
#   circle = plt.Circle((key, 13.365), 0.0025, color='red', fill=False)
#   ax.add_patch(circle)
#   plt.plot(key,13.365, marker='x', markersize=12,markeredgewidth=2, color='red')

  
# for key in xscan_keys:
#   fig = plt.gcf()
#   ax = fig.gca()
#   plt.plot(key,13.361, marker='x', markersize=12,markeredgewidth=2, color='yellow')

pcbar = plt.colorbar(photo_plot, pad=0.05, location="top", ticks=[-1, 0, 1])
pcbar.ax.tick_params(labelsize=14)
pcbar.set_label(label=r'$j_{ph}$', fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=14)
#ax = plt.gca()
#ax.text(1.13, 13.39, 'b)', fontsize=12, color='white', weight='bold')
plt.savefig("./images/Photocurrent_Map.pdf")


plt.figure(2, figsize=(6,6))
plt.xlabel(r'x-axis ($\mu$m)',fontsize=18)
plt.ylabel(r'y-axis ($\mu$m)',fontsize=18)





normed_refl_map = constants.refl_map.T - np.min(constants.refl_map.T)

refl_plot = plt.pcolormesh(x_grid, y_grid, ((constants.refl_map.T- np.min(constants.refl_map.T))/np.max(normed_refl_map)) , cmap='gray', shading='auto', rasterized=True)

cbar = plt.colorbar(refl_plot, pad=0.05, location="top",ticks=[0, 0.5, 1])
cbar.ax.tick_params(labelsize=14)
cbar.set_label(label='Reflection (a.u.)', fontsize=18)

plt.tick_params(axis='both', which='major', labelsize=14)

# for key in xscan180_keys:
#XSCAN 180
plt.plot(1.17,13.37, marker='x', markersize=16,markeredgewidth=4, color='blue')
plt.plot(1.146,13.37, marker='x',  markersize=16,markeredgewidth=4,color='red')
  
# for key in xscan_keys:
#XSCAN
# plt.plot(1.1729999999999994,13.36, marker='x', markersize=16,markeredgewidth=4, color='blue')
# plt.plot(1.148,13.36, marker='x',  markersize=16,markeredgewidth=4,color='red')


plt.tight_layout()
plt.savefig("./images/Reflectance_Map.pdf")
plt.show()

# zero needs to be at the 0 position in the map
# print("Min: " + str(np.min(current_map)))
# print("Max: " + str(np.max(current_map)))
# print("Middle: " + str(( np.max(current_map) + np.min(current_map))/2))