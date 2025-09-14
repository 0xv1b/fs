import numpy as np
import matplotlib.pyplot as  plt
import matplotlib as mpl
import h5py
import glob
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
import re
import os
from .helpers import common, detrend, constants
import warnings
warnings.filterwarnings("ignore")


mpl.rcParams['legend.title_fontsize'] = 'large'

xscan_keys = [1.148,1.1529999999999998,1.1579999999999997,1.1629999999999996,1.1679999999999995,1.1729999999999994,1.1779999999999993]

projection_phase = 113
xscans = {}
xscans_detrend = {}
currents_detrend = {}
wedge_pos_array = []
currents = {}

detrend.init_lists(constants.CEP_XSCAN, xscans, xscans_detrend, wedge_pos_array, projection_phase=projection_phase, xscan=True)
wedge_positions = detrend.get_wedge_positions(wedge_pos_array, decimals=2)
currents_detrend = detrend.get_detrended_currents(xscans, currents, deg=7)

d = 60 # This allows for basically exactly 5 peaks
pinit = [20,1,0.4,d,0.42,0]

currents_detrend_data = {}
wedge_pos_data = {}
start = [1,1,10,25,10,1,5 ]
end =   [1,1,42,16,30,1,10]
for i, key in enumerate(np.sort([*currents_detrend])):
    currents_detrend_data[key] = currents_detrend[key][start[i]:-end[i]]
    wedge_pos_data[key] = wedge_positions[start[i]:-end[i]]



fit_x = {}
fit_y = {}
wavelength = {}
IPL = {}
Theta = {}
for i, key in enumerate(np.sort([*currents_detrend])):
    fit_x[key] = np.linspace(common.trans2ins(wedge_pos_data[key][0]),common.trans2ins(wedge_pos_data[key][-1]),1000)
    popt,pcov = curve_fit(common.CEP_fit,common.trans2ins(wedge_pos_data[key]),currents_detrend_data[key]*1e-7*1e12,p0=pinit, maxfev=40000) # Should be 1e-8
    fit_y[key] = common.CEP_fit(fit_x[key],popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])
    # print(f"{key}:  {popt}")
    IPL[key] = 2*np.pi/popt[3]    
    Theta[key] = -np.pi/2 - popt[3]*popt[4]
    print(-np.pi/2 - popt[3]*popt[4])
    j_0 = max([max(fit_y[key]), abs(min(fit_y[key]))])
    #print(j_0)
    #print(2* np.pi / popt[3])
    #print(f"{key} has IPL: {2*np.pi/popt[3]}")
    wavelength[key] = 2* np.pi / popt[3]

# For j0 lets just take the largest value of the fitting fit_y



fig, axes = plt.subplots(3,3,figsize=(16, 8))
i = 0
j = 0
for keys in np.sort([*currents_detrend]):

    offset = 3.79

    def wedge2CEP(x):
        return  x/wavelength[keys] - offset

    def CEP2wedge(x):
        return x*wavelength[keys] + offset


    axes[i][j].tick_params(axis='x', direction='in')
    axes[i][j].tick_params(axis='y', direction='in')
    axes[i][j].tick_params(axis='both', which='major', labelsize=8)
    axes[i][j].tick_params(axis='both', which='minor', labelsize=16)

    
    if i==0:
        period_axis = axes[i][j].secondary_xaxis('top', functions=(wedge2CEP, CEP2wedge))
        period_axis.set_xlabel(r'Shift in CEP ($\pi$ radians)', fontsize=18)
        period_axis.tick_params(axis='x', direction='in')
        period_axis.tick_params(axis='both', which='major', labelsize=14)
        period_axis.tick_params(axis='both', which='minor', labelsize=12)

    
    tick_locations = period_axis.get_xticks()

    # 4. Iterate and draw a vertical line at each tick location
    for tick in tick_locations:
        pass#axes[i][j].axvline(x=tick, color='r', linestyle='--', linewidth=0.8) # 'r' for red color

    axes[i][j].plot(common.trans2ins(wedge_positions),currents_detrend[keys]*1e-7*1e12,"-dk", color="gray")
    axes[i][j].plot(fit_x[keys],fit_y[keys], color="purple")
    axes[i][j].tick_params(axis='both', which='major', labelsize=14)
    
    if j == 0:
        axes[i][j].set_ylabel(r'$j_{CEP}$ (pA)', fontsize = 18)
    if j != 0:
        axes[i][j].set_yticks([])
    if (i == 2) or (i == 1 and j != 0):
        axes[i][j].set_xlabel('Wedge insertion (mm)', fontsize = 18)
    if i != 2 and not (i == 1 and j != 0):
        axes[i][j].set_xticks([])

    axes[i][j].legend(title=str(np.round(keys, decimals=3)) + ' cm', loc='upper left', prop={'size': 16})

    axes[i][j].set_xlim(common.trans2ins(wedge_positions[0]), common.trans2ins(wedge_positions[-1]))
    axes[0][j].set_ylim(300, -300)
    axes[1][j].set_ylim(300, -300)
    axes[2][j].set_ylim(100, -100)

    if j == 2:
        i += 1
        j = 0
        continue
    if j != 2:
        j += 1
        continue
    
axes[2][1].axis('off')
axes[2][2].axis('off')

plt.subplots_adjust(wspace=0.04, hspace=0.14)
plt.savefig("./images/xscan.pdf")
plt.show()


ax = plt.gca()

average_ipl = 0.12

x_offset = 3.83*2

def wedge2CEP(x):
    return  2*x/(average_ipl) - x_offset

def CEP2wedge(x):
    return x*average_ipl/2 + x_offset

period_axis = ax.secondary_xaxis('top', functions=(wedge2CEP, CEP2wedge))
period_axis.tick_params(axis='x', direction='in')
period_axis.tick_params(axis='both', which='major', labelsize=14)
period_axis.tick_params(axis='both', which='minor', labelsize=12)
period_axis.set_xlabel(r'Shift in CEP ($\pi$ radians)', fontsize=18)


plt.plot(fit_x[1.1729999999999994],fit_y[1.1729999999999994], color="blue", label = "Right Contact")
plt.scatter(common.trans2ins(wedge_pos_data[1.1729999999999994]),currents_detrend_data[1.1729999999999994]*1e-7*1e12, marker="s", color="lightblue")
#plt.plot(fit_x[1.1679999999999995],fit_y[1.1679999999999995], color="gray")
plt.plot(fit_x[1.148],fit_y[1.148], color="red", label = "Left Contact")
plt.scatter(common.trans2ins(wedge_pos_data[1.148]),currents_detrend_data[1.148]*1e-7*1e12, marker="s", color="pink")
plt.xlabel('Wedge insertion (mm)', fontsize=18)
plt.ylabel(r'$ j_{CEP} $ (pA)', fontsize=18)
plt.xlim(min(fit_x[1.1729999999999994]), max(fit_x[1.1729999999999994]))
#plt.xlim(0, max(fit_x[1.1729999999999994]-x_offset))
#plt.grid()
plt.tight_layout()
plt.tick_params(axis='both', which='major', labelsize=14)
plt.legend(fontsize=18)
plt.savefig("./images/xscan_curves.pdf")
plt.show()