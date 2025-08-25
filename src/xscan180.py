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


xscan180_keys = [1.14,1.142,1.144,1.146,1.148,1.15,1.152,1.155,1.16,1.1649999999999998,1.1699999999999997,1.17,1.1749999999999998]

projection_phase = 113
xscans = {}
xscans_detrend = {}
currents_detrend = {}
wedge_pos_array = []
currents = {}

detrend.init_lists(constants.CEP_XSCAN_180, xscans, xscans_detrend, wedge_pos_array, projection_phase=projection_phase, xscan=True)
wedge_positions = detrend.get_wedge_positions(wedge_pos_array, decimals=2)
currents_detrend = detrend.get_detrended_currents(xscans, currents, deg=7)

d = 60 # This allows for basically exactly 5 peaks
pinit = [20,1,0.4,d,0.42,0]

currents_detrend_data = {}
wedge_pos_data = {}
start = [1,1,10,1,20,1,90,1,70,1 ,5 ,1,1]
end =   [1,1,60,1,20,1,20,1,1 ,90,70,1,1]
for i, key in enumerate(np.sort([*currents_detrend])):
    currents_detrend_data[key] = currents_detrend[key][start[i]:-end[i]]
    wedge_pos_data[key] = wedge_positions[start[i]:-end[i]]


# ---------------------------------------------------------------------------- #
#                                    Fitting                                   #
# ---------------------------------------------------------------------------- #
fit_x = {}
fit_y = {}
wavelength = {}
IPL = {}
Theta = {}
for i,key in enumerate(np.sort([*currents_detrend])):
    fit_x[key] = np.linspace(common.trans2ins(wedge_pos_data[key][0]),common.trans2ins(wedge_pos_data[key][-1]),1000)
    popt,pcov = curve_fit(common.CEP_fit,common.trans2ins(wedge_pos_data[key]),currents_detrend_data[key]*1e-7*1e12,p0=pinit, maxfev=40000)
    fit_y[key] = common.CEP_fit(fit_x[key],popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])
    #print(f"{key}:  {popt}")
    IPL[key] = 2*np.pi/popt[3]
    Theta[key] = -np.pi/2 - popt[3]*popt[4]
    #print(-np.pi/2 - popt[3]*popt[4])
    #print(popt[0])
    j_0 = max([max(fit_y[key]), abs(min(fit_y[key]))])
    print(j_0)
    #print(2* np.pi / popt[3])
    #print(f"{key} has IPL: {2*np.pi/popt[3]}")
    wavelength[key] = 2* np.pi / popt[3]

# -pi/2 - D*E 

# def wedge2CEP(x):
#     return  x/wavelength[1.1629999999999996]

# def CEP2wedge(x):
#     return x*wavelength[1.1629999999999996]


# fig, ax = plt.subplots()
# ax.tick_params(axis='x', direction='in')
# ax.tick_params(axis='y', direction='in')
# ax.tick_params(axis='both', which='major', labelsize=14)

# period_axis = ax.secondary_xaxis('top', functions=(wedge2CEP, CEP2wedge))
# period_axis.set_xlabel(r'Shift in CEP ($\pi$ radians)', fontsize=18)
# period_axis.tick_params(axis='x', direction='in')
# period_axis.tick_params(axis='both', which='major', labelsize=18)
# period_axis.tick_params(axis='both', which='minor', labelsize=12)

# ax.plot(common.trans2ins(wedge_positions),currents_detrend[key]*1e-7*1e12,"-dk", color="gray")
# ax.plot(fit_x[1.1629999999999996],fit_y[1.1629999999999996], color="purple")
# plt.xlabel('Wedge insertion (mm)', fontsize=18)
# plt.ylabel(r'$ I_{CEP} $ (pA)', fontsize=18)
# plt.tight_layout()
# plt.show()


fig, axes = plt.subplots(5,3,figsize=(10, 4))
i = 0
j = 0
for keys in np.sort([*currents_detrend]):

    def wedge2CEP(x):
        return  x/wavelength[keys]

    def CEP2wedge(x):
        return x*wavelength[keys]


    axes[i][j].tick_params(axis='x', direction='in')
    axes[i][j].tick_params(axis='y', direction='in')
    axes[i][j].tick_params(axis='both', which='major', labelsize=8)
    axes[i][j].tick_params(axis='both', which='minor', labelsize=16)

    
    period_axis = axes[i][j].secondary_xaxis('top', functions=(wedge2CEP, CEP2wedge))
    
    period_axis.tick_params(axis='x', direction='in')
    period_axis.tick_params(axis='both', which='major', labelsize=12)
    period_axis.tick_params(axis='both', which='minor', labelsize=12)

    tick_locations = period_axis.get_xticks()

    # 4. Iterate and draw a vertical line at each tick location
    for tick in tick_locations:
        axes[i][j].axvline(x=tick, color='r', linestyle='--', linewidth=0.8) # 'r' for red color


    axes[i][j].plot(common.trans2ins(wedge_positions),currents_detrend[keys]*1e-7*1e12,"-dk", color="gray")
    axes[i][j].plot(fit_x[keys],fit_y[keys], color="purple")
    
    if i == 0:
        period_axis.set_xlabel(r'Shift in CEP ($\pi$ radians)', fontsize=8)
    if j == 0:
        axes[i][j].set_ylabel(r'$ I_{CEP} $ (pA)', fontsize = 18)
    if j != 0:
        axes[i][j].set_yticks([])
    if (i == 4) or (i == 3 and j != 0):
        axes[i][j].set_xlabel('Wedge insertion (mm)', fontsize = 18)
    if i != 4 and not (i == 3 and j != 0):
        axes[i][j].set_xticks([])

    axes[i][j].legend(title=str(np.round(keys, decimals=3)) + ' cm', loc='upper left', prop={'size': 16})

    axes[i][j].set_xlim(common.trans2ins(wedge_positions[0]), common.trans2ins(wedge_positions[-1]))
    
    axes[0][j].set_ylim(50, -50)
    axes[1][j].set_ylim(50, -50)
    axes[2][j].set_ylim(50, -50)
    axes[3][j].set_ylim(300, -300)
    axes[4][j].set_ylim(50, -50)

    if j == 2:
        i += 1
        j = 0
        continue
    if j != 2:
        j += 1
        continue
    
axes[4][1].axis('off')
axes[4][2].axis('off')

#plt.subplots_adjust(wspace=0.04, hspace=0.08)
plt.show()

ax = plt.gca()

average_ipl = 0.113
def wedge2CEP(x):
    return  x/average_ipl

def CEP2wedge(x):
    return x*average_ipl

period_axis = ax.secondary_xaxis('top', functions=(wedge2CEP, CEP2wedge))
period_axis.tick_params(axis='x', direction='in')
period_axis.tick_params(axis='both', which='major', labelsize=12)
period_axis.tick_params(axis='both', which='minor', labelsize=12)

plt.plot(fit_x[1.17],fit_y[1.17], color="blue")
plt.plot(fit_x[1.155],fit_y[1.155], color="gray")
plt.plot(fit_x[1.146],fit_y[1.146], color="red")
plt.xlabel('Wedge insertion (mm)')
plt.ylabel(r'$ I_{CEP} $ (pA)')
plt.xlim(min(fit_x[1.17]), max(fit_x[1.17]))
plt.grid()
plt.tight_layout()
plt.show()