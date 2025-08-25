import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as  plt
import h5py
import glob
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
import re
#from obspy.signal.detrend import polynomial
from .helpers import common, detrend, constants

projection_phase = 116
powerscans = {} # Dictionary
powerscans_detrend = {}
wedge_pos_array = [] # List
wedge_positions = []
currents = {}
currents_detrend = {}
CEP_amplitudes = {}

detrend.init_lists(constants.CEP_POWER, powerscans, powerscans_detrend, wedge_pos_array, projection_phase=projection_phase)
wedge_positions = detrend.get_wedge_positions(wedge_pos_array)
currents_detrend = detrend.get_detrended_currents(powerscans, currents)


d = 60 # This allows for basically exactly 5 peaks
pinit = [20,1,0.4,d,0.42,0]


# Data plot

fig, ax = plt.subplots()
ax.tick_params(axis='x', direction='in')
ax.tick_params(axis='y', direction='in')
ax.tick_params(axis='both', which='major', labelsize=14)

# y = round(common.power2Efield(float(constants.power_powerCal_interpol(constants.angle_powerCal_interpol))*1e-3,6.5e-15,5e-6)*1e-9,2)

plt.plot(constants.angle_powerCal,common.power2Efield(constants.power_powerCal*1e3*1e-3,6.5e-15,5e-6)*1e-9,linewidth='2',label='Data')
plt.plot(constants.angle_powerCal_interpol,common.power2Efield(constants.power_powerCal_interpol(constants.angle_powerCal_interpol)*1e-3,6.5e-15,5e-6)*1e-9,linewidth='2',label='Fit')
plt.xlabel('ND filter wheel angle (°)', fontsize=14)
plt.ylabel('Peak electric field (V/nm)', fontsize=14)
plt.legend(fontsize=16)
#plt.grid()
plt.tight_layout()

plt.show()





fit_x, fit_y = common.get_CEP_fit(wedge_positions, currents_detrend, pinit, 4000)

for key in np.flip(np.sort([*currents])):
    fig, axes = plt.subplots(3, 1, figsize=(10, 4))
    for ax in axes:
        ax.tick_params(axis='x', direction='in')
        ax.tick_params(axis='y', direction='in')
    
    x = np.arange(len(currents[key]))

    axes[0].plot(x, currents[key], label='Original Data', color="purple")
    axes[1].plot(x, currents_detrend[key], label='Detrended Data', color="purple")

    axes[2].tick_params(axis='x', direction='in')
    axes[2].tick_params(axis='y', direction='in')
    axes[2].plot(common.trans2ins(wedge_positions),currents_detrend[key]*1e-7*1e12, color="purple")
    axes[2].plot(fit_x,fit_y[key], color="gray")
    
    plt.legend()
    #plt.show()
    plt.close(fig)


for key in np.sort([*currents_detrend]):
    
    #CEP_amplitudes[key] = [round(abs(max(currents_detrend[key]*1e-8*1e12)-min(currents_detrend[key]*1e-8*1e12))/2,2),round(popt[0],2)]  

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.tick_params(axis='x', direction='in')
    ax.tick_params(axis='y', direction='in')
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.tick_params(axis='both', which='minor', labelsize=12)

    ax.plot(common.trans2ins(wedge_positions),currents_detrend[key]*1e-7*1e12, color="gray", marker="o", markersize=0.8)
    ax.plot(fit_x,fit_y[key], color="purple")

    plt.ylabel(r'$ I_{CEP} $ (pA)', fontsize = 18)
    plt.xlabel('Wedge insertion (mm)', fontsize = 18)

    period_axis = ax.secondary_xaxis('top')
    period_axis.set_xlabel('angle [rad]')
    period_axis.tick_params(axis='x', direction='in')
    period_axis.tick_params(axis='both', which='major', labelsize=14)
    period_axis.tick_params(axis='both', which='minor', labelsize=12)


    plt.legend(title=str(round(common.power2Efield(float(constants.power_powerCal_interpol(key))*1e-3,6.5e-15,5e-6)*1e-9,2)) + ' V/nm', loc='upper right', prop={'size': 16})
    plt.xlim(0, 1.35)
    #plt.show()

    plt.close(fig)