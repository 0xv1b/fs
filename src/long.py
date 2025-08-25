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

detrend.init_lists(constants.CEP_LONG, powerscans, powerscans_detrend, wedge_pos_array, projection_phase=projection_phase)
wedge_positions = detrend.get_wedge_positions(wedge_pos_array)
currents_detrend = detrend.get_detrended_currents(powerscans, currents)


d = 60 # This allows for basically exactly 5 peaks
pinit = [20,1,0.4,d,0.42,0]


#fit_x, fit_y = common.get_CEP_fit(wedge_positions, currents_detrend, pinit, 4000, x_min=0, x_max=1.35)
fit_x = np.linspace(0,1.35,1000)
fit_y = {}

for key in np.sort([*currents_detrend]):
    popt,pcov = curve_fit(common.CEP_fit,common.trans2ins(wedge_positions),currents_detrend[key]*1e-7*1e12,p0=pinit, maxfev=4000)
    fit_y[key] = common.CEP_fit(fit_x,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])

wavelength = 2* np.pi / popt[3]

def wedge2CEP(x):
    return  x/wavelength

def CEP2wedge(x):
    return x*wavelength

for key in np.sort([*currents_detrend]):
    
    #CEP_amplitudes[key] = [round(abs(max(currents_detrend[key]*1e-8*1e12)-min(currents_detrend[key]*1e-8*1e12))/2,2),round(popt[0],2)]  

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.tick_params(axis='x', direction='in')
    ax.tick_params(axis='y', direction='in')
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.tick_params(axis='both', which='minor', labelsize=12)

    x_offset = 0.02
    ax.plot(common.trans2ins(wedge_positions)-x_offset,currents_detrend[key]*1e-7*1e12,'-dk', color="gray")
    ax.plot(fit_x-x_offset,fit_y[key], color="purple")

    plt.ylabel(r'$ I_{CEP} $ (pA)', fontsize = 18)
    plt.xlabel('Wedge insertion (mm)', fontsize = 18)

    period_axis = ax.secondary_xaxis('top', functions=(wedge2CEP, CEP2wedge))
    period_axis.set_xlabel(r'Shift in CEP ($\pi$ radians)', fontsize=18)
    period_axis.tick_params(axis='x', direction='in')
    period_axis.tick_params(axis='both', which='major', labelsize=18)
    period_axis.tick_params(axis='both', which='minor', labelsize=12)


    #plt.axvline(0, color='darkgrey', linestyle='--')
    #plt.axvline(CEP2wedge(10), color='darkgrey', linestyle='--')

    plt.legend(title=str(round(common.power2Efield(float(constants.power_powerCal_interpol(key))*1e-3,6.5e-15,5e-6)*1e-9,2)) + ' V/nm', loc='upper right', prop={'size': 16})
    plt.xlim(-x_offset, 1.35)
    plt.show()

    plt.close(fig)