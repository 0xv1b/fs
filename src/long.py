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

mpl.rcParams['legend.title_fontsize'] = 'xx-large'

projection_phase = 116
powerscans = {} # Dictionary
powerscans_detrend = {}
wedge_pos_array = [] # List
wedge_positions = []
currents = {}
currents_detrend = {}
CEP_amplitudes = {}

for filename in constants.CEP_LONG:
        file = h5py.File(str(filename),'r')
        # Extract powerwheel angle from filename
        pos_filename_start = []
        pos_filename_end = []

        for alles in re.finditer('left_',filename):
            pos_filename_start.append(alles.end())
        for alles in re.finditer('deg_',filename):
            pos_filename_end.append(alles.start())
        
        pos = float(filename[pos_filename_start[0]:pos_filename_end[0]])

        if pos in powerscans.keys():
            pass
        else:
            powerscans[pos] = []# Adds pos key and assigns list
            powerscans_detrend[pos] = []# Adds pos key and assigns list

        wedge_pos_array.append(file['measurement']['Zman_script']['wedge_stage:position:sweep'][0]) # Adds list to list
        current = common.func_current(file['measurement']['Zman_script']['LockIn_3:MagPhase1'][0],file['measurement']['Zman_script']['LockIn_3:MagPhase1'][1],projection_phase)
        # Calculate current from LockIn and phase data
        powerscans[pos].append(current) # Adds currents list to list

wedge_positions = np.round(np.mean(wedge_pos_array,axis=0),1)

currents_detrend = {}
for key in powerscans.keys():
    currents[key] = np.mean(powerscans[key],axis=0)

for key in np.flip(np.sort([*currents])):
    x = np.arange(len(currents[key]))
    coefficients = np.polyfit(x, currents[key], 8)
    polyfit = np.polyval(coefficients, x)
    currents_detrend[key] = currents[key] - polyfit

    fig, ax = plt.subplots(figsize=(16, 8))
    ax.plot(x,currents[key]*1e-7*1e12, color="gray")
    ax.plot(x,polyfit*1e-7*1e12, color="red")
    #ax.plot(common.trans2ins(wedge_positions),currents_detrend[key]*1e-7*1e12, color="purple")
    plt.show()

d = 60 # This allows for basically exactly 5 peaks
pinit = [20,1,0.4,d,0.42,0]


#fit_x, fit_y = common.get_CEP_fit(wedge_positions, currents_detrend, pinit, 4000, x_min=0, x_max=1.35)
fit_x = np.linspace(0,1.35,1000)
fit_y = {}

for key in np.sort([*currents_detrend]):
    popt,pcov = curve_fit(common.CEP_fit,common.trans2ins(wedge_positions),currents_detrend[key]*1e-7*1e12,p0=pinit, maxfev=4000)
    fit_y[key] = common.CEP_fit(fit_x,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])
    print(common.trans2ins(wedge_positions))

wavelength = 2* np.pi / popt[3]
IPL = 2*np.pi/popt[3]    
Theta = ((-np.pi/2 - popt[3]*popt[4]) % (2*np.pi))/np.pi

print(f"IPL: {IPL}")
print(f"Theta: {Theta} pi")
print(f"c: {popt[5]}")

def wedge2CEP(x):
    return  x/wavelength

def CEP2wedge(x):
    return x*wavelength

for key in np.sort([*currents_detrend]):
    
    #CEP_amplitudes[key] = [round(abs(max(currents_detrend[key]*1e-8*1e12)-min(currents_detrend[key]*1e-8*1e12))/2,2),round(popt[0],2)]  

    fig, ax = plt.subplots(figsize=(16, 8))
    ax.tick_params(axis='x', direction='in')
    ax.tick_params(axis='y', direction='in')
    ax.tick_params(axis='both', which='major', labelsize=18)
    ax.tick_params(axis='both', which='minor', labelsize=12)

    x_offset = -0.0015
    ax.plot(common.trans2ins(wedge_positions)-x_offset,currents_detrend[key]*1e-7*1e12,'-dk', color="gray")
    ax.plot(fit_x-x_offset,fit_y[key], color="purple")

    plt.ylabel(r'$ j_{CEP} $ (pA)', fontsize = 18)
    plt.xlabel('Wedge insertion (mm)', fontsize = 18)

    period_axis = ax.secondary_xaxis('top', functions=(wedge2CEP, CEP2wedge))
    period_axis.set_xlabel(r'Shift in CEP ($\pi$ radians)', fontsize=18)
    period_axis.tick_params(axis='x', direction='in')
    period_axis.tick_params(axis='both', which='major', labelsize=18)
    period_axis.tick_params(axis='both', which='minor', labelsize=12)


    #plt.axvline(0, color='darkgrey', linestyle='--')
    #plt.axvline(CEP2wedge(10), color='darkgrey', linestyle='--')

    plt.legend(title=r'$\varepsilon_0$ = '+str(round(common.power2Efield(float(constants.power_powerCal_interpol(key))*1e-3,6.5e-15,5e-6)*1e-9,2)) + ' V/nm', loc='upper right', prop={'size': 16}, title_fontsize=24)
    plt.xlim(0, 1.35)
    plt.savefig("./images/long.pdf")
    plt.show()

    plt.close(fig)