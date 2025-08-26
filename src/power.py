import numpy as np
import matplotlib.pyplot as  plt
import matplotlib as mpl
import h5py
import glob
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
import re
import os
from .helpers import common, detrend, constants, fits
import warnings
warnings.filterwarnings("ignore")


mpl.rcParams['legend.title_fontsize'] = 'large'


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

currents_detrend_data = {}
wedge_pos_data = {}
start = [1,1,1,1,1,1,1,1,1,1,1,1]
end =   [1,1,1,1,1,1,1,1,1,1,1,1]
for i, key in enumerate(np.sort([*currents_detrend])):
    currents_detrend_data[key] = currents_detrend[key][start[i]:-end[i]]
    wedge_pos_data[key] = wedge_positions[start[i]:-end[i]]



fit_x = {}
fit_y = {}
wavelength = {}
IPL = {}
Theta = {}
j_0 = {}
for i, key in enumerate(np.sort([*currents_detrend])):
    fit_x[key] = np.linspace(common.trans2ins(wedge_pos_data[key][0]),common.trans2ins(wedge_pos_data[key][-1]),1000)
    popt,pcov = curve_fit(common.CEP_fit,common.trans2ins(wedge_pos_data[key]),currents_detrend_data[key]*1e-7*1e12,p0=pinit, maxfev=40000)
    fit_y[key] = common.CEP_fit(fit_x[key],popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])
    IPL[key] = 2*np.pi/popt[3]    
    Theta[key] = -np.pi/2 - popt[3]*popt[4]
    j_0[key] = max([max(fit_y[key]), abs(min(fit_y[key]))])
    wavelength[key] = 2* np.pi / popt[3]
    CEP_amplitudes[key] = [round(abs(max(currents_detrend[key]*1e-8*1e12)-min(currents_detrend[key]*1e-8*1e12))/2,2),round(popt[0],2)]  

# For j0 lets just take the largest value of the fitting fit_y



fig, axes = plt.subplots(4,3,figsize=(10, 4))
i = 0
j = 0
for keys in np.flip(np.sort([*currents_detrend])):

    def wedge2CEP(x):
        return  x/wavelength[keys]

    def CEP2wedge(x):
        return x*wavelength[keys]


    axes[i][j].tick_params(axis='x', direction='in')
    axes[i][j].tick_params(axis='y', direction='in')
    axes[i][j].tick_params(axis='both', which='major', labelsize=8)
    axes[i][j].tick_params(axis='both', which='minor', labelsize=16)

    
    period_axis = axes[i][j].secondary_xaxis('top', functions=(wedge2CEP, CEP2wedge))
    period_axis.set_xlabel(r'Shift in CEP ($\pi$ radians)', fontsize=8)
    period_axis.tick_params(axis='x', direction='in')
    period_axis.tick_params(axis='both', which='major', labelsize=18)
    period_axis.tick_params(axis='both', which='minor', labelsize=12)

    
    tick_locations = period_axis.get_xticks()

    # 4. Iterate and draw a vertical line at each tick location
    for tick in tick_locations:
        axes[i][j].axvline(x=tick, color='r', linestyle='--', linewidth=0.8) # 'r' for red color

    axes[i][j].plot(common.trans2ins(wedge_positions),currents_detrend[keys]*1e-7*1e12,"-dk", color="gray") # 1e-7 is correct
    axes[i][j].plot(fit_x[keys],fit_y[keys], color="purple")
    
    if j == 0:
        axes[i][j].set_ylabel(r'$ I_{CEP} $ (pA)', fontsize = 18)
    if j != 0:
        axes[i][j].set_yticks([])
    if (i == 2) or (i == 1 and j != 0):
        axes[i][j].set_xlabel('Wedge insertion (mm)', fontsize = 18)
    if i != 2 and not (i == 1 and j != 0):
        axes[i][j].set_xticks([])

    axes[i][j].legend(title=str(round(common.power2Efield(float(constants.power_powerCal_interpol(keys))*1e-3,6.5e-15,5e-6)*1e-9,2)) + ' V/nm', loc='upper right', prop={'size': 16})

    axes[i][j].set_xlim(common.trans2ins(wedge_positions[0]), common.trans2ins(wedge_positions[-1]))
    axes[0][j].set_ylim(50, -50)
    axes[1][j].set_ylim(100, -100)
    axes[2][j].set_ylim(200, -200)
    axes[3][j].set_ylim(300, -300)

    if j == 2:
        i += 1
        j = 0
        continue
    if j != 2:
        j += 1
        continue
    

#plt.subplots_adjust(wspace=0.04, hspace=0.08)
plt.show()


beam_diameter = 5e-6
pulse_FWHM = 6.5e-15

CEP_current_minmax = []
CEP_current_fit = []
E_field = []

for keys in np.flip(np.sort([*currents_detrend])):
    CEP_current_minmax.append(CEP_amplitudes[keys][0])
    CEP_current_fit.append(CEP_amplitudes[keys][1])
    E_field.append(common.power2Efield(float(constants.power_powerCal_interpol(keys))*1e-3,pulse_FWHM,beam_diameter)*1e-9)
    
# plt.plot(E_field,CEP_current_minmax,marker='o')
# plt.xlabel('E-field in [V/nm]')
# plt.ylabel('CEP-dependent current in [pA]')
# plt.title('CEP-dependent current - Min/Max')
# #plt.yscale('log')
# #plt.xscale('log')
# plt.grid()
# plt.tight_layout()
# plt.show()

# plt.plot(E_field,list(map(abs, CEP_current_fit)),marker='o')
# plt.xlabel('E-field in [V/nm]')
# plt.ylabel('CEP-dependent current in [pA]')
# plt.title('CEP-dependent current - Fit')
# plt.grid()
# plt.tight_layout()
# plt.show()


pinit = [1.3, 1.6]
popt,pcov = curve_fit(fits.linear,np.log(E_field),np.log(CEP_current_minmax),p0=pinit)
x_values = np.linspace(min(E_field),max(E_field))
y_values = fits.linear(np.log(x_values),pinit[0],pinit[1]) # fit
#y_values = linear(np.log(x_values),popt[0],popt[1]) # fit

rect = plt.Rectangle((1.4, 3.8), 0.5, 1.8,
                     facecolor="red", alpha=0.2)
ax = plt.gca()
ax.add_patch(rect)

plt.plot(E_field,CEP_current_minmax,'o:')
plt.plot(x_values, np.exp(y_values))
plt.xlabel('Peak electric field (V/nm)')
plt.ylabel(r'$ I_{CEP} $ (pA)')
plt.yscale('log')
plt.xscale('log')
plt.grid()
plt.tight_layout()
plt.show()