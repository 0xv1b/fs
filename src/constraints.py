import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as  plt
import h5py
import glob
from scipy.optimize import curve_fit, leastsq
from scipy.interpolate import CubicSpline
import re
from pprint import pp
#from obspy.signal.detrend import polynomial
from .helpers import common, detrend, constants
from matplotlib.widgets import Slider


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

# Trim wedge_positions and currents_detrend for outliers
for key in np.sort([*currents_detrend]):
    pass

d = 60 # This allows for basically exactly 5 peaks between 0.45 annd 1.0, below 50 it becomes 4
pinit = [20,1,0.4,d,0.42,0]


fit_x, fit_y = common.get_CEP_fit(wedge_positions, currents_detrend, pinit, 4000)
sin_x, sin_y = common.get_sin_fit(wedge_positions, currents_detrend, pinit, 4000)


# Fitting with a constraint such that d cannot go below 50
def residuals(p, x, y):
    sigma_min_desirable = 0.01
    sigma_max_desirable = 0.2
    penalization = 0
    if p[3] > 70:
        penalization += 10
    if p[3] < 50:
        penalization += 10

    if p[2] < sigma_min_desirable:
        penalization += 10  
    if p[2] > sigma_max_desirable:
        penalization += 10 
    
    if np.abs(p[5]) > 0.01:
        penalization += 10
    
    fit_value = common.CEP_fit(x, p[0], p[1], p[2], p[3], p[4], p[5])

    return y - fit_value - penalization*10



cons_x = np.linspace(0.45,1.0,1000)
cons_y = {}
for key in np.sort([*currents_detrend]):
    copt, ccov = leastsq(func=residuals, x0=pinit, args=(common.trans2ins(wedge_positions),currents_detrend[key]*1e-7*1e12))
    cons_y[key] = common.CEP_fit(cons_x,copt[0],copt[1],copt[2],copt[3],copt[4],copt[5])
    #print(' ' + str(round(common.power2Efield(float(constants.power_powerCal_interpol(key))*1e-3,6.5e-15,5e-6)*1e-9,2)) + ': ' + str(copt))


for key in np.sort([*currents_detrend]):

    # TODO: fix this CEP amp shit
    #CEP_amplitudes[key] = [round(abs(max(currents_detrend[key]*1e-8*1e12)-min(currents_detrend[key]*1e-8*1e12))/2,2),round(popt[0],2)]  

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.tick_params(axis='x', direction='in')
    ax.tick_params(axis='y', direction='in')
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.tick_params(axis='both', which='minor', labelsize=12)

    ax.plot(common.trans2ins(wedge_positions),currents_detrend[key]*1e-7*1e12,'-dk', color="gray")
    #ax.plot(fit_x,fit_y[key], color="purple", label="CEP fit")
    #ax.plot(sin_x,sin_y[key], color="orange", label = "Sin fit")
    #ax.plot(cons_x,cons_y[key], color="green", label = "Constrained CEP fit")


    line, = ax.plot(fit_x, common.CEP_fit(fit_x, 90,0.75,copt[2],copt[3],copt[4], 0), lw=1, linestyle="-")

    plt.ylabel(r'$ I_{CEP} $ (pA)', fontsize = 18)
    plt.xlabel('Wedge insertion (mm)', fontsize = 18)
    plt.axhline(0, color='darkgrey', linestyle='--')

    plt.legend(title=str(round(common.power2Efield(float(constants.power_powerCal_interpol(key))*1e-3,6.5e-15,5e-6)*1e-9,2)) + ' V/nm', loc='upper right', prop={'size': 16})
    plt.xlim(0, 1.35)

    fig.subplots_adjust(left=0.25, bottom=0.25)

    def update(val):
        line.set_ydata(common.CEP_fit(fit_x, amp_slider.val, offset_slider.val, gauss_slider.val,wavelength_slider.val,copt[4], 0))
        fig.canvas.draw_idle()
    


    axamp = fig.add_axes([0.25, 0.14, 0.65, 0.03])
    amp_slider = Slider(
        ax=axamp,
        label='Amplitude',
        valmin=0.1,
        valmax=100,
        valinit=90,
    )

    axoffset = fig.add_axes([0.25, 0.1, 0.65, 0.03])
    offset_slider = Slider(
        ax=axoffset,
        label='Offset',
        valmin=0.1,
        valmax=2,
        valinit=0.75,
    )

    init_wavelength = copt[3]
    axwavelength = fig.add_axes([0.25, 0.06, 0.65, 0.03])
    wavelength_slider = Slider(
        ax=axwavelength,
        label='Wavelength',
        valmin=0.01,
        valmax=100,
        valinit=init_wavelength,
    )

    init_gauss = copt[2]
    axgauss = fig.add_axes([0.25, 0.02, 0.65, 0.03])
    gauss_slider = Slider(
        ax=axgauss,
        label='Gauss',
        valmin=0.01,
        valmax=1,
        valinit=init_gauss,
    )


    offset_slider.on_changed(update)
    gauss_slider.on_changed(update)
    wavelength_slider.on_changed(update)
    amp_slider.on_changed(update)

    # Only offset and amplitude are initially different

    plt.show()
    plt.close(fig)




# Actual format

fig, axes = plt.subplots(int(len([*currents_detrend])/3),3,figsize=(10, 4))
i = 0
j = 0
for key in np.sort([*currents_detrend]):
    
    axes[i][j].tick_params(axis='x', direction='in')
    axes[i][j].tick_params(axis='y', direction='in')
    axes[i][j].tick_params(axis='both', which='major', labelsize=8)
    axes[i][j].tick_params(axis='both', which='minor', labelsize=16)

    axes[i][j].plot(common.trans2ins(wedge_positions),currents_detrend[key]*1e-7*1e12,"-dk", color="gray")
    axes[i][j].plot(fit_x,fit_y[key], color="purple")
    # axes[i][j].plot(sin_x,sin_y[key], color="orange", label = "Sin fit")
    # axes[i][j].plot(cons_x,cons_y[key], color="green", label = "Constrained CEP fit")
    if j == 0:
        axes[i][j].set_ylabel(r'$ I_{CEP} $ (pA)', fontsize = 18)
    if j != 0:
        axes[i][j].set_yticks([])
    if i == 3:
        axes[i][j].set_xlabel('Wedge insertion (mm)', fontsize = 18)
    if i != 3:
        axes[i][j].set_xticks([])

    axes[i][j].legend(title=str(round(common.power2Efield(float(constants.power_powerCal_interpol(key))*1e-3,6.5e-15,5e-6)*1e-9,2)) + ' V/nm', loc='upper left', prop={'size': 16})
    axes[i][j].set_xlim(0, 1.35)

    if j == 2:
        i += 1
        j = 0
        continue
    if j != 2:
        j += 1
        continue
    


plt.subplots_adjust(wspace=0.04, hspace=0.08)
plt.show()
