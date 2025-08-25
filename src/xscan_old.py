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

mpl.rcParams['legend.title_fontsize'] = 'large'


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


fig = plt.figure(figsize=(15, 4))
length_keys = len([*currents])
for i,key in enumerate(np.flip(np.sort([*currents]))):
    plt.plot(common.trans2ins(wedge_positions),currents_detrend[key]*1e-8*1e12,'o-', label='xpos: '+str('{:.3f}'.format(round(key,3)))+'mm',color=[(length_keys-i)/length_keys,(length_keys-i)/length_keys,(i+1)/length_keys])

plt.xlabel('Wedge insertion (mm)',fontsize='14')
plt.ylabel('$ I_{CEP} $ (pA)',fontsize='14')
#plt.xticks([4.0, 5.0, 6.0, 7.0,8.0,9.0,10.0,11.0],fontsize='12')
#plt.yticks([-30,-20,-10,0,10,20,30],fontsize='12')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

margin = 0.005
plt.xlim(common.trans2ins(wedge_positions)[0]-margin, common.trans2ins(wedge_positions)[-1]+margin)
plt.grid()
plt.tight_layout()
#plt.show()
plt.close(fig)




CEP_amplitudes = {}

for keys in np.sort([*currents_detrend]):
    print(keys) 
    # Keys are 1.148,1.1529999999999998,1.1579999999999997,1.1629999999999996,1.1679999999999995,1.1729999999999994,1.1779999999999993

    try:
        
        
        popt,pcov = curve_fit(common.CEP_fit,common.trans2ins(wedge_positions),currents_detrend[keys]*1e-7*1e12,p0=pinit, maxfev = 4000)
        sopt,scov = curve_fit(common.sin_fit,common.trans2ins(wedge_positions),currents_detrend[keys]*1e-7*1e12,p0=pinit, maxfev = 4000)

        x_values = np.linspace(0.45,1.0,1000)
        #y_values = CEP_fit(x_values,pinit[0],pinit[1],pinit[2],pinit[3],pinit[4],pinit[5]) # rough estimate
        y_values = common.CEP_fit(x_values,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]) # fit

        y_values_sin = common.sin_fit(x_values,sopt[0],sopt[1],sopt[2],sopt[3],sopt[4],sopt[5]) 
        #1st value: min-max, 2nd value: fit amplitude

        # THIS PART WAS COMMENTED OUT BEFORE --------------


        CEP_amplitudes[keys] = [round(abs(max(currents_detrend[keys]*1e-8*1e12)-min(currents_detrend[keys]*1e-8*1e12))/2,2),round(popt[0],2)]  
        
        #Plots of CEP-currents with fits
        fig, ax = plt.subplots()
        ax.tick_params(axis='x', direction='in')
        ax.tick_params(axis='y', direction='in')
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.tick_params(axis='both', which='minor', labelsize=12)
        
        plt.plot(common.trans2ins(wedge_positions),currents_detrend[keys]*1e-7*1e12,'-dk', color="gray")
        plt.plot(x_values,y_values, color="purple", label="CEP fit")
        plt.plot(x_values,y_values_sin, color="orange", label="Sin fit")
        plt.xlabel('Wedge insertion in [mm]')
        plt.ylabel('$ I_{CEP} $ (pA)')
        plt.legend(title='x = '+str(keys)+' mm', loc='upper right', prop={'size': 16})
        plt.xlim(x_values[0], x_values[-1])

        
        
        
        #plt.show()
        plt.close(fig)

    except Exception as e:
        print("FAILED:" + e)



print(len([*currents_detrend]))

fig, axes = plt.subplots(3,3,figsize=(10, 4))
i = 0
j = 0
for keys in np.sort([*currents_detrend]):

    popt,pcov = curve_fit(common.CEP_fit,common.trans2ins(wedge_positions),currents_detrend[keys]*1e-7*1e12,p0=pinit, maxfev = 4000)
    x_values = np.linspace(0.45,1.0,1000)
    y_values = common.CEP_fit(x_values,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]) # fit

    axes[i][j].tick_params(axis='x', direction='in')
    axes[i][j].tick_params(axis='y', direction='in')
    axes[i][j].tick_params(axis='both', which='major', labelsize=8)
    axes[i][j].tick_params(axis='both', which='minor', labelsize=16)

    axes[i][j].plot(common.trans2ins(wedge_positions),currents_detrend[keys]*1e-7*1e12,"-dk", color="gray")
    axes[i][j].plot(x_values,y_values, color="purple")
    
    if j == 0:
        axes[i][j].set_ylabel(r'$ I_{CEP} $ (pA)', fontsize = 18)
    if j != 0:
        axes[i][j].set_yticks([])
    if (i == 2) or (i == 1 and j != 0):
        axes[i][j].set_xlabel('Wedge insertion (mm)', fontsize = 18)
    if i != 2 and not (i == 1 and j != 0):
        axes[i][j].set_xticks([])

    axes[i][j].legend(title=str(np.round(keys, decimals=3)) + ' cm', loc='upper left', prop={'size': 16})

    axes[i][j].set_xlim(x_values[0], x_values[-1])
    #axes[i][j].set_xlim(0, 1.35)

    if j == 2:
        i += 1
        j = 0
        continue
    if j != 2:
        j += 1
        continue
    
axes[2][1].axis('off')
axes[2][2].axis('off')

plt.subplots_adjust(wspace=0.04, hspace=0.08)
plt.show()
