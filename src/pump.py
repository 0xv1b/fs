import numpy as np
import matplotlib.pyplot as  plt
import h5py
import glob
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
import re
from .helpers import common, detrend, constants
import warnings
warnings.filterwarnings("ignore")



# -------------- MIDDLE ---------------

# Initialisation -------------------
# ---------------------------- NOT PUMPED ----------------------------------
projection_phase = 113
trace = []
wedge_pos_array = []
currents = []
currents_detrend = []
detrend.init_pump(constants.CEP_MIDDLE_WO, wedge_pos_array, trace, projection_phase)
wedge_pos = detrend.get_wedge_positions(wedge_pos_array)
currents_detrend = detrend.get_detrended_currents_pump(currents, trace, deg=8)


# --------------- PUMPED -----------------
projection_phase = 113

trace_green = []
wedge_pos_array_green = []
currents_green = []
currents_detrend_green = []
detrend.init_pump(constants.CEP_MIDDLE_GREEN, wedge_pos_array_green, trace_green, projection_phase)
wedge_pos_green = detrend.get_wedge_positions(wedge_pos_array_green)
currents_detrend_green = detrend.get_detrended_currents_pump(currents_green, trace_green)
detrend.pump_sanity_check(wedge_pos_array_green, expected_len=102)
# Initialisation -------------------


# Data ---------------------------------------------
pinit = [100,3,0.9,55,0.42,0]
popt,pcov = curve_fit(common.CEP_fit,common.trans2ins(wedge_pos),currents_detrend*1e-8*1e12,p0=pinit, maxfev=4000)


x_values = np.linspace(0.45,1.05,1000)
y_values = common.CEP_fit(x_values,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]) # fit

ipl = 2*np.pi/popt[3]   
theta = ((-np.pi/2 - popt[3]*popt[4]) % (2*np.pi))/(2*np.pi)
j_0 = max([max(y_values), abs(min(y_values))])
print(f"IPL: {ipl}\nTheta: {theta}\nj_0: {j_0}\n")
#y_values = common.CEP_fit(x_values,100,3,1,55,0.42,0) # fit

#green pumped
pinit_green = [1.3,0.7,0.3,57,0.42,0]
popt_green,pcov_green = curve_fit(common.CEP_fit,common.trans2ins(wedge_pos_green),currents_detrend_green*1e-8*1e12,p0=pinit_green, maxfev=4000)
    
x_values_green = np.linspace(0.45,1.05,1000)
y_values_green = common.CEP_fit(x_values_green,popt_green[0],popt_green[1],popt_green[2],popt_green[3],popt_green[4],popt_green[5]) # fit

ipl_green = 2*np.pi/popt[3]   
theta_green =  ((-np.pi/2 - popt[3]*popt[4]) % (2*np.pi))/(2*np.pi)
j_0_green = max([max(y_values), abs(min(y_values))])
print(f"IPL Green: {ipl_green}\nTheta Green: {theta_green}\nj_0 Green: {j_0_green}")
# Data ---------------------------------------------

print("\n\n\n\n")

#Plots of CEP-currents with fits

fig, axes = plt.subplots(2,1)

#plt.title('Photocurrent (532nm pumped/not pumped) @ middle')

axes[0].plot(x_values,y_values,'-.b', label="Unpumped Fit", color="purple")# Dashed magenta
axes[0].plot(common.trans2ins(wedge_pos),currents_detrend*1e-8*1e12,'-dk', label="Unpumped Data", color="gray") # Black Diamond

axes[1].plot(x_values_green,y_values_green, '-.r', label="Pumped Fit", color="purple") # Dashed red
axes[1].plot(common.trans2ins(wedge_pos_green),currents_detrend_green*1e-8*1e12,'-og', label="Pumped Data", color="gray") # Green circles


for ax in axes:
    ax.tick_params(axis='x', direction='in')
    ax.tick_params(axis='y', direction='in')
    ax.set_xlabel('Wedge insertion (mm)', fontsize=18)
    ax.set_ylabel(r'$ I_{CEP} $ (pA)', fontsize=18)
    ax.legend( loc='upper right', prop={'size': 16})
    ax.set_xlim(x_values[0], x_values[-1])
    ax.tick_params(axis='both', which='major', labelsize=14)

plt.show()







# RIGHT CONTACT 

# Initialisation -------------------
# NOT PUMPED
projection_phase = 116
trace = []
wedge_pos_array = []
currents = []
currents_detrend = []

detrend.init_pump(constants.CEP_RIGHT_WO, wedge_pos_array, trace, projection_phase)
wedge_pos = detrend.get_wedge_positions(wedge_pos_array)
currents_detrend = detrend.get_detrended_currents_pump(currents, trace)



# PUMPED
trace_green = []
wedge_pos_array_green = []
currents_green = []
currents_detrend_green = []
detrend.init_pump(constants.CEP_RIGHT_GREEN, wedge_pos_array_green, trace_green, projection_phase)
wedge_pos_green = detrend.get_wedge_positions(wedge_pos_array_green)
currents_detrend_green = detrend.get_detrended_currents_pump(currents_green, trace_green)
detrend.pump_sanity_check(wedge_pos_array_green, expected_len=102)
# Initialisation -------------------


# Data -------------------------------------
pinit = [43,1,0.4,55,0.42,0]
popt,pcov = curve_fit(common.CEP_fit,common.trans2ins(wedge_pos),currents_detrend*1e-7*1e12,p0=pinit, maxfev=4000)
    
x_values = np.linspace(0.45,1.05,1000)
y_values = common.CEP_fit(x_values,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]) # fit

ipl = 2*np.pi/popt[3]   
theta =  ((-np.pi/2 - popt[3]*popt[4]) % (2*np.pi))/(2*np.pi)
j_0 = max([max(y_values), abs(min(y_values))])
print(f"IPL: {ipl}\nTheta: {theta}\nj_0: {j_0}\n")

#green pumped
pinit_green = [43,1,0.4,55,0.42,0]
popt_green,pcov_green = curve_fit(common.CEP_fit,common.trans2ins(wedge_pos_green),currents_detrend_green*1e-7*1e12,p0=pinit_green, maxfev=4000)
    
x_values_green = np.linspace(0.45,1.05,1000)
y_values_green = common.CEP_fit(x_values_green,popt_green[0],popt_green[1],popt_green[2],popt_green[3],popt_green[4],popt_green[5]) # fit

ipl_green = 2*np.pi/popt[3]   
theta_green =  ((-np.pi/2 - popt[3]*popt[4]) % (2*np.pi))/(2*np.pi)
j_0_green = max([max(y_values), abs(min(y_values))])
print(f"IPL Green: {ipl_green}\nTheta Green: {theta_green}\nj_0 Green: {j_0_green}")
# Data -------------------------------------  



#Plots of CEP-currents with fits ---------------------
fig, axes = plt.subplots(2,1)

#plt.title('Photocurrent (532nm pumped/not pumped) @contact')

axes[0].plot(x_values,y_values,'-.m', label="Unpumped Fit", color="purple")# Dashed magenta
axes[0].plot(common.trans2ins(wedge_pos),currents_detrend*1e-7*1e12,'-dk', label="Unpumped Data", color="gray") # Black Diamond

axes[1].plot(x_values_green,y_values_green, '-.r', label="Pumped Fit", color="purple") # Dashed red
axes[1].plot(common.trans2ins(wedge_pos_green),currents_detrend_green*1e-7*1e12,'-og', label="Pumped Data", color="gray") # Green circles


for ax in axes:
    ax.tick_params(axis='x', direction='in')
    ax.tick_params(axis='y', direction='in')
    ax.set_xlabel('Wedge insertion (mm)', fontsize=18)
    ax.set_ylabel(r'$ I_{CEP} $ (pA)', fontsize=18)
    ax.legend( loc='upper right', prop={'size': 16})
    ax.set_xlim(x_values[0], x_values[-1])

plt.show()
#Plots of CEP-currents with fits ---------------------




