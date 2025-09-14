import numpy as np
import matplotlib.pyplot as  plt
import h5py
import glob
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
import re
import matplotlib as mpl
from .helpers import common, detrend, constants, fits
mpl.rcParams['axes.linewidth'] = 1.5

filename200 = r"data/0605/20250605_13_5V-SampleI_I-V_curve_BF_0.3_V_underill200deg_right.h5"
h5_200= h5py.File(filename200, 'r')
voltage200 = h5_200['measurement']['Zman_script']['yokogawa:V:sweep'][:][0]
current200 = h5_200['measurement']['Zman_script']['multimeter_3:voltage'][:][:][0]

filename220 = r"data/0605/20250605_14_5V-SampleI_I-V_curve_BF_0.3_V_underill220deg_right.h5"
h5_220= h5py.File(filename220, 'r')
voltage220 = h5_200['measurement']['Zman_script']['yokogawa:V:sweep'][:][0]
current220 = h5_220['measurement']['Zman_script']['multimeter_3:voltage'][:][:][0]

filename240 = r"data/0605/20250605_16_5V-SampleI_I-V_curve_BF_0.3_V_underill240deg_right.h5"
h5_240= h5py.File(filename240, 'r')
voltage240 = h5_240['measurement']['Zman_script']['yokogawa:V:sweep'][:][0]
current240 = h5_240['measurement']['Zman_script']['multimeter_3:voltage'][:][:][0]

filename190 = r"data/0605/20250605_17_5V-SampleI_I-V_curve_BF_0.3_V_underill190deg_right.h5"
h5_190= h5py.File(filename190, 'r')
voltage190 = h5_190['measurement']['Zman_script']['yokogawa:V:sweep'][:][0]
current190 = h5_190['measurement']['Zman_script']['multimeter_3:voltage'][:][:][0]

filename180 = r"data/0605/20250605_15_5V-SampleI_I-V_curve_BF_0.3_V_underill180deg_right.h5"
h5_180= h5py.File(filename180, 'r')
voltage180 = h5_180['measurement']['Zman_script']['yokogawa:V:sweep'][:][0]
current180 = h5_180['measurement']['Zman_script']['multimeter_3:voltage'][:][:][0]

filename0 = r"data/0605/20250605_18_5V-SampleI_I-V_curve_BF_0.3_V_dark_right.h5"
h5_0= h5py.File(filename0, 'r')
voltage0 = h5_0['measurement']['Zman_script']['yokogawa:V:sweep'][:][0]
current0 = h5_0['measurement']['Zman_script']['multimeter_3:voltage'][:][:][0]

# Plot I-V curve
#plt.figure(figsize=(8, 6))


fig, ax = plt.subplots()
ax.tick_params(axis='x', direction='in')
ax.tick_params(axis='y', direction='in')


beam_diameter = 5e-6
pulse_FWHM = 6.5e-15

start = 30
end = 92
plt.plot(voltage0[start:end], current0[start:end]*1e-6*1e6, 'k-', linewidth=1.5, label='Dark') # Black
plt.plot(voltage240[start:end], current240[start:end]*1e-6*1e6, color='#8B0000', linewidth=1.5, label=str(round(common.power2Efield(float(constants.power_powerCal_interpol(240))*1e-3,6.5e-15,5e-6)*1e-9,2))  + ' V/nm')  # Dark Red
# plt.plot(voltage220, current220, color='#B22222', linewidth=1, label='220')  # Firebrick
# plt.plot(voltage200, current200, color='#FF4500', linewidth=1, label='200')  # OrangeRed
# plt.plot(voltage190, current190, color='#FF6347', linewidth=1, label='190')  # Tomato
# plt.plot(voltage180, current180, color='#FFA500', linewidth=1, label='180')  # Orange

plt.legend(fontsize=14)
plt.xlabel('Voltage (V)',fontsize=18)
plt.ylabel(r'Current ($\mu A$)',fontsize=18) # Why arbitrary units?
plt.tick_params(axis='both', which='major', labelsize=14)
plt.savefig("./images/IV.pdf")
plt.show()