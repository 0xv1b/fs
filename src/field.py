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


power = [0.28,0.47,0.72,0.95,1.28,1.51,1.76,2.06,2.23,2.38,2.47,2.59]
phase = [-21.102432439149688,-20.3886048487173,-19.72613141668551,-17.08560367438448,
         -21.8051336398769,-24.24189232729099,-23.46911634882678,-25.722341723213482,-20.612788776857542,-22.76083290924663,-27.61705659461053,-40.00807325012532]
j0 = [5.070370786753486,21.55871981725074,16.8865962072211,54.23287979811419,
      44.11455982906513,25.72390932432658,21.623560558008403,67.49597940123074,88.89865880517732,114.26508344839799,149.27380568123354,120.75431841136434]

fig, ax = plt.subplots()
ax.tick_params(axis='x', direction='in')
ax.tick_params(axis='y', direction='in')
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xlabel(r'$\varepsilon_{0}$ (V/nm)', fontsize=14)
ax.set_ylabel(r'$\Theta$ ($2\pi$ radians)', fontsize=14)

# def j2theta(x):
#     return  2*x/max(j0)

# def theta2j(x):
#     return x*max(j0)/2

# period_axis = ax.secondary_yaxis('right', functions=(j2theta, theta2j), color="red")
# period_axis.set_ylabel(r'Theta (2$\pi$ radians)', fontsize=12)
# period_axis.tick_params(axis='y', direction='in')
# period_axis.tick_params(axis='both', which='major', labelsize=12)
# period_axis.tick_params(axis='both', which='minor', labelsize=12)
# ax.plot(power, j0, marker=".", linestyle="--", markersize=8, color="black")
# ax.plot(power, [(x % (2*np.pi))*max(j0)/(2*np.pi) for x in phase], marker="x", linestyle="", markeredgewidth=2,markersize=8, color="red")

ax.plot(power, [(x % (2*np.pi))/(2*np.pi) for x in phase], marker="x", linestyle="", markeredgewidth=2,markersize=8, color="red")
ax.set_ylim(0,1)
plt.tight_layout()
plt.savefig("./images/power_phase.pdf")
plt.show()