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


xscan_keys = [1.148,1.1529999999999998,1.1579999999999997,1.1629999999999996,1.1679999999999995,1.1729999999999994,1.1779999999999993]
xscan_ipl = [0.11574724965210058,0.12552673884479312,0.1138344625198837,0.14332755563773583,0.10149822657938015,0.11918298218851806,0.11972916302992587] # Changed negative value to positive
xscan_theta = [-23.008359106593396,-20.062554042460626,3.5189999063273367,-17.817791094417586,-27.179657761455807,-22.126972672615643,-22.166268016571046]
xscan_j0 = [201.05294242176504,30.292548463526956,122.70158036062176,230.51397993483187,214.04398823427226,231.93641657727756,88.89718404682124]
xscan_j0 = [x/10 for x in xscan_j0]


xscan180_keys = [1.14,1.142,1.144,1.146,1.148,1.15,1.152,1.155,1.16,1.1649999999999998,1.1699999999999997,1.17,1.1749999999999998]
xscan180_ipl = [0.11379581923671946,0.11524709865669037,0.12064651975203011,0.11377574612685347,0.11483895658606236,0.1160681600926213,0.11647405989454629,
                0.11346541617336277,0.11348984443762429,0.11475829886946332,0.05073649591698381,0.11525522614970271,0.11833683309186913]
# do modulo 2 pi
xscan180_theta = [-23.221106044677356,-29.02769327091837,-22.113574695256105,-23.616499551540038,-23.775024784657226,-23.85319172915859,19.610521391089875,
                  -22.028408121269443,-23.844956269473137,-24.356662353150725,-29.07091844282327,-57.417648186603415,-19.00252995004415] 
xscan180_j0 = [12.622341813010472,21.01546452044434,36.34438146293278,16.341659214909868,4.177746321105634,4.705219111838885,17.526132323491964,
               17.993766591456392,21.811901793716814,505.51998784435705,137.77937716942273,171.56244274060458,10.69313163264829]
xscan180_j0 = [x/10 for x in xscan180_j0]

# ---------------------------------------------------------------------------- #
#                                   Plotting                                   #
# ---------------------------------------------------------------------------- #
fig, axes = plt.subplots(2,1)
for ax in axes:
    ax.tick_params(axis='x', direction='in')
    ax.tick_params(axis='y', direction='in')
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlim(1.14,1.18) 
    ax.set_xlabel('X position (mm)', fontsize=14)
    ax.set_ylabel(r'$j_0$ (pA)', fontsize=14)

def j2theta(x):
    return  2*x/max(xscan180_j0)

def theta2j(x):
    return x*max(xscan180_j0)/2

period_axis = axes[0].secondary_yaxis('right', functions=(j2theta, theta2j), color="red")
period_axis.set_ylabel(r'Theta (2$\pi$ radians)', fontsize=12)
period_axis.tick_params(axis='y', direction='in')
period_axis.tick_params(axis='both', which='major', labelsize=12)
period_axis.tick_params(axis='both', which='minor', labelsize=12)

def j2theta(x):
    return  2*x/max(xscan_j0)

def theta2j(x):
    return x*max(xscan_j0)/2

deriod_axis = axes[1].secondary_yaxis('right', functions=(j2theta, theta2j), color="red")
deriod_axis.set_ylabel(r'Theta (2$\pi$ radians)', fontsize=12)
deriod_axis.tick_params(axis='y', direction='in')
deriod_axis.tick_params(axis='both', which='major', labelsize=12)
deriod_axis.tick_params(axis='both', which='minor', labelsize=12)
# ---------------------------------------------------------------------------- #
#                                 Photocurrent                                 #
# ---------------------------------------------------------------------------- #
axes[1].errorbar(xscan_keys, xscan_j0, marker=".", linestyle="--", markersize=8, color="black",xerr=0.0025)
axes[1].axvspan(1.16752, 1.17649, facecolor='blue', alpha=0.5)
axes[1].axvspan(1.14451, 1.15051, facecolor='red', alpha=0.5)

axes[0].errorbar(xscan180_keys, xscan180_j0, marker=".", linestyle="--", markersize=8, color="black",xerr=0.0025)
axes[0].axvspan(1.16952, 1.17547, facecolor='blue', alpha=0.5)
axes[0].axvspan(1.14646, 1.14949, facecolor='red', alpha=0.5)


# ---------------------------------------------------------------------------- #
#                                     Theta                                    #
# ---------------------------------------------------------------------------- #
axes[1].plot(xscan_keys, [(x % (2*np.pi))*max(xscan_j0)/(2*np.pi) for x in xscan_theta], marker="x", linestyle="", markeredgewidth=2,markersize=8, color="red")
axes[0].plot(xscan180_keys, [(x % (2*np.pi))*max(xscan180_j0)/(2*np.pi) for x in xscan180_theta], marker="x", linestyle="", markeredgewidth=2,markersize=8, color="red")


plt.tight_layout()
plt.show()