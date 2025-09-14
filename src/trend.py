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

print(round(common.power2Efield(float(constants.power_powerCal_interpol(130))*1e-3,6.5e-15,5e-6)*1e-9,2))
print(round(common.power2Efield(float(constants.power_powerCal_interpol(180))*1e-3,6.5e-15,5e-6)*1e-9,2))