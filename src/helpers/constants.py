import h5py
import glob
import numpy as np
from scipy.interpolate import CubicSpline

POWER_CAL_FILE = h5py.File(r'data/powercl.h5','r')
angle_powerCal = POWER_CAL_FILE['measurement']['Zman_script']['pump_powerwheel:position:sweep'][0]
power_powerCal = POWER_CAL_FILE['measurement']['Zman_script']['powermeter:power'][0]

power_powerCal_interpol = CubicSpline(angle_powerCal[50:600],power_powerCal[50:600]*1e3)
angle_powerCal_interpol = np.linspace(30,200)


MAP_FILE = h5py.File(r'data/Maps/20250610_31_5V-SampleI_2Dmap_327Hz_180deg_BF_0.0V_wedge6mm.h5', 'r')
current_magnitude = MAP_FILE['measurement']['Zman_script']['LockIn_3:MagPhase1'][0]
current_phase = MAP_FILE['measurement']['Zman_script']['LockIn_3:MagPhase1'][1]
refl_map = MAP_FILE['measurement']['Zman_script']['logitech_webcam:intensity'][0]

CEP_POWER = sorted(glob.glob('data/CEP_power/*'))
CEP_LONG = sorted(glob.glob('data/CEP_long/*'))
CEP_XSCAN_180 = sorted(glob.glob('data/CEP_xscan_180/*'))
CEP_XSCAN = sorted(glob.glob('data/CEP_xscan/*'))

CEP_RIGHT_WO =  sorted(glob.glob('data/CEP_right_wo/*'))
CEP_RIGHT_GREEN =  sorted(glob.glob('data/CEP_right_green/*'))
CEP_MIDDLE_WO = sorted(glob.glob('data/CEP_middle_wo/*'))
CEP_MIDDLE_GREEN =  sorted(glob.glob('data/CEP_middle_green/*'))

OSC =sorted(glob.glob('data/Osc/*.CSV'))

