from . import common
import numpy as np
import re
import h5py
import matplotlib.pyplot as plt

def init_lists(filenames, scans, scans_detrend, wedge_pos_array, projection_phase=116, xscan=False):
    for filename in filenames:
        file = h5py.File(str(filename),'r')
        # Extract powerwheel angle from filename
        pos_filename_start = []
        pos_filename_end = []

        if xscan:
            # Workaround for xscan
            for alles in re.finditer('xpos_',filename):
                pos_filename_start.append(alles.end())
            for alles in re.finditer('_BF',filename):
                pos_filename_end.append(alles.start())

        else:
            for alles in re.finditer('left_',filename):
                pos_filename_start.append(alles.end())
            for alles in re.finditer('deg_',filename):
                pos_filename_end.append(alles.start())
        
        pos = float(filename[pos_filename_start[0]:pos_filename_end[0]])

        if pos in scans.keys():
            pass
        else:
            scans[pos] = []# Adds pos key and assigns list
            scans_detrend[pos] = []# Adds pos key and assigns list

        wedge_pos_array.append(file['measurement']['Zman_script']['wedge_stage:position:sweep'][0]) # Adds list to list
        current = common.func_current(file['measurement']['Zman_script']['LockIn_3:MagPhase1'][0],file['measurement']['Zman_script']['LockIn_3:MagPhase1'][1],projection_phase)
        # Calculate current from LockIn and phase data
        scans[pos].append(current) # Adds currents list to list

def get_wedge_positions(wedge_pos_array, decimals=1):
    return np.round(np.mean(wedge_pos_array,axis=0),decimals)

def get_detrended_currents(scans, currents, deg=8):
    currents_detrend = {}
    for key in scans.keys():
        currents[key] = np.mean(scans[key],axis=0)
    
    for key in np.flip(np.sort([*currents])):
        x = np.arange(len(currents[key]))
        coefficients = np.polyfit(x, currents[key], deg)
        polyfit = np.polyval(coefficients, x)
        currents_detrend[key] = currents[key] - polyfit

        #fig, ax = plt.subplots()
        #ax.plot(scans[key], currents[key])
        #plt.show()
    


    return currents_detrend

def init_pump(filenames, wedge_pos_array, trace,projection_phase=116):
    for filename in filenames:
        file = h5py.File(str(filename),'r')
        
        wedge_pos_array.append(file['measurement']['Zman_script']['wedge_stage:position:sweep'][0])
        
        current = common.func_current(file['measurement']['Zman_script']['LockIn_3:MagPhase1'][0],
                            file['measurement']['Zman_script']['LockIn_3:MagPhase1'][1],
                            projection_phase)
        
        trace.append(current)

def get_detrended_currents_pump(currents, trace, deg=8, positions=[]):
    currents_detrend = []
    currents = np.mean(trace,axis=0)
    #currents_detrend = polynomial(currents,order=8,plot=True)

    x = np.arange(len(currents))
    coefficients = np.polyfit(x, currents, deg)
    polyfit = np.polyval(coefficients, x)
    currents_detrend = currents - polyfit

    return currents_detrend

def pump_sanity_check(wedge_pos_array, expected_len=102):
    for i, arr in enumerate(wedge_pos_array):
        if len(arr) != expected_len:
            print(f"⚠️ Entry {i} has length {len(arr)}, expected {expected_len}")