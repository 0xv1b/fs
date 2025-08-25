import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from scipy.signal import find_peaks
from .helpers import common, constants


for j, filename in enumerate(constants.OSC):


    df = pd.read_csv(filename, header=2, usecols=[3,4])
    df.columns = ['time', 'current']
    fig, ax = plt.subplots()
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.tick_params(axis='x', direction='in')
    ax.tick_params(axis='y', direction='in')

    up_zeroes = []
    down_zeroes = []

    for i, val in enumerate(np.diff(np.sign(df.current))):
        if val == 2 or val == -2 or df.current[i] == 0:
            if df.current[i+3] > 0:
                #up_zeroes[df.time[i]] = i
                #ax.axvline(df.time[i], color='darkgrey', linestyle='--')
                #if(np.abs(up_zeroes[0] - df.time[i]) < 0.0005):
                    #continue
                up_zeroes.append(df.time[i])
            if df.current[i+3] < 0:
                down_zeroes.append(df.time[i])
                #ax.axvline(df.time[i], color='red', linestyle='--')
    


    if j == 0:
        print("hardcoded fix")
        up_zeroes.pop(0)

    

    for i in range(0, len(up_zeroes),1):
        # Edge cases are if it starts with a downzero or if it ends in an upzero
        # Case 1: starts with downzero ends in downzero : downzero array is 1 bigger than upzero
        # Case 2: starts with upzero ends in upzero : upzero array is 1 bigger than downzero
        # Case 3: starts with downzero and ends in upzero: same size arrays.
        # Knowing the data they all start with down zeroes
        ax.axvspan(up_zeroes[i], down_zeroes[i+1], color='0.9')

    peaks, _ = find_peaks(df.current)
    valley_peaks = [peak for peak in peaks if df.current[peak] < 0]
    valley_data = [df.current[i] for i in valley_peaks]
    valley_offset = np.mean(valley_data)    
    

    plt.text((up_zeroes[0] + down_zeroes[1])/2, 1.1 - valley_offset, 'On', horizontalalignment='center', fontsize=14)
    plt.text((up_zeroes[1] + down_zeroes[1])/2, 1.1 - valley_offset, 'Off', horizontalalignment='center', fontsize=14)
    plt.axhline(0, color='darkgrey', linestyle='--')


    plt.xlabel('Time (ms)', fontsize='18')
    plt.ylabel('Photocurrent (pA)', fontsize='18')
    #plt.title(str(filename), fontsize='16')
    plt.plot(df.time,df.current - valley_offset)
    #plt.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout()
    plt.show()

    #plt.close(fig)
