import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from matplotlib import rcParams

import os

rcParams['axes.linewidth'] = 1.5

# --- Lorentzian function ---
def lorentzian(x, A, x0, gamma):
    return A * gamma**2 / ((x - x0)**2 + gamma**2)

# List of sample letters
samples = ['A', 'C', 'E', 'G', 'H', 'I']
#samples = ['A', 'B', 'C', 'D', 'E', 'I']
file_template = r"data/Vibhor_MoS2/PL_Sample{}.txt"


plt.figure(figsize=(8, 6))

# Process each PL file
for sample in samples:
    file_path = file_template.format(sample)
    if not os.path.exists(file_path):
        print(f"Missing: {file_path}")
        continue

    # Manually extract data from [Data] section
    with open(file_path, 'r', encoding='latin1') as f:
        lines = f.readlines()
    data_start = next(i for i, line in enumerate(lines) if line.strip() == '[Data]')
    data_lines = lines[data_start + 3:]
    data = [line.strip().split() for line in data_lines if len(line.strip().split()) == 2]

    # Convert to DataFrame
    df = pd.DataFrame(data, columns=['Energy', 'Intensity']).astype(float)
    x = df['Energy'].values
    y = df['Intensity'].values

    # Background subtraction via 5th-degree polynomial
    background = y[-1]
    y_corrected = y - background
    #y_corrected /= np.max(y_corrected)  # normalize

    # Plot background-subtracted spectrum
    plt.plot(x, y_corrected) #, label=f'Sample {sample}')

# Final plot formatting
plt.xlabel('Photon Energy (eV)', fontsize=18)
plt.xlim(1.70,2.10)
plt.ylabel('Counts (a.u.)', fontsize=18)
plt.legend(["Sample A","Sample B","Sample C","Sample D","Sample E","Sample F"], fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.tight_layout()
plt.show()