import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from matplotlib import rcParams

rcParams['axes.linewidth'] = 1.5

import os




# --- Lorentzian function ---
def lorentzian(x, A, x0, gamma):
    return A * gamma**2 / ((x - x0)**2 + gamma**2)

# Load data (adjust path as needed)
data = pd.read_csv(r"data/Vibhor_MoS2/Raman_SampleC_long.txt", sep='\t', comment="'", skiprows=13,
                   names=["Wavenumber", "Intensity"], encoding='latin1')

# --- Clean data ---
data = data[pd.to_numeric(data["Wavenumber"], errors='coerce').notnull()]
data["Wavenumber"] = pd.to_numeric(data["Wavenumber"])
data["Intensity"] = pd.to_numeric(data["Intensity"])

x = data["Wavenumber"].values
y = data["Intensity"].values

# --- Filter region 300–450 cm⁻¹ ---
mask = (x >= 320) & (x <= 450)
x = x[mask]
y = y[mask]

# --- Background subtraction using 5th-degree polynomial ---
background_coeffs = np.polyfit(x, y, deg=1)
background = np.polyval(background_coeffs, x)
y_corrected = y - background

# --- Peak finding ---
peaks, _ = find_peaks(y_corrected, height=50, distance=20, prominence=20)
peak_positions = x[peaks]

# --- Fit Lorentzian to each peak ---
fitted_curves = []
fitted_params = []

for peak in peak_positions:
    mask = (x > peak - 35) & (x < peak + 35)
    x_peak = x[mask]
    y_peak = y_corrected[mask]
    p0 = [max(y_peak), peak, 5]
    try:
        popt, _ = curve_fit(lorentzian, x_peak, y_peak, p0=p0)
        fitted_curves.append((x_peak, lorentzian(x_peak, *popt)))
        fitted_params.append(popt)
    except RuntimeError:
        continue

# --- Plot ---
plt.figure(figsize=(8, 6))
plt.plot(x, y_corrected, label='Background-subtracted', color='black', linewidth=1.5)

# plt.axvspan(383.7, 389.6, color='0.9')
# plt.axvspan(403.1, 410, color='0.9')

colors = plt.cm.viridis(np.linspace(0, 1, len(fitted_curves)))
for (x_fit, y_fit), params, color in zip(fitted_curves, fitted_params, colors):
    plt.plot(x_fit, y_fit, '--', color=color, linewidth=1.2)
    plt.text(params[1], max(y_fit) + 2, f"{params[1]:.1f}", 
             ha='center', va='bottom', fontsize=10, color="gray")
    plt.fill(x_fit, y_fit, 'gray')

plt.xlabel('Raman Shift (cm⁻¹)', fontsize=18)
plt.ylabel('Intensity (a.u.)', fontsize=18)
plt.tight_layout()
plt.ylim(0,300)

plt.tick_params(axis='both', which='major', labelsize=14)
plt.xlim(360,440)
plt.show()



# --- Setup ---
samples = ['A', 'C', 'E', 'G', 'H', 'I']
file_template = r"data/Vibhor_MoS2/Raman_Sample{}_long.txt"

# Plot style for Nature-like figure
rcParams.update({
    'font.size': 12,
    'font.family': 'serif',
    'axes.linewidth': 1.2,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.width': 1.1,
    'ytick.major.width': 1.1,
    'legend.frameon': False
})

plt.figure(figsize=(10, 6))

# --- Loop through each sample ---
for sample in samples:
    file_name = file_template.format(sample)
    
    if not os.path.exists(file_name):
        print(f"File not found: {file_name}")
        continue

    # Load and clean data
    data = pd.read_csv(file_name, sep='\t', comment="'", skiprows=13,
                       names=["Wavenumber", "Intensity"], encoding='latin1')
    data = data[pd.to_numeric(data["Wavenumber"], errors='coerce').notnull()]
    data["Wavenumber"] = pd.to_numeric(data["Wavenumber"])
    data["Intensity"] = pd.to_numeric(data["Intensity"])
    
    x = data["Wavenumber"].values
    y = data["Intensity"].values

    # Filter PlotWindow
    mask = (x >= 370) & (x <= 430)
    x_crop = x[mask]
    y_crop = y[mask]

    # Filter for BKG
    maskBKG = (x >= 370) & (x <= 373)
    bkg = np.mean(y[maskBKG])

    # Background subtraction (5th-order polynomial)
    background = bkg#np.polyval(np.polyfit(x, y, deg=1), x)
    y_corrected = y_crop - background

    # Normalize (optional, for fair comparison)
    y_corrected /= np.max(y_corrected)

    # Plot
    plt.plot(x_crop, y_corrected, label=f"Sample {sample}")

# --- Final touches ---
plt.xlabel('Raman Shift (cm⁻¹)')
plt.ylabel('Counts (a.u.)')
plt.legend()
plt.tight_layout()
plt.show()



