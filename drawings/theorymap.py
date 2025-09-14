import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import labellines
mpl.rcParams['axes.linewidth'] = 2

# Data

phase_x = [1.145, 1.15, 1.17, 1.175]
phase = [0, np.pi, 0, np.pi]

#x = np.linspace(1.14, 1.18, 1000)
x = np.linspace(1.14, 1.18, 1000)

mu1, sigma1 = 1.15, 0.005

# Bump 2: centered at x=2 with a standard deviation of 1.0
mu2, sigma2 = 1.17, 0.005

# Calculate the y values by summing two Gaussian functions
y1 = 1*np.exp(-0.5 * ((x - mu1) / sigma1)**2)
y2 = 1*np.exp(-0.5 * ((x - mu2) / sigma2)**2)
y = y1 + y2

print(y)
# ---------------------------------------------------------------------------- #
# Plot 1 (for xscan180_keys and xscan180_j0)                                  #
# ---------------------------------------------------------------------------- #
fig1, ax1 = plt.subplots()

ax1.tick_params(axis='x', direction='in')
ax1.tick_params(axis='y', direction='in')
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_xlim(1.14, 1.18) 
ax1.set_ylim(0,1)
ax1.set_xlabel('x position', fontsize=18)
ax1.set_ylabel(r'$j_0$ (normalised)', fontsize=18)

def j2theta(x):
    return  x

def theta2j(x):
    return x

# Create secondary y-axis for theta
period_axis1 = ax1.secondary_yaxis('right', functions=(j2theta, theta2j), color="red")
period_axis1.set_ylabel(r'$\Theta_{CEP,ins}$ ($\pi$ radians)', fontsize=18)
period_axis1.tick_params(axis='y', direction='in')
period_axis1.tick_params(axis='both', which='major', labelsize=12)
period_axis1.tick_params(axis='both', which='minor', labelsize=12)

# Photocurrent Plot
ax1.plot(x, y, markersize=8, color="black", lw=4, label="Photocurrent") #, xerr=0.0025)

ax1.axvspan(1.145, 1.15, facecolor='red', alpha=0.5)
ax1.axvspan(1.17, 1.175, facecolor='blue', alpha=0.5)

# Theta Plot
ax1.plot(phase_x, [(val)/(np.pi) for val in phase], markeredgewidth=2, markersize=8, color="red", lw=4, label="Phase")

labellines.labelLines(ax1.get_lines(), fontsize=16)
plt.xticks([])
plt.tight_layout()
plt.show()
