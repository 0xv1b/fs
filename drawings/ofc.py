import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from brokenaxes import brokenaxes
import labellines
mpl.rcParams['axes.linewidth'] = 1.5

# Parameters
offset = 10
f0 = 1

N = 20                 # number of comb lines
f_rep = 1.0             # spacing (arbitrary units)
frequencies = (np.arange(N) * f_rep) + offset
gray_frequencies = np.arange(N) * f_rep

# Gaussian spectral envelope for amplitudes
center = (N/2) + offset
sigma = N/6
amplitudes = np.exp(-0.5 * ((frequencies - center*f_rep)/(sigma*f_rep))**2)

# Normalize amplitudes
amplitudes /= amplitudes.max()

# Colormap (rainbow-like)
colors = cm.jet(np.linspace(0, 1, N))  # from red → green → blue

# Plot
fig = plt.figure(figsize=(10,4))
bax = brokenaxes(xlims=((0, 5), (10, N+offset+1)), hspace=0.05)

for i, (f, amp, c) in enumerate(zip(gray_frequencies, amplitudes, colors)):
    bax.vlines(f+f0, 0, 1, color="gray", linewidth=6, linestyles="dashed")

for i, (f, amp, c) in enumerate(zip(frequencies, amplitudes, colors)):
    bax.vlines(f+f0, 0, amp, color=c, linewidth=6)

# Labels and formatting
bax.set_xlabel("Frequency", fontsize=16)
bax.set_ylabel("Optical Power", fontsize=16)

bax.tick_params(axis='both', which='major', labelsize=14)
bax.tick_params(axis='both', which='minor', labelsize=12)

bax.set_ylim(0, 1.1)
#bax.set_xlim(-1, N+offset)

plt.savefig("../images/OFC.pdf")
plt.show()
