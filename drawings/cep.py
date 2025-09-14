import numpy as np
import matplotlib.pyplot as plt
import labellines
import matplotlib as mpl
from brokenaxes import brokenaxes

mpl.rcParams['axes.linewidth'] = 1.5

# Parameters
f_rep = 1e12        # Repetition rate (Hz)
T_pulse = 1 / f_rep # Pulse period
f_carrier = 5e12    # Carrier frequency (Hz)
omega = 2 * np.pi * f_carrier
tau = 50e-15        # Pulse duration (s)

n_pulses = 4        # Number of pulses to plot
t = np.linspace(-0.5*T_pulse, n_pulses*T_pulse, 5000)

# Function to generate a Gaussian-modulated sinusoidal pulse with CEP
def pulse(t, t0, phi0):
    envelope = np.exp(-((t - t0) ** 2) / (2 * tau ** 2))
    carrier = np.cos(omega * (t - t0) + phi0)
    return envelope * carrier, envelope

# Generate pulse train with CEP shift of pi/2 each pulse
signal = np.zeros_like(t)
envelope_total = np.zeros_like(t)
phi0 = 0
for n in range(n_pulses):
    t0 = n * T_pulse
    p, env = pulse(t, t0, phi0)
    signal += p
    envelope_total += env  # store envelopes separately (non-overlapping in sum)
    phi0 += np.pi / 2  # CEP shift

# Plot

bax = brokenaxes(xlims=((-200, 200), (800, 1200), (1800, 2200), (2800, 3200)), hspace=0.05)
start = 0
end = 4470
bax.plot(t[start:end] * 1e15, signal[start:end], 'k', label="Pulse", color="purple", linewidth=2)  # black field
bax.plot(t * 1e15, envelope_total, '--', alpha=0.8, label="Envelope", color="black")  # dotted Gaussian
bax.plot(t * 1e15, -envelope_total, '--', alpha=0.8, color="black")  # negative dotted Gaussian
# start = 1450
# end = 3100
# bax.plot(t[start:-end] * 1e15, envelope_total[start:-end], '--', alpha=0.8, label="Envelope", color="black")  # dotted Gaussian
# bax.plot(t[start:-end] * 1e15, -envelope_total[start:-end], '--', alpha=0.8, color="black")  # negative dotted Gaussian

# # 1st
# start = 0
# end = 4000
# bax.plot(t[start:end] * 1e15, envelope_total[start:end], '--', alpha=0.8, label="Envelope", color="black")  # dotted Gaussian
# bax.plot(t[start:end] * 1e15, -envelope_total[start:end], '--', alpha=0.8, color="black")  # negative dotted Gaussian


bax.tick_params(axis='both', which='major', labelsize=14)
#bax.tick_params(axis='both', which='minor', labelsize=12)
bax.set_xlabel("Time (fs)", fontsize=16)
bax.set_ylabel(r'$\varepsilon_{0}$ (V/nm)', fontsize=16)
# ax.axvspan(-505, 505, facecolor='red', alpha=0.3)
# ax.axvspan(505, 1515, facecolor='blue', alpha=0.3)
# ax.axvspan(1515, 2525, facecolor='green', alpha=0.3)
# ax.axvspan(2525, 3530, facecolor='purple', alpha=0.3)
#labellines.labelLines(ax.get_lines(),xvals=(100, 1500), fontsize=16)

#plt.savefig("../images/CEP.pdf")
plt.show()
