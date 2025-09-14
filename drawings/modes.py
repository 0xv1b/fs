import numpy as np
import matplotlib.pyplot as plt

# Time axis
t = np.linspace(-2*np.pi, 2*np.pi, 4000)

# Parameters
N_modes = 15              # number of frequency modes
omega_rep = 2*np.pi / 2   # mode spacing (repetition rate)
phi0 = 0                  # same phase for mode locking

# Build electric field: sum of modes
field = np.zeros_like(t)
for n in range(-(N_modes//2), N_modes//2 + 1):
    field += np.cos(n * omega_rep * t + phi0)

# Normalize
field /= np.max(np.abs(field))

# Plot
fig, axs = plt.subplots(2, 1, figsize=(10,6), sharex=True)

# --- Top: a few representative modes (like sinusoids) ---
for n in range(1, 6):
    axs[0].plot(t, np.cos(n * omega_rep * t + phi0), label=f"mode {n}")
axs[0].set_ylabel("Amplitude")
axs[0].set_title("Individual Modes (locked in phase)")
axs[0].legend(loc="upper right")

# --- Bottom: sum = pulse train ---
axs[1].plot(t, field, 'k', linewidth=2)
axs[1].set_xlabel("Time (a.u.)")
axs[1].set_ylabel("Electric field")
axs[1].set_title("Mode-Locked Pulse Train")

plt.savefig("../images/Mode_Locking.pdf")
plt.tight_layout()
plt.show()
