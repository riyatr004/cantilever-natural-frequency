import numpy as np
import matplotlib.pyplot as plt

E = 210e9; rho = 7850

# Reference values (the ones NOT being varied stay fixed at these)
D0 = 6.0; t0 = 0.07; L0 = 60.0

def first_frequency(D, t, L):
    Di = D - 2*t
    I = np.pi/64 * (D**4 - Di**4)
    A = np.pi/4 * (D**2 - Di**2)
    m = rho * A
    return (1.875**2)/(2*np.pi) * np.sqrt(E*I/(m*L**4))

L_values = np.linspace(40, 140, 30)
D_values = np.linspace(3, 12, 30)
t_values = np.linspace(0.03, 0.15, 30)

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# 1) f vs L
axes[0].plot(L_values, first_frequency(D0, t0, L_values), "o-", ms=3)
axes[0].set_xlabel("Length L (m)")
axes[0].set_title("f vs L  (D = 6 m, t = 0.07 m)")

# 2) f vs D
axes[1].plot(D_values, first_frequency(D_values, t0, L0), "o-", ms=3)
axes[1].set_xlabel("Outer diameter D (m)")
axes[1].set_title("f vs D  (L = 60 m, t = 0.07 m)")

# 3) f vs t
axes[2].plot(t_values, first_frequency(D0, t_values, L0), "o-", ms=3)
axes[2].set_xlabel("Wall thickness t (m)")
axes[2].set_title("f vs t  (L = 60 m, D = 6 m)")

for ax in axes:
    ax.set_ylabel("First natural frequency (Hz)")
    ax.grid(True)

plt.tight_layout()
plt.savefig("L_D_t_vs_f.png", dpi=200)   # saves the figure for your report
plt.show()
