import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt

# ---------- 1. Inputs ----------
E = 210e9; D = 6.0; t = 0.07; L = 60.0; rho = 7850
Di = D - 2*t
I = np.pi/64 * (D**4 - Di**4)
A = np.pi/4 * (D**2 - Di**2)
m = rho * A

# ---------- 2. Beam element (same as before) ----------
def beam_element(E, I, m, l):
    k = E*I/l**3 * np.array([[ 12,   6*l,   -12,   6*l],
                             [6*l, 4*l**2, -6*l, 2*l**2],
                             [-12,  -6*l,    12,  -6*l],
                             [6*l, 2*l**2, -6*l, 4*l**2]])
    M = m*l/420 * np.array([[156,   22*l,   54,  -13*l],
                            [22*l, 4*l**2, 13*l, -3*l**2],
                            [54,    13*l,  156,  -22*l],
                            [-13*l, -3*l**2, -22*l, 4*l**2]])
    return k, M

# ---------- 3. Solver with tip mass and base springs ----------
def solve_beam(n_el=20, n_modes=3, M_tip=0.0, k_x=1e15, k_r=1e18):
    ndof = 2*(n_el + 1)
    l = L / n_el
    K = np.zeros((ndof, ndof)); M = np.zeros((ndof, ndof))
    ke, me = beam_element(E, I, m, l)
    for e in range(n_el):
        d = [2*e, 2*e+1, 2*e+2, 2*e+3]
        K[np.ix_(d, d)] += ke
        M[np.ix_(d, d)] += me
    K[0, 0] += k_x          # lateral soil spring at the base (N/m)
    K[1, 1] += k_r          # rotational soil spring at the base (Nm/rad)
    M[-2, -2] += M_tip      # nacelle + rotor mass at the top node (kg)
    lam = eigh(K, M, eigvals_only=True)
    return np.sqrt(lam[:n_modes]) / (2*np.pi)

# ---------- 4. Experiment A: tip mass ----------
print("Fixed base, no tip mass:", np.round(solve_beam(), 3))
for Mt in [0, 200e3, 400e3, 800e3]:
    print(f"Tip mass {Mt/1e3:5.0f} t -> f1 = {solve_beam(M_tip=Mt)[0]:.3f} Hz")

# ---------- 5. Experiment B: soil stiffness ----------
kr_values = np.logspace(9, 13, 30)                  # rotational stiffness range
f1_soil = [solve_beam(M_tip=400e3, k_r=kr)[0] for kr in kr_values]

plt.semilogx(kr_values, f1_soil, "o-")
plt.xlabel("Base rotational stiffness k_r (Nm/rad)")
plt.ylabel("First natural frequency (Hz)")
plt.title("Effect of soil stiffness (tip mass 400 t)")
plt.grid(True, which="both")
plt.show()
