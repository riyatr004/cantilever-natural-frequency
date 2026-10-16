import numpy as np
from scipy.optimize import brentq

# ---------- 1. Section and material (same as your 1.py) ----------
E = 210e9; D = 6.0; t = 0.07; L = 60.0; rho = 7850
Di = D - 2*t
I = np.pi/64 * (D**4 - Di**4)
A = np.pi/4 * (D**2 - Di**2)
m = rho * A

# ---------- 2. Characteristic equations f(x) = 0, with x = beta*L ----------
boundary_conditions = {
    "Fixed-free (cantilever)": lambda x: np.cos(x)*np.cosh(x) + 1,
    "Pinned-pinned":           lambda x: np.sin(x),
    "Fixed-pinned":            lambda x: np.sin(x)*np.cosh(x) - np.cos(x)*np.sinh(x),
    "Fixed-fixed":             lambda x: np.cos(x)*np.cosh(x) - 1,
}

# ---------- 3. Find the first n roots by scanning for sign changes ----------
def find_roots(func, n_roots=3, x_start=0.5, x_end=15, steps=3000):
    xs = np.linspace(x_start, x_end, steps)
    ys = func(xs)
    roots = []
    for i in range(len(xs) - 1):
        if ys[i] * ys[i+1] < 0:                  # sign change -> root in between
            roots.append(brentq(func, xs[i], xs[i+1]))
            if len(roots) == n_roots:
                break
    return roots

# ---------- 4. Frequency from beta*L ----------
def frequency(beta_L):
    return beta_L**2 / (2*np.pi) * np.sqrt(E*I / (m * L**4))

# ---------- 5. Run and print ----------
for name, func in boundary_conditions.items():
    roots = find_roots(func)
    print(f"\n{name}")
    for k, r in enumerate(roots, start=1):
        print(f"  Mode {k}: beta*L = {r:.4f}   f = {frequency(r):.3f} Hz")
