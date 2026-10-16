import numpy as np

# Example steel tube cantilever (made-up values, just for practice)
E = 210e9          # Pa, steel stiffness
D = 6.0            # m, outer diameter
t = 0.07           # m, wall thickness
L = 60.0           # m, length
rho = 7850         # kg/m3, steel density

Di = D - 2*t
I = np.pi/64 * (D**4 - Di**4)       # second moment of area
A = np.pi/4 * (D**2 - Di**2)        # cross-section area
m = rho * A                          # mass per metre

f1 = (1.875**2) / (2*np.pi) * np.sqrt(E*I / (m * L**4))
print(f"First natural frequency: {f1:.3f} Hz")
