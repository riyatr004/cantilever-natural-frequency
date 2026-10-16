"""
Phase A (continued): realistic reference monopile on distributed soil springs.

Euler-Bernoulli beam FE model (2 DOF per node: lateral deflection w, rotation theta)
- varying section along height (pile, transition piece, tapered tower)
- distributed Winkler soil springs  k(z) = n_h * depth   (linear with depth)
- optional hydrodynamic added mass on the submerged part
- four installation stages: bare pile, +TP, +tower, +nacelle/rotor
- all values ILLUSTRATIVE (generic ~8 MW class), not site data

Coordinate z [m]: 0 = mudline, +z up, -z into the seabed.
Run:  python phase_a_monopile.py      (writes figures + CSV to ./results)
"""
import os
import numpy as np
import scipy.linalg as sla
import matplotlib.pyplot as plt

OUT = "results"
os.makedirs(OUT, exist_ok=True)

# ----------------------------------------------------------------------------
# 1. PARAMETER TABLE  (edit here, document every change)
# ----------------------------------------------------------------------------
P = dict(
    E=210e9,            # Pa, Young's modulus of steel
    rho_s=7850.0,       # kg/m3, steel
    rho_w=1025.0,       # kg/m3, sea water
    Ca=1.0,             # added-mass coefficient (circular cylinder)
    # --- vertical levels (z = 0 at mudline) ---
    z_toe=-30.0,        # m, pile toe  -> embedded length 30 m
    z_msl=30.0,         # m, mean sea level -> water depth 30 m
    z_pile_top=40.0,    # m, pile head (10 m above MSL)
    z_tp_top=55.0,      # m, top of transition piece
    z_hub=140.0,        # m, hub (110 m above MSL)
    # --- sections ---
    D_pile=8.0, t_pile=0.090,
    D_tp=8.4,   t_tp=0.080,
    D_tow_base=6.5, D_tow_top=4.5,
    t_tow_base=0.045, t_tow_top=0.022,
    # --- lumped masses ---
    m_platform=200e3,   # kg, TP platform, boat landing, secondary steel (at TP top)
    m_rna=450e3,        # kg, rotor-nacelle assembly (at hub)
    dz_max=1.0,         # m, maximum element length
)

# soil: lateral subgrade modulus gradient n_h  [N/m3]; spring per length = n_h*depth [N/m/m]
NH = {"soft": 5e6, "medium": 15e6, "stiff": 40e6}

STAGES = {1: "Bare monopile", 2: "+ Transition piece", 3: "+ Tower", 4: "+ Nacelle & rotor"}


# ----------------------------------------------------------------------------
# 2. FE BUILDING BLOCKS
# ----------------------------------------------------------------------------
def beam_matrices(EI, mu, l):
    """Euler-Bernoulli element: stiffness and consistent mass (mu = mass per metre)."""
    k = EI / l**3 * np.array([[12, 6*l, -12, 6*l],
                              [6*l, 4*l**2, -6*l, 2*l**2],
                              [-12, -6*l, 12, -6*l],
                              [6*l, 2*l**2, -6*l, 4*l**2]])
    m = mu * l / 420 * np.array([[156, 22*l, 54, -13*l],
                                 [22*l, 4*l**2, 13*l, -3*l**2],
                                 [54, 13*l, 156, -22*l],
                                 [-13*l, -3*l**2, -22*l, 4*l**2]])
    return k, m


def assemble(z, D, t, E, rho_s, madd=None):
    """Assemble global K, M for nodes z and per-element D, t (and added mass/m)."""
    n = len(z)
    K = np.zeros((2*n, 2*n)); M = np.zeros((2*n, 2*n))
    for e in range(n - 1):
        l = z[e+1] - z[e]
        A = np.pi * (D[e]**2 - (D[e] - 2*t[e])**2) / 4
        I = np.pi * (D[e]**4 - (D[e] - 2*t[e])**4) / 64
        mu = rho_s * A + (0.0 if madd is None else madd[e])
        ke, me = beam_matrices(E * I, mu, l)
        idx = slice(2*e, 2*e + 4)
        K[idx, idx] += ke; M[idx, idx] += me
    return K, M


def section(zm, p):
    """Outer diameter and wall thickness at height zm."""
    if zm < p["z_pile_top"]:
        return p["D_pile"], p["t_pile"]
    if zm < p["z_tp_top"]:
        return p["D_tp"], p["t_tp"]
    s = (zm - p["z_tp_top"]) / (p["z_hub"] - p["z_tp_top"])       # linear taper
    return (p["D_tow_base"] + s*(p["D_tow_top"] - p["D_tow_base"]),
            p["t_tow_base"] + s*(p["t_tow_top"] - p["t_tow_base"]))


def build_mesh(stage, p):
    z_top = {1: p["z_pile_top"], 2: p["z_tp_top"], 3: p["z_hub"], 4: p["z_hub"]}[stage]
    bps = [b for b in [p["z_toe"], 0.0, p["z_msl"], p["z_pile_top"], p["z_tp_top"], p["z_hub"]]
           if b <= z_top]
    z = [bps[0]]
    for a, b in zip(bps[:-1], bps[1:]):
        n = int(np.ceil((b - a) / p["dz_max"]))
        z.extend(np.linspace(a, b, n + 1)[1:])
    return np.array(z)


def build_model(stage, nh, p=P, added_mass=True, fixed_mudline=False):
    """Return z, K, M (with soil springs and lumped masses) and the free-DOF mask."""
    z = build_mesh(stage, p)
    zm = 0.5 * (z[:-1] + z[1:])
    DT = np.array([section(x, p) for x in zm]); D, t = DT[:, 0], DT[:, 1]
    madd = None
    if added_mass:                      # submerged part: mudline .. MSL
        sub = (zm > 0.0) & (zm < p["z_msl"])
        madd = np.where(sub, p["Ca"] * p["rho_w"] * np.pi * D**2 / 4, 0.0)
    K, M = assemble(z, D, t, p["E"], p["rho_s"], madd)

    # distributed Winkler springs on embedded nodes (tributary length x n_h x depth)
    n = len(z)
    for i in range(n):
        if z[i] < 0.0:
            lt = 0.5 * ((z[i] - z[i-1]) if i > 0 else 0.0) + 0.5 * ((z[i+1] - z[i]) if i < n-1 else 0.0)
            K[2*i, 2*i] += nh * (-z[i]) * lt

    # lumped masses
    if stage >= 2:
        M[2*np.argmin(abs(z - p["z_tp_top"])), 2*np.argmin(abs(z - p["z_tp_top"]))] += p["m_platform"]
    if stage == 4:
        M[2*(n-1), 2*(n-1)] += p["m_rna"]

    free = np.ones(2*n, bool)
    if fixed_mudline:                   # clamp everything at/below mudline (upper bound)
        for i in range(n):
            if z[i] <= 0.0:
                free[2*i:2*i+2] = False
    return z, K, M, free


def modes(stage, nh, nmodes=3, **kw):
    z, K, M, free = build_model(stage, nh, **kw)
    w2, V = sla.eigh(K[np.ix_(free, free)], M[np.ix_(free, free)])
    f = np.sqrt(np.abs(w2[:nmodes])) / (2*np.pi)
    full = np.zeros((len(free), nmodes)); full[free] = V[:, :nmodes]
    W = full[0::2, :]                                    # deflection part only
    for j in range(nmodes):                              # normalise: max |w| = 1, tip positive
        W[:, j] /= np.max(np.abs(W[:, j]))
        if W[-1, j] < 0: W[:, j] *= -1
    return f, z, W


# ----------------------------------------------------------------------------
# 3. VALIDATION (must still reproduce the Phase A cantilever result)
# ----------------------------------------------------------------------------
def validate_cantilever():
    L, D, t, ne = 60.0, 6.0, 0.07, 20
    z = np.linspace(0, L, ne + 1)
    K, M = assemble(z, np.full(ne, D), np.full(ne, t), P["E"], P["rho_s"])
    w2 = sla.eigh(K[2:, 2:], M[2:, 2:], eigvals_only=True)[:3]
    A = np.pi*(D**2-(D-2*t)**2)/4; I = np.pi*(D**4-(D-2*t)**4)/64
    ana = [(bl**2)/(2*np.pi)*np.sqrt(P["E"]*I/(P["rho_s"]*A*L**4)) for bl in (1.8751, 4.6941, 7.8548)]
    print("Validation (cantilever, 20 el.):  FE =", np.round(np.sqrt(w2)/(2*np.pi), 3),
          " analytical =", np.round(ana, 3))


# ----------------------------------------------------------------------------
# 4. RESULTS
# ----------------------------------------------------------------------------
def main():
    validate_cantilever()

    # mass / geometry summary of the full structure (hand-check input)
    z, *_ = build_model(4, NH["medium"])
    zm = 0.5*(z[:-1] + z[1:]); mass = 0.0
    for a, b, x in zip(z[:-1], z[1:], zm):
        D, t = section(x, P); mass += P["rho_s"]*np.pi*(D**2-(D-2*t)**2)/4*(b-a)
    print(f"\nSteel mass of structure: {mass/1e3:.0f} t  (+ {P['m_platform']/1e3:.0f} t platform, "
          f"+ {P['m_rna']/1e3:.0f} t RNA)")

    # ---- frequency table: stages x soil x (added mass on/off) ----
    rows = []
    print("\nFirst natural frequency f1 [Hz]  (with added mass | without)")
    print(f"{'Stage':22s}" + "".join(f"{k:>16s}" for k in NH) + f"{'fixed@mudline':>16s}")
    for s, name in STAGES.items():
        line = f"{name:22s}"; r = [name]
        for key, nh in NH.items():
            f_w = modes(s, nh)[0][0]; f_n = modes(s, nh, added_mass=False)[0][0]
            line += f"{f_w:8.3f} |{f_n:6.3f}"; r += [f_w, f_n]
        f_fix = modes(s, 1.0, fixed_mudline=True)[0][0]
        line += f"{f_fix:16.3f}"; r.append(f_fix); rows.append(r)
        print(line)
    hdr = "stage," + ",".join(f"f1_{k}_addedmass,f1_{k}_noaddedmass" for k in NH) + ",f1_fixed_mudline"
    np.savetxt(f"{OUT}/frequency_table.csv", [r[1:] for r in rows], delimiter=",", header=hdr, comments="")

    print("\nFirst three modes, medium soil, with added mass [Hz]")
    for s, name in STAGES.items():
        print(f"  {name:22s}", np.round(modes(s, NH['medium'])[0], 3))

    # ---- soil spring stiffness table ----
    print("\nSoil spring table (k = n_h * depth):")
    for key, nh in NH.items():
        print(f"  {key:7s} n_h = {nh/1e6:5.1f} MN/m3 | k at 15 m depth = {nh*15/1e6:6.0f} MN/m/m"
              f" | at 30 m = {nh*30/1e6:6.0f} MN/m/m")

    # ---- Figure: mode shapes per stage (medium soil) ----
    fig, axs = plt.subplots(1, 4, figsize=(15, 6), sharey=False)
    for ax, (s, name) in zip(axs, STAGES.items()):
        f, z, W = modes(s, NH["medium"], 3)
        for j in range(3):
            ax.plot(W[:, j], z, label=f"mode {j+1}: {f[j]:.2f} Hz")
        ax.axhline(0, color="saddlebrown", ls="--", lw=1); ax.axhline(P["z_msl"], color="tab:blue", ls=":", lw=1)
        ax.text(-1.0, 1.5, "mudline", fontsize=8, color="saddlebrown")
        ax.text(-1.0, P["z_msl"]+1.5, "MSL", fontsize=8, color="tab:blue")
        ax.set_title(name); ax.set_xlabel("normalised deflection [-]"); ax.grid(alpha=.3)
        ax.legend(fontsize=8, loc="lower right")
    axs[0].set_ylabel("height z above mudline [m]")
    plt.tight_layout(); plt.savefig(f"{OUT}/mode_shapes_stages.png", dpi=200); plt.close()

    # ---- Figure: f1 versus soil stiffness gradient, per stage ----
    nhs = np.logspace(5.5, 9, 25)
    plt.figure(figsize=(6.5, 4.5))
    for s, name in STAGES.items():
        plt.semilogx(nhs/1e6, [modes(s, n, 1)[0][0] for n in nhs], label=name)
    for key, nh in NH.items():
        plt.axvline(nh/1e6, color="gray", ls=":", lw=1)
    plt.xlabel("subgrade modulus gradient n_h [MN/m$^3$]"); plt.ylabel("first natural frequency [Hz]")
    plt.title("Effect of soil stiffness (distributed springs)"); plt.grid(alpha=.3, which="both")
    plt.legend(fontsize=8); plt.tight_layout(); plt.savefig(f"{OUT}/f1_vs_soil_stiffness.png", dpi=200)
    plt.close()
    print(f"\nSaved figures and CSV to ./{OUT}/")


if __name__ == "__main__":
    main()
