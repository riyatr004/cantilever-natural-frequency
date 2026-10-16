# Fatigue Risk of Bare Offshore Monopiles During Installation: A Python-Based Study

Self-study portfolio project | Riya Thripakkal Rajesh, M.Sc. ACCESS, TU Dresden

## Why this project

Offshore wind foundations (monopiles) are often installed months before the turbine is
connected and monitored. In this period the structure is bare or only partly built. It has a
different mass, stiffness and damping than in its final state and is loaded by waves, current
and wind. Vibrations such as vortex-induced vibration (VIV) can cause fatigue damage that
nobody measures.

**Research question:** how much vibration and fatigue damage can a bare monopile accumulate,
and what should a low-power logger record to capture it?

All geometry, mass, soil and sea-state values are **illustrative, not site data**. Results are
for learning and demonstration, not for design use.

## Project plan and status

| Phase | Content | Main outputs | Status |
|---|---|---|---|
| A | Structure and soil model | Beam FE model, natural frequencies and mode shapes for 4 installation stages | Complete |
| B | Loading | JONSWAP wave spectra, Morison force spectrum, VIV check | Planned |
| C | Response and fatigue | Response time series, stress at the mudline weld, rainflow counting, S-N curve, Miner damage for a 6-month wait | Planned |
| D | Logger specification and uncertainty | Sensor range, noise limit, sampling rate, logging strategy comparison, Monte Carlo study | Planned |

## Phase A: structure and soil model (complete)

### A-1: cantilever column (validation)

A self-written Euler-Bernoulli beam model of an upright, fixed-base circular tube, validated
against the analytical solution.

![L, D and t vs first natural frequency](figures/L_D_t_vs_f.png)

1. **Length:** frequency falls roughly with 1/L^2.
2. **Diameter:** frequency rises roughly in proportion to D.
3. **Wall thickness:** almost no effect on a bare tube (about 2 % over t = 0.03 to 0.15 m).

Method: f = (beta L)^2 / (2 pi) x sqrt(EI / (m L^4)); FE with consistent mass matrix and a
generalised eigenproblem. With 20 elements the first three modes match the analytical values to
three decimals.

### A-2: realistic monopile, four installation stages

Pile in the seabed (30 m embedded, 30 m water depth), transition piece, tapered tower and
nacelle/rotor. Soil is modelled with distributed springs (k = n_h x depth), and the submerged
part carries hydrodynamic added mass.

![Mode shapes per stage](results/mode_shapes_stages.png)

| Stage | First natural frequency (Hz), medium soil |
|---|---|
| Bare monopile | 1.22 |
| + Transition piece | 0.86 |
| + Tower | 0.51 |
| + Nacelle and rotor | 0.26 |

Key findings:
1. Every installed part lowers the natural frequency.
2. Soft soil and the water around the pile both lower it. Both matter a lot for the bare pile
   and little for the finished turbine.
3. Mode 1 is the important bending shape; higher modes are much faster and far from the wave band.
4. No frequency limit can be chosen; each stage that stays standing must be checked against wave
   (about 0.08 to 0.17 Hz) and wind-vortex (VIV) excitation. The tower without the rotor is the
   stage to watch. This is a hypothesis that Phases B and C will test.

![f1 vs soil stiffness](results/f1_vs_soil_stiffness.png)

Write-ups: [Phase A-1 report](reports/Phase_A_Results_Monopile_Frequency.pdf),
[Phase A-2 report](reports/Phase_A2_Results_Realistic_Monopile.pdf).

## Phase B: loading (planned)

- JONSWAP wave spectra for typical North Sea sea states (Hs 1 to 4 m, Tp 6 to 12 s).
- Linearised Morison wave force and a force spectrum along the pile:
  F = 0.5 rho Cd D u|u| + rho Cm (pi D^2 / 4) du/dt.
- VIV check: critical wind speed from the Strouhal relation, U = f D / St with St about 0.2,
  compared with the natural frequencies from Phase A.

## Phase C: response and fatigue (planned)

- Frequency-domain response, S_x(f) = |H(f)|^2 S_F(f), and time-domain simulation.
- Bending stress at the mudline weld from the mode shape.
- Rainflow counting, an S-N curve (for example the DNV D-curve) and Miner's rule,
  D = sum(n_i / N_i), for damage per sea state.
- Accumulated damage over a 6-month waiting period using a simple sea-state scatter table.

## Phase D: logger specification and uncertainty (planned)

- Define sensor range, low-frequency noise limit and sampling rate (likely 10 to 50 Hz) from the
  response results.
- Compare logging strategies by memory and energy: raw data, statistics, rainflow cycle counts
  and event-triggered recording.
- Monte Carlo study on soil stiffness, damping and wave height to show the spread of the
  damage estimate.

Constants and S-N curve parameters will be verified against the current standard or a textbook
before any numbers are reported.
