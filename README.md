# Phase A:Natural Frequency of a Fixed Circular Cantilever Column

A Python study of how geometry controls the first natural frequency of an
upright, fixed-base circular tube (Euler-Bernoulli beam theory), validated
against the analytical solution with a self-written finite element model.

## Observations

![L, D and t vs first natural frequency](figures/L_D_t_vs_f.png)

1. **Length:** frequency decreases as length increases, roughly with 1/L^2.
   A longer column is much softer (stiffness falls with L^3) and also heavier.
2. **Diameter:** frequency increases with diameter, roughly in proportion to D.
   Bending stiffness grows with about D^3 while mass grows only with D.
3. **Wall thickness:** almost no effect on a bare tube. A thicker wall adds
   stiffness and mass in nearly equal proportion, so they cancel out
   (frequency changes by about 2 % over t = 0.03 to 0.15 m).

## Method

- Analytical solution: f = (beta L)^2 / (2 pi) x sqrt(EI / (m L^4)), beta L = 1.8751 for the first mode
- FE model: Euler-Bernoulli beam elements, consistent mass matrix, generalised eigenproblem
- Validation: 20 elements reproduce the analytical first three modes to three decimals

## Run

```
pip install numpy scipy matplotlib
python src/05_L_D_t_vs_f.py
```

## Notes

Dimensions and material values are illustrative.# riyaTR04.github.io
