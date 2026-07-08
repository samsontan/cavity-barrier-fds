# Grid-Convergence + Stationary-Series Statistics (v2)

Cavity-barrier no-barrier series, FDS. Meshes are **anisotropic**: gap-normal (Y) cells 16.7 / 10 / 5 mm (3/5/10 across the 50 mm gap); X,Z cells 20 / 10 / 5 mm. Refinement ratios are therefore NOT constant. All numbers computed from the FDS output CSVs by `gci_stats_v2.py`.

Steady windows: 60-120 s (nb20, nb10), 60-100 s (nb05). DT_DEVC = 0.5 s.

## Task 1 - Stationary-series uncertainty

Per case per parameter: steady-window mean, sample SD, integral time scale `tau` (ACF integrated to first zero crossing), independent-sample count `N_ind = T_window/(2*tau)`, `SEM = SD/sqrt(N_ind)`, and 95% CI = 1.96*SEM.

### 1a. Mean cross-check against published Table 1

| Case | Param | This script | Table 1 | abs diff | OK? |
|---|---|---|---|---|---|
| nb20 | T @1.35 m (C) | 196.1 | 196.1 | 0.0 | yes |
| nb20 | T @1.65 m (C) | 148.8 | 148.8 | 0.0 | yes |
| nb20 | mdot @1.5 m (kg/s) | 0.0274 | 0.0274 | 0.0000 | yes |
| nb20 | q'' @1.5 m (kW/m2) | 1.109 | 1.11 | 0.001 | yes |
| nb20 | HRR (kW) | 6.109 | 6.11 | 0.001 | yes |
| nb10 | T @1.35 m (C) | 200.1 | 200.1 | 0.0 | yes |
| nb10 | T @1.65 m (C) | 148.1 | 148.1 | 0.0 | yes |
| nb10 | mdot @1.5 m (kg/s) | 0.0261 | 0.0261 | 0.0000 | yes |
| nb10 | q'' @1.5 m (kW/m2) | 0.944 | 0.94 | 0.004 | yes |
| nb10 | HRR (kW) | 6.194 | 6.19 | 0.004 | yes |
| nb05 | T @1.35 m (C) | 179.1 | 179.1 | 0.0 | yes |
| nb05 | T @1.65 m (C) | 135.4 | 135.4 | 0.0 | yes |
| nb05 | mdot @1.5 m (kg/s) | 0.0240 | 0.024 | 0.0000 | yes |
| nb05 | q'' @1.5 m (kW/m2) | 0.974 | 0.97 | 0.004 | yes |
| nb05 | HRR (kW) | 6.207 | 6.21 | 0.003 | yes |

**Verdict:** All means reproduce Table 1 within rounding. Cross-check PASSED - downstream statistics trusted.

### 1b. Full stationary-series statistics

| Case | Param | Mean | SD | tau (s) | N_ind | SEM | 95% CI (+/-) |
|---|---|---|---|---|---|---|---|
| nb20 | T @1.35 m (C) | 196.1 | 39.9 | 0.25 | 120.0 | 3.6 | 7.1 |
| nb20 | T @1.65 m (C) | 148.8 | 30.8 | 0.25 | 120.0 | 2.8 | 5.5 |
| nb20 | mdot @1.5 m (kg/s) | 0.0274 | 0.0014 | 0.25 | 120.0 | 0.0001 | 0.0003 |
| nb20 | q'' @1.5 m (kW/m2) | 1.109 | 0.334 | 0.25 | 120.0 | 0.030 | 0.060 |
| nb20 | HRR (kW) | 6.109 | 0.351 | 0.50 | 60.0 | 0.045 | 0.089 |
| nb10 | T @1.35 m (C) | 200.1 | 48.2 | 0.65 | 45.9 | 7.1 | 14.0 |
| nb10 | T @1.65 m (C) | 148.1 | 29.8 | 0.25 | 120.0 | 2.7 | 5.3 |
| nb10 | mdot @1.5 m (kg/s) | 0.0261 | 0.0012 | 0.25 | 120.0 | 0.0001 | 0.0002 |
| nb10 | q'' @1.5 m (kW/m2) | 0.944 | 0.281 | 0.63 | 47.5 | 0.041 | 0.080 |
| nb10 | HRR (kW) | 6.194 | 0.330 | 0.50 | 60.0 | 0.043 | 0.084 |
| nb05 | T @1.35 m (C) | 179.1 | 58.3 | 0.25 | 80.0 | 6.5 | 12.8 |
| nb05 | T @1.65 m (C) | 135.4 | 38.3 | 0.25 | 80.0 | 4.3 | 8.4 |
| nb05 | mdot @1.5 m (kg/s) | 0.0240 | 0.0013 | 0.25 | 80.0 | 0.0001 | 0.0003 |
| nb05 | q'' @1.5 m (kW/m2) | 0.974 | 0.415 | 0.25 | 80.0 | 0.046 | 0.091 |
| nb05 | HRR (kW) | 6.207 | 0.519 | 0.50 | 40.0 | 0.082 | 0.161 |

**Verdict:** SEM uses the ACF-based independent-sample count, not the raw 121/81 samples, so it honestly reflects temporal correlation in the steady window.

### 1c. Two-halves window-length-bias test (nb05: 60-80 vs 80-100 s)

| Param | Mean 60-80 s | Mean 80-100 s | Diff | Combined 95% band (+/-) | Biased? |
|---|---|---|---|---|---|
| T @1.35 m (C) | 178.9 | 179.3 | +0.4 | 25.7 | no |
| T @1.65 m (C) | 138.8 | 132.1 | -6.7 | 16.8 | no |
| mdot @1.5 m (kg/s) | 0.0239 | 0.0241 | +0.0001 | 0.0006 | no |
| q'' @1.5 m (kW/m2) | 1.009 | 0.940 | -0.069 | 0.240 | no |
| HRR (kW) | 6.190 | 6.224 | +0.035 | 0.327 | no |

**Verdict:** No parameter's half-window shift exceeds its combined 95% band - the nb05 60-100 s window is statistically stationary; window-length bias is negligible.

## Task 2 - Resolved-difference verdicts

For each parameter, is each grid-to-grid step larger than the combined 95% CI of the two grids? `combined = sqrt(CI_a^2 + CI_b^2)`. If a step is inside the band it is sampling noise, not a resolved grid effect.

| Param | 20->10 step | vs band | 10->5 step | vs band | Monotone? | Verdict |
|---|---|---|---|---|---|---|
| T @1.35 m (C) | +4.0 | in band (15.7) | -21.0 | >band (18.9) | NO | non-monotone but within noise |
| T @1.65 m (C) | -0.7 | in band (7.7) | -12.7 | >band (9.9) | yes | monotone, partly within noise |
| mdot @1.5 m (kg/s) | -0.0013 | >band (0.0003) | -0.0021 | >band (0.0004) | yes | monotone & resolved |
| q'' @1.5 m (kW/m2) | -0.165 | >band (0.100) | +0.030 | in band (0.121) | NO | non-monotone but within noise |
| HRR (kW) | +0.086 | in band (0.122) | +0.013 | in band (0.181) | yes | monotone, partly within noise |

**Adjudication of the referee-flagged sequences:**

- **T@1.35 m** (196.1 -> 200.1 -> 179.1): the reversal is NOT statistically resolved - at least one step is inside the combined 95% CI, so within sampling noise. Cannot claim a real non-monotone grid trend.
- **T@1.65 m** (148.8 -> 148.1 -> 135.4): monotone but at least one step is within the 95% CI - trend direction is plausible but not fully resolved.
- **q''@1.5 m** (1.11 -> 0.94 -> 0.97): the reversal is NOT statistically resolved - at least one step is inside the combined 95% CI, so within sampling noise. Cannot claim a real non-monotone grid trend.
- **mdot@1.5 m** (0.0274 -> 0.0261 -> 0.0240): monotone AND both steps exceed the combined 95% CI -> a real, grid-resolved trend.

## Task 3 - GCI redo with variable refinement ratio (Celik 2008)

GCI is computed for every parameter that is **monotone**, using three grid-size definitions: (iso) naive constant r=2 from the X/Z cell 20/10/5 mm; (gap-normal) the Y cell 16.7/10/5 mm giving r32=1.667, r21=2.0; (h-based) the Celik effective size h=(dx dy dz)^(1/3)=18.82/10/5 mm giving r32=1.882, r21=2.0. Apparent order p is solved from the Celik iterative equation with the q(p) term (see appendix).

### T @1.35 m (C)  (coarse=196.1, med=200.1, fine=179.1)

Non-monotone across the three grids -> formal GCI is not defined (solution not in the asymptotic range). Reported as a finding, not a GCI. See Task 2 for whether the non-monotonicity is statistically real.

### T @1.65 m (C)  (coarse=148.8, med=148.1, fine=135.4)

| Grid def | r32 | r21 | p | phi_ext | GCI_fine | e_a | asymptotic? |
|---|---|---|---|---|---|---|---|
| iso r=2 | 2.000 | 2.000 | 4.143 | 134.6420 | 0.7% | 9.4% | NO |
| gap-normal | 1.670 | 2.000 | 3.172 | 133.8213 | 1.5% | 9.4% | NO |
| h-based | 1.882 | 2.000 | 3.782 | 134.4089 | 0.9% | 9.4% | NO |

### mdot @1.5 m (kg/s)  (coarse=0.0274, med=0.0261, fine=0.0240)

| Grid def | r32 | r21 | p | phi_ext | GCI_fine | e_a | asymptotic? |
|---|---|---|---|---|---|---|---|
| iso r=2 | 2.000 | 2.000 | 0.724 | 0.0207 | 17.1% | 8.9% | NO |
| gap-normal | 1.670 | 2.000 | 0.255 | 0.0129 | 57.5% | 8.9% | NO |
| h-based | 1.882 | 2.000 | 0.565 | 0.0195 | 23.2% | 8.9% | NO |

**Comparison to old constant-r=2 numbers** (mass flow): old p=0.72, old GCI=17.0%, old phi_ext=0.0207. The iso row here should reproduce the old figure; the gap-normal / h-based rows show how the anisotropy correction moves p and GCI.

### q'' @1.5 m (kW/m2)  (coarse=1.109, med=0.944, fine=0.974)

Non-monotone across the three grids -> formal GCI is not defined (solution not in the asymptotic range). Reported as a finding, not a GCI. See Task 2 for whether the non-monotonicity is statistically real.

### HRR (kW)  (coarse=6.109, med=6.194, fine=6.207)

| Grid def | r32 | r21 | p | phi_ext | GCI_fine | e_a | asymptotic? |
|---|---|---|---|---|---|---|---|
| iso r=2 | 2.000 | 2.000 | 2.709 | 6.2097 | 0.0% | 0.2% | yes |
| gap-normal | 1.670 | 2.000 | 3.815 | 6.2083 | 0.0% | 0.2% | yes |
| h-based | 1.882 | 2.000 | 3.015 | 6.2092 | 0.0% | 0.2% | yes |

**Verdict:** Because r21 = 2.0 in every definition but r32 differs (1.667 / 1.882 / 2.0), the q(p) correction shifts the apparent order p and the extrapolated value. The honest, anisotropy-consistent numbers are the gap-normal (physically the controlling direction for a 50 mm gap) and h-based rows, not the iso r=2 row the old analysis used.

## Task 4 - Monte-Carlo uncertainty on p and GCI

Each grid mean is perturbed ~ N(mean, SEM) (SEM = 95%CI/1.96), 10,000 draws. For every draw p and GCI_fine are recomputed; non-monotone or degenerate draws are DISCARDED (the discard fraction is itself a convergence-robustness diagnostic). Reported: median and 2.5/97.5 percentiles.

### T @1.65 m (C)

| Grid def | discard frac | p median [2.5, 97.5] | GCI% median [2.5, 97.5] |
|---|---|---|---|
| iso r=2 | 44.0% | 2.05 [0.11, 6.70] | 3.4 [0.1, 78.3] |
| gap-normal | 44.0% | 1.48 [0.09, 5.36] | 5.8 [0.2, 118.8] |
| h-based | 44.0% | 1.84 [0.09, 6.19] | 4.1 [0.1, 95.1] |

### mdot @1.5 m (kg/s)

| Grid def | discard frac | p median [2.5, 97.5] | GCI% median [2.5, 97.5] |
|---|---|---|---|
| iso r=2 | 0.0% | 0.73 [0.23, 1.28] | 17.0 [8.6, 55.3] |
| gap-normal | 0.0% | 0.27 [0.02, 0.75] | 53.8 [18.1, 926.2] |
| h-based | 0.0% | 0.57 [0.11, 1.10] | 23.1 [10.8, 126.6] |

### HRR (kW)

| Grid def | discard frac | p median [2.5, 97.5] | GCI% median [2.5, 97.5] |
|---|---|---|---|
| iso r=2 | 49.0% | 1.22 [0.05, 5.51] | 1.0 [0.0, 43.7] |
| gap-normal | 49.0% | 1.42 [0.06, 6.67] | 0.8 [0.0, 47.9] |
| h-based | 49.0% | 1.26 [0.05, 5.73] | 1.0 [0.0, 46.7] |

**Verdict:** A large discard fraction means the monotonicity itself flips under sampling noise -> the GCI for that parameter is fragile and should be reported with the caveat that the grids are not reliably in the asymptotic range. A small discard fraction means the convergence direction is robust and the p/GCI credible interval can be quoted directly.

## Task 5 - Peak outer-skin heat flux over height

For each nb grid, the steady-window mean of `HFLUX_OUT_Zxxx` is computed at every height; the peak-over-height and its location are reported, with SEM from the ACF method.

| Case | Peak q'' (kW/m2) | SEM | 95% CI (+/-) | Height of peak |
|---|---|---|---|---|
| nb20 | 11.35 | 0.17 | 0.34 | Z0.40 m (40 cm) |
| nb10 | 15.79 | 0.27 | 0.53 | Z0.40 m (40 cm) |
| nb05 | 17.83 | 0.33 | 0.64 | Z0.40 m (40 cm) |

Data sequence (20->10->5 mm): 11.3 -> 15.8 -> 17.8 kW/m2. Old figure claim: 11.4 -> 15.8 -> 17.8 kW/m2.

**Verdict:** peak-over-height is monotone RISING as claimed; data matches the old figure within 0.6 kW/m2. The peak sits at the lowest resolved height near the fire base (see height column), consistent with the expectation that the outer skin is hottest where the plume attaches. This flux is NOT grid converged (it keeps rising as cells shrink), which must be stated as a limitation.

## Appendix - Celik (2008) variable-r formulas implemented

Let f1 = fine, f2 = medium, f3 = coarse grid solutions; grid sizes h1<h2<h3; refinement ratios r21 = h2/h1, r32 = h3/h2 (allowed to differ).

```
eps21 = f2 - f1,   eps32 = f3 - f2
s     = sign(eps32 / eps21)              # +1 monotone, -1 oscillatory

Apparent order p (fixed-point iteration):
    p     = (1/ln(r21)) * | ln|eps32/eps21| + q(p) |
    q(p)  = ln( (r21^p - s) / (r32^p - s) )
  (for constant r, r21=r32 -> q=0 -> p = ln|eps32/eps21| / ln(r), the old formula)

Richardson extrapolated value (fine):
    phi_ext21 = (r21^p * f1 - f2) / (r21^p - 1)

Error estimators:
    e_a21   = |(f1 - f2) / f1|                     # approximate relative
    e_ext21 = |(phi_ext21 - f1) / phi_ext21|       # extrapolated relative
    GCI_fine = 1.25 * e_a21 / (r21^p - 1)          # 1.25 = Fs (3-grid)
```

Grid-size definitions used:
- iso (old, constant r=2): X/Z cell 20 / 10 / 5 mm -> r21=2.0, r32=2.0
- gap-normal (Y cell): 16.7 / 10 / 5 mm -> r21=2.0, r32=1.667
- h-based (dx*dy*dz)^(1/3): 18.82 / 10 / 5 mm -> r21=2.0, r32=1.882

Integral time scale (Task 1): tau = dt * (0.5 + sum_{k=1}^{k0-1} rho(k)), where rho(k) is the normalised autocovariance and k0 the first lag with rho<=0; floored at 0.5*dt. N_ind = T_window / (2*tau); SEM = SD/sqrt(N_ind); 95% CI = 1.96*SEM.

Reference: Celik, Ghia, Roache, Freitas, Coleman, Raad (2008), "Procedure for Estimation and Reporting of Uncertainty Due to Discretization in CFD Applications," J. Fluids Eng. 130(7):078001.
