# cavb_op05if postprocess results

Outcome: **COMPLETED** (last Total Time 100.00 s)

## BRANCH A -- ran stably to completion

| param | mean | 95% CI | two-halves resolved? |
|---|---|---|---|
| T135 [C] | 201.6 | +/-6.21 | no (stationary) |
| T165 [C] | 162 | +/-4.66 | no (stationary) |
| mdot150 [kg/s] | 0.02305 | +/-0.000142 | no (stationary) |
| q150 [kW/m2] | 1.596 | +/-0.0811 | no (stationary) |
| HRR [kW] | 6.207 | +/-0.104 | no (stationary) |
| qout_40 [kW/m2] | 18.28 | +/-0.394 | no (stationary) |
| qout_80 [kW/m2] | 5.674 | +/-0.357 | no (stationary) |
| qout_100 [kW/m2] | 2.768 | +/-0.206 | no (stationary) |
| qout_135 [kW/m2] | 1.969 | +/-0.0972 | no (stationary) |
| qout_150 [kW/m2] | 1.596 | +/-0.0811 | no (stationary) |
| qout_170 [kW/m2] | 1.212 | +/-0.0499 | no (stationary) |

### vs cavb_op10 (middle grid) -- 'within combined sampling bands' test

| param | op05if - op10 | combined 95% band | within? |
|---|---|---|---|
| T135 | +11.09 | 7.68 | NO |
| T165 | -14.18 | 6.42 | NO |
| mdot150 | -0.001485 | 0.000201 | NO |
| q150 | +0.1999 | 0.0986 | NO |
| HRR | +0.004767 | 0.125 | YES |
| qout_40 | +2.541 | 0.436 | NO |
| qout_80 | +0.5267 | 0.451 | NO |
| qout_100 | +0.1191 | 0.263 | YES |
| qout_135 | +0.1575 | 0.12 | NO |
| qout_150 | +0.1999 | 0.0986 | NO |
| qout_170 | +0.2241 | 0.0614 | NO |

### vs cavb_nb05if (fine-grid no-barrier, same decomposition)

| param | op05if | nb05if | delta |
|---|---|---|---|
| T135 | 201.6 | 163.1 | +38.5 |
| T165 | 162 | 120.9 | +41.19 |
| mdot150 | 0.02305 | 0.02377 | -0.0007209 |
| q150 | 1.596 | 0.9325 | +0.6636 |
| HRR | 6.207 | 6.186 | +0.02087 |
| qout_40 | 18.28 | 18.22 | +0.06103 |
| qout_80 | 5.674 | 4.64 | +1.034 |
| qout_100 | 2.768 | 2.118 | +0.6501 |
| qout_135 | 1.969 | 1.068 | +0.9015 |
| qout_150 | 1.596 | 0.9325 | +0.6636 |
| qout_170 | 1.212 | 0.8149 | +0.3969 |

Branch A slot map: HRR -> [X.XX] kW; T165 -> [XXX.X] C (vs nb05if T165);
mdot150 -> [0.0XXX] kg/s; q150 (and peak qout_40) -> [X.XX] kW/m2.
