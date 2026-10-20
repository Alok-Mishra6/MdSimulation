# Essential dynamics of apo_nt_core_ca_10-100ns.pdb

Principal component analysis of a conformational ensemble (Amadei, Linssen & Berendsen 1993), computed with this project's own
`enm-pca`. Produced by:

```
quiver/build/enm-pca apo_nt_core_ca_10-100ns.pdb --all-chains --modes 20 --compare-anm quiver/validation/edcr/modes/nt_anm15.txt --rmsip 5,10,20 --boltzmann --report apo_10-100_report.md
```

## 1. The ensemble

| quantity | value |
|---|---|
| frames, M | **901** |
| nodes, N | 358 |
| coordinates, 3N | 1074 |
| chains | AB |
| superposed | yes, to the running mean |
| superposition iterations | 3 |
| final mean shift | 0.000000221 A |
| RMSD to the mean | 0.5622 - 1.8355 A (mean 1.0134) |

**The rank bound is M - 1 = 900.** A PCA of 901 conformations has at most 900 non-zero components; asking for more
is meaningless, so every number below is capped by it.

## 2. Principal components

Covariance convention **1/M**. This is not cosmetic: the two conventions differ by a factor M/(M-1) = 1.0011 in
every eigenvalue, so a variance quoted without it is not interpretable.

Method `gram` -- the dual form: an 901 x 901 eigenproblem instead of a 1074 x 1074 one.

| PC | variance (A^2) | fraction | cumulative | rigid-body leak |
|---|---|---|---|---|
| 1 | 130.860806 | 34.07% | 34.07% | 0.000000 |
| 2 | 31.507773 | 8.20% | 42.27% | 0.000000 |
| 3 | 27.309854 | 7.11% | 49.38% | 0.000000 |
| 4 | 17.384253 | 4.53% | 53.90% | 0.000000 |
| 5 | 14.285340 | 3.72% | 57.62% | 0.000000 |
| 6 | 10.385213 | 2.70% | 60.33% | 0.000000 |
| 7 | 9.001298 | 2.34% | 62.67% | 0.000000 |
| 8 | 7.863977 | 2.05% | 64.72% | 0.000000 |
| 9 | 5.958491 | 1.55% | 66.27% | 0.000000 |
| 10 | 5.134316 | 1.34% | 67.60% | 0.000000 |
| 11 | 4.796570 | 1.25% | 68.85% | 0.000000 |
| 12 | 4.311784 | 1.12% | 69.97% | 0.000000 |
| 13 | 3.894846 | 1.01% | 70.99% | 0.000000 |
| 14 | 3.828082 | 1.00% | 71.98% | 0.000000 |
| 15 | 3.780905 | 0.98% | 72.97% | 0.000000 |
| 16 | 2.859361 | 0.74% | 73.71% | 0.000000 |
| 17 | 2.644311 | 0.69% | 74.40% | 0.000000 |
| 18 | 2.585254 | 0.67% | 75.07% | 0.000000 |
| 19 | 2.447375 | 0.64% | 75.71% | 0.000000 |
| 20 | 2.251784 | 0.59% | 76.30% | 0.000000 |

| | |
|---|---|
| total fluctuation (trace of the covariance) | **384.141265 A^2** |
| captured by the 20 components reported | 293.091594 A^2 (76.30%) |
| components above the zero threshold | 900 |
| discarded as null space | 1 |
| constant vector inside the Gram null space | 1.000000000000 (must be 1) |
| orthonormality defect, max \|U^T U - I\| | 3.42651e-15 |

The orthonormality defect is a diagnostic, not a flaw. The dual form lifts each
component as `u = D w / sqrt(denom * mu)`, so rounding in `D w` is amplified by
`1/sqrt(mu)` and orthogonality degrades roughly as `sqrt(mu_max/mu)`. A large
value means more components were asked for than the ensemble supports.

**The rigid-body leak column is a check no reference implementation performs.**
After superposition to the mean the deviations are orthogonal to global
translation and rotation at first order, so a leading component with a large
leak means the superposition was wrong and everything here is meaningless.

## 3. Against the elastic network

The canonical test of an elastic network model: does it predict the motion
the ensemble measures? Three quantities, which are **not** the same and are
never collapsed here:

| name | definition |
|---|---|
| overlap | `\|u_i . v_k\|`, a cosine between PC i and ANM mode k |
| cumulative overlap | `sqrt(sum_{k<=K} (u_i.v_k)^2)` -- how much of PC i the first K modes span |
| RMSIP | `sqrt((1/K) sum_i sum_k (u_i.v_k)^2)` -- agreement of the two subspaces |

### 3a. RMSIP

| subspace size K | RMSIP | |
|---|---|---|
| 5 | **0.3649** |  |
| 10 | **0.4661** |  |
| 20 | **0.5000** | K exceeds the 13 real component(s) -- not a statement about the model |

The ensemble has **13** component(s) carrying more than 1% of the total fluctuation. RMSIP
compares the *leading* K of each set, so at a K above that the PCA side is
padded with noise directions and the value falls for reasons unrelated to
the elastic network. Read section 3b at those sizes.

### 3b. Per component

| PC | best ANM mode | overlap | cumulative over all 40 modes | note |
|---|---|---|---|---|
| 1 | 2 | 0.3886 | 0.7578 |  |
| 2 | 1 | 0.3381 | 0.7175 |  |
| 3 | 22 | 0.2554 | 0.6477 |  |
| 4 | 9 | 0.2423 | 0.5889 |  |
| 5 | 3 | 0.3293 | 0.5850 |  |
| 6 | 4 | 0.3063 | 0.5903 |  |
| 7 | 5 | 0.3484 | 0.6540 |  |
| 8 | 1 | 0.2816 | 0.5110 |  |
| 9 | 2 | 0.4165 | 0.6892 |  |
| 10 | 13 | 0.2955 | 0.5308 |  |
| 11 | 12 | 0.2125 | 0.5315 |  |
| 12 | 11 | 0.1975 | 0.5234 |  |
| 13 | 12 | 0.3398 | 0.6582 | **degenerate with PC 14** -- only the subspace is meaningful |
| 14 | 12 | 0.1520 | 0.4103 | **degenerate with PC 15** -- only the subspace is meaningful |
| 15 | 3 | 0.1983 | 0.5340 | **degenerate with PC 14** -- only the subspace is meaningful |
| 16 | 8 | 0.1538 | 0.4966 |  |
| 17 | 10 | 0.2437 | 0.5808 | **degenerate with PC 18** -- only the subspace is meaningful |
| 18 | 13 | 0.1595 | 0.4276 | **degenerate with PC 17** -- only the subspace is meaningful |
| 19 | 13 | 0.2081 | 0.4972 |  |
| 20 | 10 | 0.2692 | 0.5053 |  |

### 3c. The spring constant

| quantity | value |
|---|---|
| MSF correlation, predicted vs measured | **0.7226** |
| fitted kT/gamma | 9.671658 A^2 |
| **fitted gamma** | **0.103395 kT/A^2** |

The elastic network predicts fluctuation only up to the spring constant;
the ensemble measures it in angstroms. The least-squares slope between the
two per-node profiles is that constant. This is the standard ANM
calibration, normally done against crystallographic B-factors and done here
against an ensemble instead.

## 4. Warnings

- the ensemble has only 13 component(s) above 1% of the total fluctuation, but RMSIP was requested at larger subspace sizes. RMSIP compares the leading K of each set, so padding the PCA side with noise directions lowers it for reasons unrelated to the model -- read the per-component overlap column instead at those sizes.
- 5 principal component(s) are near-degenerate in variance; individual vectors are ill-defined there and only the subspace they span is meaningful -- compare by RMSIP, not by per-mode overlap

---

*Produced by `enm-pca`. Timing: PCA 137.878 ms.*
