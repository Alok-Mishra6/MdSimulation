# Essential dynamics of holo_nt_core_ca.pdb

Principal component analysis of a conformational ensemble (Amadei, Linssen & Berendsen 1993), computed with this project's own
`enm-pca`. Produced by:

```
quiver/build/enm-pca holo_nt_core_ca.pdb --all-chains --modes 20 --compare-anm quiver/validation/edcr/modes/nt_anm15.txt --rmsip 5,10,20 --boltzmann --report holo_full_report.md
```

## 1. The ensemble

| quantity | value |
|---|---|
| frames, M | **1001** |
| nodes, N | 358 |
| coordinates, 3N | 1074 |
| chains | AB |
| superposed | yes, to the running mean |
| superposition iterations | 3 |
| final mean shift | 0.000000512 A |
| RMSD to the mean | 0.5770 - 1.7220 A (mean 0.8852) |

**The rank bound is M - 1 = 1000.** A PCA of 1001 conformations has at most 1000 non-zero components; asking for more
is meaningless, so every number below is capped by it.

## 2. Principal components

Covariance convention **1/M**. This is not cosmetic: the two conventions differ by a factor M/(M-1) = 1.0010 in
every eigenvalue, so a variance quoted without it is not interpretable.

Method `gram` -- the dual form: an 1001 x 1001 eigenproblem instead of a 1074 x 1074 one.

| PC | variance (A^2) | fraction | cumulative | rigid-body leak |
|---|---|---|---|---|
| 1 | 50.381812 | 17.40% | 17.40% | 0.000000 |
| 2 | 38.923525 | 13.45% | 30.85% | 0.000000 |
| 3 | 26.654625 | 9.21% | 40.06% | 0.000000 |
| 4 | 15.545365 | 5.37% | 45.43% | 0.000000 |
| 5 | 11.501471 | 3.97% | 49.40% | 0.000000 |
| 6 | 9.810967 | 3.39% | 52.79% | 0.000000 |
| 7 | 8.665165 | 2.99% | 55.78% | 0.000000 |
| 8 | 7.288160 | 2.52% | 58.30% | 0.000000 |
| 9 | 5.849543 | 2.02% | 60.32% | 0.000000 |
| 10 | 5.337215 | 1.84% | 62.16% | 0.000000 |
| 11 | 4.881276 | 1.69% | 63.85% | 0.000000 |
| 12 | 3.894115 | 1.35% | 65.19% | 0.000000 |
| 13 | 3.585422 | 1.24% | 66.43% | 0.000000 |
| 14 | 3.302939 | 1.14% | 67.57% | 0.000000 |
| 15 | 2.675630 | 0.92% | 68.50% | 0.000000 |
| 16 | 2.540697 | 0.88% | 69.37% | 0.000000 |
| 17 | 2.330417 | 0.80% | 70.18% | 0.000000 |
| 18 | 2.129003 | 0.74% | 70.92% | 0.000000 |
| 19 | 2.078326 | 0.72% | 71.63% | 0.000000 |
| 20 | 1.890319 | 0.65% | 72.29% | 0.000000 |

| | |
|---|---|
| total fluctuation (trace of the covariance) | **289.497495 A^2** |
| captured by the 20 components reported | 209.265990 A^2 (72.29%) |
| components above the zero threshold | 1000 |
| discarded as null space | 1 |
| constant vector inside the Gram null space | 1.000000000000 (must be 1) |
| orthonormality defect, max \|U^T U - I\| | 2.66454e-15 |

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
| 5 | **0.4252** |  |
| 10 | **0.4498** |  |
| 20 | **0.5001** | K exceeds the 14 real component(s) -- not a statement about the model |

The ensemble has **14** component(s) carrying more than 1% of the total fluctuation. RMSIP
compares the *leading* K of each set, so at a K above that the PCA side is
padded with noise directions and the value falls for reasons unrelated to
the elastic network. Read section 3b at those sizes.

### 3b. Per component

| PC | best ANM mode | overlap | cumulative over all 40 modes | note |
|---|---|---|---|---|
| 1 | 2 | 0.4565 | 0.8022 |  |
| 2 | 14 | 0.2732 | 0.5771 |  |
| 3 | 3 | 0.3264 | 0.7012 |  |
| 4 | 5 | 0.2662 | 0.6290 |  |
| 5 | 4 | 0.2630 | 0.4820 |  |
| 6 | 4 | 0.2036 | 0.5180 |  |
| 7 | 6 | 0.2282 | 0.5071 |  |
| 8 | 1 | 0.3086 | 0.5730 |  |
| 9 | 11 | 0.2580 | 0.6077 |  |
| 10 | 4 | 0.2323 | 0.5458 |  |
| 11 | 6 | 0.2832 | 0.5577 |  |
| 12 | 5 | 0.2365 | 0.6524 |  |
| 13 | 1 | 0.2412 | 0.5156 |  |
| 14 | 1 | 0.2586 | 0.6460 |  |
| 15 | 13 | 0.1997 | 0.5063 |  |
| 16 | 9 | 0.2839 | 0.5388 |  |
| 17 | 13 | 0.2378 | 0.5665 |  |
| 18 | 4 | 0.1506 | 0.3963 | **degenerate with PC 19** -- only the subspace is meaningful |
| 19 | 8 | 0.1764 | 0.4035 | **degenerate with PC 18** -- only the subspace is meaningful |
| 20 | 21 | 0.2234 | 0.5076 |  |

### 3c. The spring constant

| quantity | value |
|---|---|
| MSF correlation, predicted vs measured | **0.4278** |
| fitted kT/gamma | 7.078611 A^2 |
| **fitted gamma** | **0.141271 kT/A^2** |

The elastic network predicts fluctuation only up to the spring constant;
the ensemble measures it in angstroms. The least-squares slope between the
two per-node profiles is that constant. This is the standard ANM
calibration, normally done against crystallographic B-factors and done here
against an ensemble instead.

## 4. Warnings

- the ensemble has only 14 component(s) above 1% of the total fluctuation, but RMSIP was requested at larger subspace sizes. RMSIP compares the leading K of each set, so padding the PCA side with noise directions lowers it for reasons unrelated to the model -- read the per-component overlap column instead at those sizes.
- 2 principal component(s) are near-degenerate in variance; individual vectors are ill-defined there and only the subspace they span is meaningful -- compare by RMSIP, not by per-mode overlap
- the fitted gamma has a weak basis: MSF correlation r = 0.427821. The slope is only meaningful where the shapes agree.

---

*Produced by `enm-pca`. Timing: PCA 159.403 ms.*
