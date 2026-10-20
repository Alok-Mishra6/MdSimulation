# Essential dynamics of holo_nt_core_ca_10-100ns.pdb

Principal component analysis of a conformational ensemble (Amadei, Linssen & Berendsen 1993), computed with this project's own
`enm-pca`. Produced by:

```
quiver/build/enm-pca holo_nt_core_ca_10-100ns.pdb --all-chains --modes 20 --compare-anm quiver/validation/edcr/modes/nt_anm15.txt --rmsip 5,10,20 --boltzmann --report holo_10-100_report.md
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
| final mean shift | 0.000000468 A |
| RMSD to the mean | 0.5675 - 1.6589 A (mean 0.8717) |

**The rank bound is M - 1 = 900.** A PCA of 901 conformations has at most 900 non-zero components; asking for more
is meaningless, so every number below is capped by it.

## 2. Principal components

Covariance convention **1/M**. This is not cosmetic: the two conventions differ by a factor M/(M-1) = 1.0011 in
every eigenvalue, so a variance quoted without it is not interpretable.

Method `gram` -- the dual form: an 901 x 901 eigenproblem instead of a 1074 x 1074 one.

| PC | variance (A^2) | fraction | cumulative | rigid-body leak |
|---|---|---|---|---|
| 1 | 51.222554 | 18.25% | 18.25% | 0.000000 |
| 2 | 33.806799 | 12.05% | 30.30% | 0.000000 |
| 3 | 26.443268 | 9.42% | 39.73% | 0.000000 |
| 4 | 15.411103 | 5.49% | 45.22% | 0.000000 |
| 5 | 11.148719 | 3.97% | 49.19% | 0.000000 |
| 6 | 9.223968 | 3.29% | 52.48% | 0.000000 |
| 7 | 8.226183 | 2.93% | 55.41% | 0.000000 |
| 8 | 7.173581 | 2.56% | 57.97% | 0.000000 |
| 9 | 5.781142 | 2.06% | 60.03% | 0.000000 |
| 10 | 5.434555 | 1.94% | 61.97% | 0.000000 |
| 11 | 4.897594 | 1.75% | 63.71% | 0.000000 |
| 12 | 3.879362 | 1.38% | 65.09% | 0.000000 |
| 13 | 3.459400 | 1.23% | 66.33% | 0.000000 |
| 14 | 3.209136 | 1.14% | 67.47% | 0.000000 |
| 15 | 2.712227 | 0.97% | 68.44% | 0.000000 |
| 16 | 2.497584 | 0.89% | 69.33% | 0.000000 |
| 17 | 2.177624 | 0.78% | 70.10% | 0.000000 |
| 18 | 2.138152 | 0.76% | 70.86% | 0.000000 |
| 19 | 1.898024 | 0.68% | 71.54% | 0.000000 |
| 20 | 1.833309 | 0.65% | 72.19% | 0.000000 |

| | |
|---|---|
| total fluctuation (trace of the covariance) | **280.596061 A^2** |
| captured by the 20 components reported | 202.574285 A^2 (72.19%) |
| components above the zero threshold | 900 |
| discarded as null space | 1 |
| constant vector inside the Gram null space | 1.000000000000 (must be 1) |
| orthonormality defect, max \|U^T U - I\| | 3.84935e-15 |

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
| 5 | **0.4166** |  |
| 10 | **0.4444** |  |
| 20 | **0.4972** | K exceeds the 14 real component(s) -- not a statement about the model |

The ensemble has **14** component(s) carrying more than 1% of the total fluctuation. RMSIP
compares the *leading* K of each set, so at a K above that the PCA side is
padded with noise directions and the value falls for reasons unrelated to
the elastic network. Read section 3b at those sizes.

### 3b. Per component

| PC | best ANM mode | overlap | cumulative over all 40 modes | note |
|---|---|---|---|---|
| 1 | 2 | 0.4167 | 0.8088 |  |
| 2 | 14 | 0.2827 | 0.5524 |  |
| 3 | 7 | 0.2706 | 0.7211 |  |
| 4 | 3 | 0.3228 | 0.6376 |  |
| 5 | 4 | 0.3104 | 0.5286 |  |
| 6 | 3 | 0.2354 | 0.5404 |  |
| 7 | 5 | 0.2116 | 0.5244 |  |
| 8 | 1 | 0.2881 | 0.5373 |  |
| 9 | 11 | 0.2652 | 0.5579 |  |
| 10 | 4 | 0.1978 | 0.5476 |  |
| 11 | 6 | 0.2329 | 0.5259 |  |
| 12 | 15 | 0.2561 | 0.6373 |  |
| 13 | 1 | 0.3545 | 0.5683 |  |
| 14 | 6 | 0.2170 | 0.5751 |  |
| 15 | 13 | 0.2274 | 0.5139 |  |
| 16 | 4 | 0.2172 | 0.5913 |  |
| 17 | 13 | 0.1973 | 0.4300 | **degenerate with PC 18** -- only the subspace is meaningful |
| 18 | 8 | 0.1639 | 0.4764 | **degenerate with PC 17** -- only the subspace is meaningful |
| 19 | 21 | 0.1917 | 0.4360 | **degenerate with PC 20** -- only the subspace is meaningful |
| 20 | 13 | 0.2676 | 0.5349 | **degenerate with PC 19** -- only the subspace is meaningful |

### 3c. The spring constant

| quantity | value |
|---|---|
| MSF correlation, predicted vs measured | **0.4748** |
| fitted kT/gamma | 6.987316 A^2 |
| **fitted gamma** | **0.143116 kT/A^2** |

The elastic network predicts fluctuation only up to the spring constant;
the ensemble measures it in angstroms. The least-squares slope between the
two per-node profiles is that constant. This is the standard ANM
calibration, normally done against crystallographic B-factors and done here
against an ensemble instead.

## 4. Warnings

- the ensemble has only 14 component(s) above 1% of the total fluctuation, but RMSIP was requested at larger subspace sizes. RMSIP compares the leading K of each set, so padding the PCA side with noise directions lowers it for reasons unrelated to the model -- read the per-component overlap column instead at those sizes.
- 4 principal component(s) are near-degenerate in variance; individual vectors are ill-defined there and only the subspace they span is meaningful -- compare by RMSIP, not by per-mode overlap
- the fitted gamma has a weak basis: MSF correlation r = 0.474807. The slope is only meaningful where the shapes agree.

---

*Produced by `enm-pca`. Timing: PCA 109.156 ms.*
