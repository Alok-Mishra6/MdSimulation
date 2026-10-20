# Essential dynamics of apo_nt_core_ca.pdb

Principal component analysis of a conformational ensemble (Amadei, Linssen & Berendsen 1993), computed with this project's own
`enm-pca`. Produced by:

```
quiver/build/enm-pca apo_nt_core_ca.pdb --all-chains --modes 20 --compare-anm quiver/validation/edcr/modes/nt_anm15.txt --rmsip 5,10,20 --boltzmann --report apo_full_report.md
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
| final mean shift | 0.000000111 A |
| RMSD to the mean | 0.5647 - 1.8459 A (mean 1.0233) |

**The rank bound is M - 1 = 1000.** A PCA of 1001 conformations has at most 1000 non-zero components; asking for more
is meaningless, so every number below is capped by it.

## 2. Principal components

Covariance convention **1/M**. This is not cosmetic: the two conventions differ by a factor M/(M-1) = 1.0010 in
every eigenvalue, so a variance quoted without it is not interpretable.

Method `gram` -- the dual form: an 1001 x 1001 eigenproblem instead of a 1074 x 1074 one.

| PC | variance (A^2) | fraction | cumulative | rigid-body leak |
|---|---|---|---|---|
| 1 | 122.542806 | 31.30% | 31.30% | 0.000000 |
| 2 | 41.361791 | 10.56% | 41.86% | 0.000000 |
| 3 | 25.543826 | 6.52% | 48.39% | 0.000000 |
| 4 | 17.280944 | 4.41% | 52.80% | 0.000000 |
| 5 | 13.404369 | 3.42% | 56.23% | 0.000000 |
| 6 | 10.601303 | 2.71% | 58.93% | 0.000000 |
| 7 | 9.571972 | 2.44% | 61.38% | 0.000000 |
| 8 | 7.873935 | 2.01% | 63.39% | 0.000000 |
| 9 | 7.335518 | 1.87% | 65.26% | 0.000000 |
| 10 | 5.921363 | 1.51% | 66.78% | 0.000000 |
| 11 | 4.924243 | 1.26% | 68.03% | 0.000000 |
| 12 | 4.783933 | 1.22% | 69.26% | 0.000000 |
| 13 | 4.223619 | 1.08% | 70.33% | 0.000000 |
| 14 | 4.129893 | 1.05% | 71.39% | 0.000000 |
| 15 | 3.910112 | 1.00% | 72.39% | 0.000000 |
| 16 | 3.506123 | 0.90% | 73.28% | 0.000000 |
| 17 | 2.943718 | 0.75% | 74.04% | 0.000000 |
| 18 | 2.751406 | 0.70% | 74.74% | 0.000000 |
| 19 | 2.405852 | 0.61% | 75.35% | 0.000000 |
| 20 | 2.336116 | 0.60% | 75.95% | 0.000000 |

| | |
|---|---|
| total fluctuation (trace of the covariance) | **391.511802 A^2** |
| captured by the 20 components reported | 297.352841 A^2 (75.95%) |
| components above the zero threshold | 1000 |
| discarded as null space | 1 |
| constant vector inside the Gram null space | 1.000000000000 (must be 1) |
| orthonormality defect, max \|U^T U - I\| | 2.44249e-15 |

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
| 5 | **0.3715** |  |
| 10 | **0.4724** |  |
| 20 | **0.4984** | K exceeds the 14 real component(s) -- not a statement about the model |

The ensemble has **14** component(s) carrying more than 1% of the total fluctuation. RMSIP
compares the *leading* K of each set, so at a K above that the PCA side is
padded with noise directions and the value falls for reasons unrelated to
the elastic network. Read section 3b at those sizes.

### 3b. Per component

| PC | best ANM mode | overlap | cumulative over all 40 modes | note |
|---|---|---|---|---|
| 1 | 2 | 0.3902 | 0.7677 |  |
| 2 | 1 | 0.3703 | 0.6875 |  |
| 3 | 22 | 0.2670 | 0.6527 |  |
| 4 | 9 | 0.2647 | 0.6262 |  |
| 5 | 3 | 0.3097 | 0.5921 |  |
| 6 | 4 | 0.3367 | 0.6155 |  |
| 7 | 5 | 0.3170 | 0.5834 |  |
| 8 | 1 | 0.2483 | 0.5226 |  |
| 9 | 9 | 0.2182 | 0.5469 |  |
| 10 | 2 | 0.4526 | 0.7013 |  |
| 11 | 12 | 0.3241 | 0.5296 | **degenerate with PC 12** -- only the subspace is meaningful |
| 12 | 13 | 0.2859 | 0.5411 | **degenerate with PC 11** -- only the subspace is meaningful |
| 13 | 15 | 0.2333 | 0.6323 | **degenerate with PC 14** -- only the subspace is meaningful |
| 14 | 3 | 0.2748 | 0.5165 | **degenerate with PC 13** -- only the subspace is meaningful |
| 15 | 12 | 0.2162 | 0.4297 |  |
| 16 | 12 | 0.2851 | 0.5617 |  |
| 17 | 2 | 0.1538 | 0.5270 |  |
| 18 | 9 | 0.2038 | 0.5687 |  |
| 19 | 1 | 0.1715 | 0.3837 | **degenerate with PC 20** -- only the subspace is meaningful |
| 20 | 13 | 0.2677 | 0.4717 | **degenerate with PC 19** -- only the subspace is meaningful |

### 3c. The spring constant

| quantity | value |
|---|---|
| MSF correlation, predicted vs measured | **0.7226** |
| fitted kT/gamma | 9.790291 A^2 |
| **fitted gamma** | **0.102142 kT/A^2** |

The elastic network predicts fluctuation only up to the spring constant;
the ensemble measures it in angstroms. The least-squares slope between the
two per-node profiles is that constant. This is the standard ANM
calibration, normally done against crystallographic B-factors and done here
against an ensemble instead.

## 4. Warnings

- the ensemble has only 14 component(s) above 1% of the total fluctuation, but RMSIP was requested at larger subspace sizes. RMSIP compares the leading K of each set, so padding the PCA side with noise directions lowers it for reasons unrelated to the model -- read the per-component overlap column instead at those sizes.
- 6 principal component(s) are near-degenerate in variance; individual vectors are ill-defined there and only the subspace they span is meaningful -- compare by RMSIP, not by per-mode overlap

---

*Produced by `enm-pca`. Timing: PCA 195.707 ms.*
