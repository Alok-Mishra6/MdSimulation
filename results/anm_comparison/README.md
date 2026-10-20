# Does the elastic network predict the MD motion?

Essential dynamics (principal component analysis) of the two 100 ns
trajectories, compared with the anisotropic network model (ANM) modes Quiver
computes for EdcR, using Quiver's own `enm-pca`.

## Inputs

| | |
|---|---|
| ensembles | `results/<system>/traj_100ps.xtc`: 1,001 frames, one per 100 ps |
| ANM modes | Quiver `validation/edcr/modes/nt_anm15.txt` (commit `98b9cd7`): cutoff 15 Å, 40 modes |
| ANM structure | Quiver `validation/edcr/nt_core.pdb`: residues 23–201 of both chains, 358 Cα |

`nt_core.pdb` has the same residue numbering and sequence as both MD
starting structures, and lies 0.70 Å (apo model) and 0.99 Å (holo model)
Cα RMSD from them. The MD construct is longer (residues 15–208), and
`enm-pca` refuses to compare mode sets of different dimension, so the
trajectories were cut to exactly the ANM's 358 Cα first.

## Reproducing it

```sh
# 1. multi-model PDB of the solute, every 100 ps
printf 'System\n' | gmx trjconv -s results/apo/solute.tpr \
    -f results/apo/traj_100ps.xtc -o apo_traj.pdb -dt 100

# 2. the ANM's atoms: CA of residues 23-201, labelled chain A and B
#    (trjconv leaves the chain column blank)
scripts/match_nt_core.py apo_traj.pdb apo_nt_core_ca.pdb

# 3. PCA and the comparison (enm-pca built from Quiver)
enm-pca apo_nt_core_ca.pdb --all-chains --modes 20 \
    --compare-anm nt_anm15.txt --rmsip 5,10,20 --boltzmann --report apo_full_report.md
```

The 10–100 ns variants drop the first 100 models
(`awk '/^MODEL/{m++; keep=(m>100)} keep || /^END$/'`). `--modes 20` is
needed for RMSIP at K = 20. `--boltzmann` states that the ensemble is an
equilibrium sample, which an MD trajectory is; it only affects the fitted
spring constant and the fluctuation comparison, not RMSIP. All four PDBs
are attached to the release
[`edcr-anm-pca-2026-10-20`](https://github.com/Alok-Mishra6/MdSimulation/releases/tag/edcr-anm-pca-2026-10-20).

## Result

| | apo | holo | random subspaces |
|---|---|---|---|
| RMSIP(5) | 0.37 | 0.43 | 0.07 |
| RMSIP(10) | 0.47 | 0.45 | 0.10 |
| RMSIP(20) | 0.50 | 0.50 | 0.14 |
| ANM mode best matching PC1 | mode 2, overlap 0.39 | mode 2, overlap 0.46 | |
| PC1 spanned by all 40 ANM modes | 77 % | 80 % | |
| per-residue fluctuation, MD vs ANM | r = 0.72 | r = 0.43 | |
| total Cα fluctuation | 392 Å² | 289 Å² | |
| PC1 share of it | 31 % | 17 % | |

Full 100 ns. Dropping the first 10 ns moves every number by 0.01 or
less, except holo's fluctuation correlation (0.43 → 0.47); see the
`*_10-100ns_*` files. Random subspaces give RMSIP ≈ √(K/1074).

- The ANM captures the main MD motion, but loosely. RMSIP is 4–7 times
  chance, but at about 0.45 for ten modes it sits below the 0.6–0.7
  usually taken as good agreement.
- The dominant motion of both trajectories lines up best with ANM mode 2,
  which is symmetric across the dimer: the kind of motion that changes the
  separation of the two recognition helices. 77–80 % of PC1 lies inside the
  ANM's 40 lowest modes.
- The ANM reproduces apo's fluctuation profile well (r = 0.72) and holo's
  poorly (r = 0.43). The network is protein-only, and estradiol is not in it.
- Holo fluctuates 26 % less than apo in total, and less of its motion is in
  one collective mode. With the head-separation result in `../README.md`,
  this fits estradiol stiffening the dimer.

## Read with care

- One 100 ns trajectory per system, so the PCA is not converged: only 14
  components carry more than 1 % of the fluctuation. RMSIP(20) is therefore
  padded with noise directions, as `enm-pca` itself warns; for K > 14
  read the per-component overlaps in the reports instead.
- Several PCs are near-degenerate, so individual vectors beyond PC1 are
  ill-defined; compare subspaces (RMSIP), not single PCs.
- The apo–holo differences carry the caveat of `../README.md`: the two
  runs start from different AlphaFold models.

## Files

| file | |
|---|---|
| `<system>_full_report.md` | `enm-pca`'s full report, 0–100 ns: variances, overlap tables, spring constant, warnings |
| `<system>_10-100ns_report.md` | the same, first 10 ns dropped |
| `<system>_*_summary.txt` | the two summary lines `enm-pca` prints |
