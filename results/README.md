# Results: 100 ns of apo and estradiol-bound EdcR

Produced with this repository's pipeline on a CPU-only workstation (Ryzen 9
9950X, 16 cores, GROMACS 2026.3, no GPU), built with

```sh
NTERM_TRIM=14 HMR=yes ./run.sh <apo|holo>
```

i.e. residues 1–14 removed and ACE-capped, ~80,300 atoms per system, and a 4 fs
production step with hydrogen mass repartitioning. Both systems used the same
settings (`settings.sh` and the `*_run.mdp` files in each folder are the exact
inputs) and ran at ~172–178 ns/day with `-ntmpi 4 -ntomp 8 -npme 1`.

| system | status |
|---|---|
| `holo/` | **complete**, 100 ns |
| `apo/` | running; added when it finishes |

## What is in each folder

| file | what it is |
|---|---|
| `summary.txt` | the summary `05_analyse.sh` prints |
| `head_separation.xvg` | recognition-helix (residues 56–63) centre-of-mass distance, every 10 ps, nm |
| `rmsd.xvg`, `rmsf.xvg`, `gyrate.xvg` | backbone RMSD against the start of production, per-residue RMSF, radius of gyration |
| `mindist_EST_*.xvg` | holo only: minimum distance of each estradiol copy to the protein |
| `energy.xvg` | temperature, pressure and density over the run |
| `start.pdb`, `end.pdb` | the solute at 0 and 100 ns, dimer in one piece |
| `traj_100ps.xtc` | the solute trajectory every 100 ps (1,001 frames), dimer in one piece |
| `solute.tpr` | topology matching the trajectory, for GROMACS tools |
| `settings.sh`, `*_run.mdp` | the build options and the exact `.mdp` of each stage |

The full 10 ps trajectory (234 MB per system) is over GitHub's file size
limit and is not included. Every measurement above was made on it, so the
`.xvg` files are at full resolution; only the coordinates are thinned.
Rerunning the pipeline regenerates it.

```sh
gmx rms -s holo/solute.tpr -f holo/traj_100ps.xtc     # any GROMACS tool
vmd holo/start.pdb holo/traj_100ps.xtc                 # or PyMOL: load, then load_traj
```

## Holo, 100 ns

| | |
|---|---|
| head separation | **46.70 ± 1.06 Å** (start 44.00 Å; 10–100 ns: 46.71 ± 0.11 Å, s.e. over 5 ns blocks) |
| backbone RMSD | 0.27 ± 0.05 nm, plateau from ~60 ns |
| radius of gyration | 2.317 ± 0.010 nm |
| estradiol, both copies | 1.75 Å from the protein throughout, never above 2.26 Å |
| temperature / density | 300.0 K / 1,014 kg/m³ |

The heads opened by ~2.7 Å within the first nanosecond after the restraints
came off and stayed there for 100 ns. The closest they came was 43.0 Å, far
from TetR induced (38.88 Å) or TetR on DNA (36.13 Å). Both ligands stayed
bound.

## Read with care

- **Apo and holo start from different AlphaFold predictions**
  (`fold_edcr_alpha3_seed37` model 0 and `fold_edcr_est_seed39` model 2; see
  `input/PROVENANCE.md`), already 1.15 Å apart in head separation. Any
  apo–holo difference mixes the ligand effect with the difference between
  the two models. Simulating the holo model's protein without estradiol would
  separate the two.
- One trajectory per system. Block errors understate the uncertainty of
  slow domain motions; replicas are needed for a claim about estradiol.
- The construct lacks residues 1–14, and HMR slows dynamics somewhat: compare
  conformations and distances, not rates. See the main README.
