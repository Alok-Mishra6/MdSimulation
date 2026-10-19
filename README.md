# EdcR molecular dynamics — a self-contained GROMACS setup

Part of **Quiver**, Team IISER-Kolkata's iGEM 2026 software project.


Two systems, built from AlphaFold 3 predictions of the EdcR homodimer:

| system | what it is | atoms |
|---|---|---|
| **apo** | the dimer alone, ligand-free | 381,383 |
| **holo** | the same dimer with **two copies of 17β-estradiol** bound | 306,178 |

Hand this folder to anyone with a GPU and GROMACS. Everything needed is in
here — structures, force-field parameters for the ligand, run parameters and
scripts. Nothing is downloaded at run time.

```sh
./scripts/00_check_env.sh      # does this machine have what it needs?
./run.sh both                  # build, equilibrate, run, analyse
```

---

## Why two systems

EdcR is a repressor: it grips its DNA operator when no estrogen is present and
lets go when estradiol binds. The pair is the experiment. Anything that is
genuinely caused by the ligand has to show up as a *difference* between these
two runs; anything that appears in both is the force field, the water, or the
starting structure.

The two starting structures have the **identical sequence and identical residue
numbering (1–208, two chains)**, so they can be compared residue by residue
without any alignment step.

## What you need

| | |
|---|---|
| GROMACS | **2023 or newer**, built with CUDA. Developed and tested on 2026.3 |
| a GPU | NVIDIA, any modern card. ~1.5 GB of VRAM is enough. Optional: see "Without a GPU" |
| disk | ~400 MB per 100 ns per system (the trajectory stores the solute only; see `mdp/md.mdp`) |
| Python 3 | only for three small helper scripts; no packages needed |

GROMACS is **not** bundled — it is a 45 MB source tree and a 20-minute build,
and your cluster almost certainly has it. Point the scripts at it either way:

```sh
module load gromacs                        # most HPC systems
GROMACS_ROOT=/path/to/gromacs ./run.sh apo # or name it explicitly
```

The force field itself needs no download: **amber99sb-ildn ships inside
GROMACS**. That was the deciding reason for choosing it.

## Running it

```sh
./run.sh apo                 # one system, 100 ns
./run.sh holo
./run.sh both                # both, in sequence
NS=10 ./run.sh both          # 10 ns instead of 100
./run.sh apo --build         # stop after building, inspect, continue later
./run.sh apo --analyse       # re-measure an existing trajectory
```

Each stage can also be run on its own, in order:

| script | what it does | output |
|---|---|---|
| `01_build.sh` | pdb2gmx → box → solvate → ions → index groups | `system.gro`, `topol.top` |
| `02_minimise.sh` | steepest descent until Fmax < 1000 | `em.gro` |
| `03_equilibrate.sh` | 100 ps NVT at 300 K, then 100 ps NPT at 1 bar | `npt.gro` |
| `04_production.sh` | unrestrained production MD | **`md.xtc`** |
| `05_analyse.sh` | RMSD, RMSF, Rg, head separation, ligand drift | `analysis/` |

Everything lands in `runs/apo/` and `runs/holo/`. The trajectory you want is
`runs/<system>/md.xtc`, and `runs/<system>/analysis/` holds the measurements.

**If a run is interrupted, just run `04_production.sh` again** — it passes
`-cpi` and picks up from the last checkpoint rather than starting over.

## How long it takes

Production is 100 ns by default. Rough throughput for a ~300,000-atom system:

| GPU | ns/day | 100 ns takes |
|---|---|---|
| A100 / H100 | 100–200 (estimated) | under a day |
| RTX 3090 / 4090 | 60–120 (estimated) | 1–2 days |
| **RTX 2050, laptop** | **16.5 (measured)** | **~6 days** |

The RTX 2050 figure is measured, not guessed: 16.483 ns/day for the 306,178-atom
holo system, from the equilibration log on the machine this was built on. That
is a low-end laptop GPU, so treat it as the floor. The other rows are scaled
from it and are estimates.

Both systems at once is twice that. `NS=10` first is a sensible way to confirm
the setup on your hardware before committing a week to it.

The runs use GROMACS' **GPU-resident path** — non-bonded, PME, bonded *and* the
integration all run on the GPU (`-nb gpu -pme gpu -bonded gpu -update gpu`),
which is typically worth a factor of two on a single card. The barostat is
`c-rescale` specifically because it is the one GROMACS supports alongside
`-update gpu`; it samples the correct NPT ensemble, so the speed costs nothing
physically.

### Making it 12–15x faster: `NTERM_TRIM` and `HMR`

```sh
NTERM_TRIM=14 HMR=yes ./run.sh both
```

Two build-time options, both off by default so the full-length setup verified
below stays the default:

- **`NTERM_TRIM=14`** removes residues 1–14 of each chain and caps residue 15
  with an acetyl group (`scripts/trim_nterm.py`). Those residues are a
  low-confidence tail (pLDDT 42–56) that AlphaFold left pointing more than 4 nm
  away from the core, and they alone set the box: the solute is 14.9 nm across
  with them and 8.0 nm without. The systems shrink from 381,383 / 306,178 atoms
  to **80,330 / 80,298**. The recognition helices (56–63) keep their residue
  numbers, and both estradiol copies sit ≥ 26 Å from any removed atom. It is a
  change of construct, though: residues 12–19 are Arg-rich and might contact
  DNA in the real complex, which this box does not contain anyway.
- **`HMR=yes`** sets `mass-repartition-factor = 3`: hydrogens become 3x heavier,
  the mass taken from the atom each is bonded to, and production runs at
  **4 fs**. Equilibrium properties do not depend on the masses; dynamics are
  somewhat slowed, so do not read rates off these trajectories. GROMACS applies
  it to water as well. Equilibration stays at 2 fs. Expect LINCS warnings at
  4 fs, all on the methionine S-CH3 methyls (about 740 an hour for this dimer,
  constraint error under 1.5 %); `04_production.sh` raises mdrun's abort limit
  of 1000 warnings so the run does not stop every 80 minutes.

The options are written to `runs/<system>/settings.sh` by `01_build.sh`, and
every later stage reads them from there, so a run cannot be half one and half
the other. The `.mdp` each stage actually ran is kept as `runs/<system>/*_run.mdp`.

### Without a GPU

The scripts check how GROMACS was built. A CPU-only build gets no offload flags
and runs everything on the CPU; `USE_GPU=no` forces that path on a GPU build.
Measured on a 16-core Ryzen 9 9950X (AVX-512, GROMACS 2026.3, no GPU), ns/day:

| system | layout | apo | holo |
|---|---|---|---|
| full length, 2 fs | all 32 threads, mdrun's own choice | 12.1 | 14.4 |
| full length, 2 fs | both at once, one per 8-core CCD | 6.2 | 7.0 |
| **trimmed + HMR, 4 fs** | **all 32 threads, `-ntmpi 4 -ntomp 8 -npme 1`** | **177** | **170** |
| trimmed + HMR, 4 fs | all 32 threads, mdrun's own choice | 139 | 137 |
| trimmed + HMR, 4 fs | both at once, one per CCD | 76 | 73 |

So on this CPU the full-length pair takes about 15 days for 100 ns each, and
the trimmed HMR pair about 14 hours each. One system at a time on the whole
machine beat one per CCD for the small systems; mdrun's default layout leaves
~25 % on the table, so pass the layout explicitly through `GMX_MDRUN_EXTRA`.

## Simulation parameters

Chosen to be unremarkable, because a comparison between two runs is only
meaningful if neither does anything unusual.

| | |
|---|---|
| protein force field | amber99sb-ildn |
| ligand | GAFF2 with AM1-BCC charges (see `ligand/README.md`) |
| water | TIP3P |
| box | rhombic dodecahedron, 1.2 nm to the edge |
| ions | neutralising Na⁺/Cl⁻ plus 0.15 M NaCl |
| temperature | 300 K, v-rescale thermostat, τ = 0.1 ps |
| pressure | 1 bar, c-rescale barostat, τ = 2 ps |
| cutoffs | PME, 1.0 nm Coulomb and van der Waals, dispersion correction on |
| time step | 2 fs with h-bond constraints (LINCS) |
| frames | every 10 ps |

The cutoffs are the values AMBER was parameterised with, and they are identical
in every stage. Changing them changes the force field.

## What comes out, and the one measurement that matters here

`05_analyse.sh` writes the usual diagnostics — backbone RMSD, per-residue RMSF,
radius of gyration — and then one project-specific number: **the distance
between the two recognition helices**, the DNA-reading heads.

That is the quantity this whole EdcR analysis turns on:

| | separation |
|---|---|
| AlphaFold 3's prediction, across 55 models | **43.97 ± 0.72 Å** |
| a real induced TetR (2TRT) | 38.88 Å |
| a real TetR gripping its operator (1QPI) | 36.13 Å |

AlphaFold puts the heads about 8 Å wider than a repressor actually bound to
DNA, and the ligand moved them by 0.007 Å — a flat result that AlphaFold's own
error estimate (6.6 Å on that very geometry) cannot resolve. **The open
question is whether the heads close when the structure is free to move.** That
is what these trajectories are for.

The measurement uses residues **56–63**, which is the EdcR segment that
structurally superposes onto TetR's α3 recognition helix. This is not the
project's original 55–60 guess from secondary structure — the matched range was
established in `validation/af3/06_convention.md` and the difference is about
1 Å, which matters for any comparison against the crystal structures above.

`scripts/head_groups.py` builds the `HEAD_A` / `HEAD_B` index groups, and it
refuses to run rather than guess if the residue numbering is not what it
expects.

## What has actually been verified

The pipeline was run end to end on the holo system — build, minimise,
equilibrate, 100 ps of production, analysis — on a laptop RTX 2050. What came
out:

| check | result |
|---|---|
| system charge | **0.000** |
| minimisation | converged, Fmax < 1000 in 2,632 steps |
| temperature after NVT | 296.9 K (target 300) |
| water density after NPT | **991.6 kg/m³** — correct for TIP3P |
| GPU path | non-bonded, PME, bonded **and** integration all on the GPU |
| throughput | 16.6 ns/day |
| backbone RMSD over 100 ps | 0.172 nm |
| radius of gyration | 7.118 ± 0.008 nm — stable |
| estradiol, copy 1 | **1.76 ± 0.09 Å** from the protein — in contact |
| estradiol, copy 2 | **1.80 ± 0.09 Å** — in contact |

Both systems were run: apo at 12.8 ns/day (381,383 atoms, density 990.4 kg/m³)
and holo at 16.6 ns/day (306,178 atoms, density 991.6 kg/m³).

Three of those deserve comment.

**The water density.** 990–992 kg/m³ is the single best evidence that the force
field, the water model and the barostat are all behaving. It is not a number you
hit by accident.

**The ligand.** Sitting at 1.8 Å van der Waals contact and unmoved across every
frame says the GAFF2 parameters and the predicted pose are at least
self-consistent — a badly parameterised ligand drifts out of its pocket almost
at once. Each copy is checked separately, so one escaping cannot be averaged
away by the other.

**The head separation agrees with the project's own analysis code.** This matters
most, because it is the measurement everything else rests on, and the two
implementations share nothing: `gmx distance` over an index group here, versus
`scripts/convention_check.py` in pure Python there.

| | project's script, raw prediction | this pipeline, after minimisation + equilibration |
|---|---|---|
| apo | 45.12 Å | 45.30 Å |
| holo | 43.98 Å | 44.02 Å |

They agree to **0.04 Å on the holo system and 0.18 Å on the apo one**, and the
residual is real relaxation, not disagreement: hydrogens were added and the
structure was minimised and equilibrated in between. The structure, the chain
split, the 56–63 selection and the centroid convention are therefore all wired
up the same way as in the published analysis — so a change seen in these
trajectories can be compared against it directly.

## Honest limits

- **No *full-length* production run has been done.** Everything above is from a
  100 ps test. The 100 ns runs are what this folder is *for*, not what it
  contains.
- **100 ns probably is not long enough** to see a domain rearrangement of ~8 Å.
  If the heads do not close, that is not evidence they cannot; it is a run that
  was too short. Treat 100 ns as a first look and expect to need microseconds,
  or enhanced sampling, for a real answer.
- **Two runs are not statistics.** For any claim about an estradiol effect,
  run three or more independent replicas of each system (change `gen_seed` in
  `mdp/nvt.mdp` from `-1` to different integers) and compare distributions, not
  single trajectories.
- **The starting structures are predictions, not experiments.** The apo model's
  own DNA-binding helix H1 has the lowest pLDDT of the nine helices, and the
  inter-domain PAE is 6.6 Å. The heads may simply have been placed arbitrarily.
  MD started from an uncertain structure inherits that uncertainty.
- **The ligand parameters are GAFF2/AM1-BCC**, which is a general small-molecule
  force field, not a steroid-specific one. `ligand/README.md` lists exactly what
  was substituted by analogy.
- **No membrane, no DNA, no crowding.** A dimer in a box of water.

## Provenance

> **Where the references to `validation/af3/...` point.** This repository is
> the molecular-dynamics setup only. Those paths refer to the main project
> repository, where the AlphaFold analysis that produced these structures
> lives: <https://gitlab.igem.org/2026/software/iiser-kolkata/quiver>. Nothing
> here needs that repository in order to run.

See `input/PROVENANCE.md` for exactly which AlphaFold job and model each
structure came from, with checksums, and `ligand/README.md` for how the
estradiol parameters were generated and how to regenerate them.
