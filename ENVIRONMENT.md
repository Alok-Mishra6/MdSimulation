# The environment this was built and tested in

Recorded so that a difference between your results and ours is traceable to
something, rather than a mystery. **None of this is required to run the
simulation** — any GROMACS 2023+ with CUDA will do.

## GROMACS

Not in the Arch repositories, so it was built from source. No `sudo` and no
system files touched: it installs under `$HOME`.

```sh
curl -O https://ftp.gromacs.org/gromacs/gromacs-2026.3.tar.gz
tar xzf gromacs-2026.3.tar.gz && cd gromacs-2026.3

cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=gcc-15 -DCMAKE_CXX_COMPILER=g++-15 \
  -DGMX_GPU=CUDA \
  -DCMAKE_CUDA_ARCHITECTURES=86 \
  -DCMAKE_CUDA_HOST_COMPILER=g++-15 \
  -DGMX_MPI=OFF -DGMX_HWLOC=ON -DGMX_BUILD_OWN_FFTW=OFF -DGMX_DOUBLE=OFF \
  -DCMAKE_INSTALL_PREFIX=$HOME/opt/gromacs/2026.3

cmake --build build -j"$(nproc)" && cmake --install build
source $HOME/opt/gromacs/2026.3/bin/GMXRC
```

`sha256` of the tarball:
`1094b7bbc6a3960223827114626657110b40096cdf9598a727935fc84ebf8aa0`

GROMACS publishes no checksum file next to its release tarballs, so that hash
is recorded here from the download rather than verified against upstream. The
transfer itself was over HTTPS to `ftp.gromacs.org`.

**Two flags worth copying.** `-DCMAKE_CUDA_ARCHITECTURES=86` must match your
card (`nvidia-smi --query-gpu=compute_cap --format=csv`). And the host compiler
is pinned to **g++-15**: the system default here is GCC 16.2.1, which postdates
GROMACS 2026.3 and is past the GCC 15 ceiling in CUDA 13.3's
`crt/host_config.h`.

Resulting build, from `gmx --version`:

```
GROMACS version:  2026.3
Precision:        mixed
GPU support:      CUDA
SIMD instructions: AVX2_256
CUDA targets:     86
CUDA driver:      13.30
CUDA runtime:     13.30
```

## Machine

| | |
|---|---|
| OS | Arch Linux, kernel 7.1.9 |
| CPU | 12 logical cores |
| GPU | **NVIDIA GeForce RTX 2050, 4 GB**, driver 610.57.04, compute capability 8.6 |
| CUDA | 13.3 |
| compiler | GCC 15.3.0 (CUDA host), GCC 16.2.1 (system default, unused here) |
| CMake | 4.4.2 |

A 4 GB laptop GPU is a *small* card for a 300,000-atom system. It still managed
**16.5 ns/day**, which is the floor in the README's throughput table.

## Ligand parameterisation (not needed to run)

Only required to regenerate `ligand/`:

| | |
|---|---|
| rdkit | 2026.09.1 (PyPI wheel) |
| acpype | 2026.10.15 (PyPI wheel; bundles `antechamber`, `sqm`, `parmchk2`) |
| Python | 3.13 |

`openbabel` is deliberately **not** used. `acpype` wants it only for PDB or
SMILES input; feeding it an MDL molfile instead avoids the dependency and also
passes the bond orders explicitly rather than letting them be re-derived from
interatomic distances.

## What was actually tested here

| | |
|---|---|
| apo system built | yes — 381,383 atoms |
| holo system built | yes — 306,178 atoms, 2 × estradiol, total charge 0.000 |
| minimisation | yes — converged, Fmax < 1000 in 2,632 steps |
| NVT / NPT equilibration | yes — 296.9 K, water density 991.6 kg/m³ |
| GPU-resident path | yes — non-bonded, PME, bonded and update all on the GPU |
| production + analysis | yes, both systems, but **only a 100 ps test run** |
| throughput measured | apo 12.8 ns/day, holo 16.6 ns/day |
| head separation vs the project's own script | agrees to 0.04 Å (holo), 0.18 Å (apo) |
| **a full 100 ns production run** | **no** |

One caveat on the runs in `runs/` as they stand on this machine: the apo
equilibration there was done before `refcoord_scaling = com` was added to
`npt.mdp`, so it carries the artifact that setting fixes. The delivered
`mdp/npt.mdp` has the fix and grompp now passes it with `-maxwarn 0`; anyone
starting from this folder gets the corrected version. `runs/` is not part of
the delivered folder.

The last line is the important one. This folder is a verified *pipeline*, not a
finished simulation.
