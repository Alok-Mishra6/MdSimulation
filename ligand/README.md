# 17β-estradiol parameters — how these were made

The protein force field (amber99sb-ildn) ships with GROMACS. Estradiol does
not, so its parameters were generated once, here, and are shipped ready to use.
**Nothing in this directory needs to be regenerated to run the simulation.**

| file | what it is |
|---|---|
| `est_atomtypes.itp` | the 9 GAFF2 atom types the ligand uses |
| `est.itp` | the molecule: atoms, charges, bonds, angles, dihedrals, impropers |
| `posre_est.itp` | position restraints, used during equilibration |
| `est_pose_C.pdb`, `est_pose_D.pdb` | the two ligand copies, with hydrogens, in the predicted pose |
| `est_pose_C.mdl`, `est_pose_D.mdl` | the same, as MDL molfiles (what acpype consumes) |
| `EST_ccd.cif` | the PDB Chemical Component Dictionary entry for EST |
| `est_gaff2.frcmod` | the parameters parmchk2 filled in by analogy |
| `est_am1bcc_sqm.out` | the AM1 calculation the charges came from |
| `build_est.py` | regenerates the hydrogen-added poses from the CCD + prediction |

`est_atomtypes.itp` **must be included before any `[ moleculetype ]`** — that
is a GROMACS requirement, and it is why the acpype output was split in two.
`scripts/01_build.sh` splices both into `topol.top` in the right places.

## The method

| | |
|---|---|
| atom types | GAFF2 |
| charges | AM1-BCC, via `antechamber` and `sqm` |
| tool | `acpype` 2026.10.15 (bundles the AmberTools binaries) |
| net charge | **0** — estradiol is neutral; both hydroxyls stay protonated at pH 7 |

GAFF2 + AM1-BCC is the standard partner for an AMBER protein force field, which
is the other reason amber99sb-ildn was chosen for the protein.

## Two things that were checked rather than assumed

**1. The chemistry was not guessed from the coordinates.** AlphaFold returns
heavy atoms with no bonds. Rather than let a tool infer bond orders from
interatomic distances, `build_est.py` reads them from the **PDB Chemical
Component Dictionary entry for EST**, which is the authoritative definition.
Only the coordinates come from the prediction. The CCD's 20 heavy-atom names
match the AlphaFold output exactly, in the same order, so the mapping is
unambiguous.

**2. The stereochemistry was verified.** Estradiol has five stereocentres, and
17β-estradiol differs from its 17α epimer at one of them. After hydrogens are
added, `build_est.py` computes the InChIKey of the result and compares it to
the CCD's own:

```
InChIKey      : VOXZDWNPVJITMN-ZBRFXRBCSA-N
CCD reference : VOXZDWNPVJITMN-ZBRFXRBCSA-N   -> match
```

Both copies pass. The script **exits with an error** rather than continuing if
they do not, and it distinguishes a wrong stereoisomer from a wrong molecule.
Stereochemistry is read back off the predicted 3D coordinates, so this is a
genuine check on the prediction, not a restatement of an assumption.

The hydrogens are placed geometrically on the predicted heavy-atom frame and
then relaxed with MMFF94 **with every heavy atom held fixed**, so the predicted
pose is preserved exactly and only the hydrogen positions are optimised.

## What was substituted by analogy — the honest part

`parmchk2` filled in **40 parameters** that GAFF2 does not list explicitly.
Most carry a penalty score of 0.0, meaning an exact analogue was found. Three
do not:

| parameter | penalty | what it is |
|---|---|---|
| `ca-ca-ca-oh` | 41.9 | improper keeping the phenol OH coplanar with the ring |
| `c6-ca-ca-ca` | 41.9 | improper at the aromatic/aliphatic ring fusion |
| `ca-ca-ca-ha` | 6.0 | improper keeping an aromatic CH coplanar |

All three are **impropers whose only job is to keep an aromatic ring flat**, and
all three fall back to GAFF's general aromatic planarity term (10.5 kcal/mol,
180°, periodicity 2). A high penalty on a general aromatic improper is routine
for any substituted aromatic in GAFF and does not mean the value is wrong — but
it is the loosest part of this parameter set, and anyone pushing on the phenol's
behaviour specifically should consider CGenFF or QM-derived parameters instead.

`acpype` also warns `No Openbabel python module, no chiral groups`, meaning it
did not add explicit chirality restraints. That is harmless here: every
stereocentre in estradiol sits in a fused ring, so inverting one would require
breaking a covalent bond, which classical MD cannot do.

## Regenerating, if you want to

Needs `rdkit` and `acpype` (neither is needed to *run* the simulation):

```sh
python -m venv .venv
.venv/bin/pip install rdkit acpype

# hydrogen-added poses, with the stereochemistry check.
# The input is the unmodified AlphaFold output, shipped as
# ../input/edcr_holo_with_est_ORIGINAL.pdb -- it still has the HETATM records,
# which ../input/edcr_holo_protein.pdb does not.
.venv/bin/python build_est.py EST_ccd.cif \
    ../input/edcr_holo_with_est_ORIGINAL.pdb C est_pose_C.pdb
.venv/bin/python build_est.py EST_ccd.cif \
    ../input/edcr_holo_with_est_ORIGINAL.pdb D est_pose_D.pdb

# parameters (acpype takes .mdl or .mol2, never .pdb)
.venv/bin/acpype -i est_pose_C.mdl -b EST -c bcc -a gaff2 -n 0
```

Then split `EST.acpype/EST_GMX.itp` at `[ moleculetype ]` into
`est_atomtypes.itp` and `est.itp`, and rename the residue `MOL` to `EST`.

AM1-BCC charges depend slightly on conformation, so regenerating from a
different pose gives slightly different charges. The shipped set came from the
chain C pose.
