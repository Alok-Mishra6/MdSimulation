# Where these structures came from

Both are **AlphaFold 3 server predictions**, taken unmodified except for
stripping everything that is not protein. Neither has been energy-minimised,
refined, or altered by hand, and neither is one of this project's own derived
or engineered models.

## `edcr_apo.pdb` — the ligand-free dimer

| | |
|---|---|
| AlphaFold 3 job | `fold_edcr_alpha3_seed37` |
| model | 0 |
| source file | `validation/af3/best/EdcR_CHOSEN_apo_seed37_m0.pdb` |
| sha256 of source | `5fb83824201c20a436c95e00c9256e7cc09f7e4d9d601c896352c58b8fdd4193` |
| sha256 as delivered | `078d25ee5c3cab14e0f28c78c6f005d7a5bc76376178cf68353ef7667925fa47` |

**This is the project's chosen reference structure**, selected by the project
owner on 2026-10-16 (`validation/af3/best/README.md`). Its elastic-network
modes are in `validation/af3/modes/seed37_m0_anm15.txt` and its full analysis
in `validation/af3/modes/seed37_m0_nma.md`, so MD results can be compared
directly against the normal-mode predictions for the same coordinates.

## `edcr_holo_protein.pdb` — the estradiol-bound dimer

| | |
|---|---|
| AlphaFold 3 job | `fold_edcr_est_seed39` |
| model | 2 |
| source file | `validation/af3/best/EdcR_dimer_EST_best.pdb` |
| sha256 of source | `d3676646680a7d05be920755767f18c5ea8cb67661444578cbe6181e7daf4a55` |
| sha256 as delivered | `98523556586fbaede42b5c5da8153bac9263d211e90d0d645dc93b859b65d11d` |

Ranked **1st of 87 structures** on combined confidence, C2 symmetry,
typicality and backbone bond quality. Its estradiol sits in the majority pose
(0.16 Å from the consensus of 27 of 30 models; the three outliers flip the
steroid end for end and are exactly the three whose ligand ipTM drops to 0.62).
Ligand ipTM 0.75, no clashes.

The two estradiol copies were removed from this file and are supplied
separately, with hydrogens, as `../ligand/est_pose_C.pdb` and
`../ligand/est_pose_D.pdb` — their coordinates are untouched from this
prediction.

`edcr_holo_with_est_ORIGINAL.pdb` is that same AlphaFold output **unmodified,
with the two estradiol molecules still in it as HETATM records**. The pipeline
does not read it: it is here so the ligand parameters can be regenerated
without reaching back into the main repository, and as the reference anyone
checking provenance should diff against. Its sha256 is the "sha256 of source"
value in the table above.

## Both files

- Protein only: `ATOM` records, chains A and B, residues 1–208, 1605 atoms per
  chain, 3210 total.
- **Heavy atoms only — no hydrogens.** AlphaFold does not predict them;
  `pdb2gmx` builds them during stage 1.
- No alternate locations, no insertion codes, no partial occupancies.
- **Sequence-identical to each other** (verified by hashing the ordered residue
  list of chain A), with identical numbering. No alignment is needed to compare
  the two runs.

## What was deliberately not used

`validation/af3/best/` also holds three **closed** models, built by this project
by pushing the apo structure along its head-closing normal mode to the
operator-bound spacing. They are not used here, for two reasons: they are
**Cα-only backbone hypotheses** with no side chains, which would need
rebuilding before any force field could touch them; and they are this project's
own constructions rather than predictions, so starting MD from one would be
simulating our own hypothesis instead of testing it.

The two files marked `control` in that directory are closures of
estradiol-bound structures and are, biologically, the wrong state — see the
warning at the top of `validation/af3/best/README.md`.
