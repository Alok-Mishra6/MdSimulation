#!/usr/bin/env python3
"""Cut a GROMACS multi-model PDB of the EdcR dimer down to the atoms of
Quiver's validation/edcr/nt_core.pdb: CA of residues 23-201, chain A then B.

trjconv leaves the chain column blank, so the chain is assigned by order:
chain A until the residue number drops, chain B after it, per model.

Usage:  match_nt_core.py <in.pdb> <out.pdb> [first] [last]
"""
import sys

src, dst = sys.argv[1], sys.argv[2]
first = int(sys.argv[3]) if len(sys.argv) > 3 else 23
last = int(sys.argv[4]) if len(sys.argv) > 4 else 201

models = 0
per_model = None
with open(src) as fin, open(dst, "w") as fout:
    chain, prev, kept, serial = "A", None, 0, 0
    for line in fin:
        if line.startswith("MODEL"):
            chain, prev, kept, serial = "A", None, 0, 0
            fout.write(line)
        elif line.startswith("ENDMDL"):
            if chain != "B":
                sys.exit("model %d: only one chain found" % (models + 1))
            if per_model is None:
                per_model = kept
            elif kept != per_model:
                sys.exit("model %d has %d CA, model 1 had %d" % (models + 1, kept, per_model))
            models += 1
            fout.write(line)
        elif line.startswith("ATOM") and line[12:16] == " CA ":
            resid = int(line[22:26])
            if prev is not None and resid < prev:
                chain = "B"
            prev = resid
            if first <= resid <= last:
                serial += 1
                kept += 1
                fout.write("%s%5d%s%s%s" % (line[:6], serial, line[11:21], chain, line[22:]))
    fout.write("END\n")
print("%s: %d models, %d CA each (residues %d-%d, chains A and B)" % (dst, models, per_model, first, last))
