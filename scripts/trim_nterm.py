#!/usr/bin/env python3
"""Cut the N-terminal tail off every chain and cap the new N-terminus with an
acetyl group (ACE), so the cut leaves no charge the real protein does not have.

Why: in both AlphaFold models the first residues of each EdcR chain are a
low-confidence tail (pLDDT 42-56 over residues 1-14) pointing more than 4 nm
away from the core. A box is sized by the longest dimension of the solute, so
that tail alone sets it: 381,383 atoms for the full apo dimer, against roughly
80,000 without it, and the cost of every MD step scales with the atom count.

No coordinates are invented for the cap. Residues 1..n are removed as amino
acids, but the CA, C and O of residue n stay behind as the ACE (CA renamed
CH3), which is exactly the CH3-C(=O) the force field's ACE is made of, already
bonded to residue n+1 at the right geometry. Residue numbers are not changed,
so residue 56 is still residue 56 and head_groups.py selects the same helix.

Usage:  trim_nterm.py <in.pdb> <out.pdb> <n>
        residues 1..n removed, residue n+1 is the first amino acid, ACE-capped
"""
import sys


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    src, dst, cap = sys.argv[1], sys.argv[2], int(sys.argv[3])
    if cap < 1:
        sys.exit("n must be at least 1")

    out, capped, kept = [], {}, {}
    for line in open(src):
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            if rec in ("TER", "END"):
                out.append(line)
            continue
        chain, resid, name = line[21], int(line[22:26]), line[12:16].strip()
        if resid < cap:
            continue
        if resid == cap:
            rename = {"CA": " CH3", "C": " C  ", "O": " O  "}
            if name not in rename:
                continue
            line = line[:12] + rename[name] + line[16:17] + "ACE" + line[20:]
            capped.setdefault(chain, set()).add(name)
        else:
            kept[chain] = kept.get(chain, 0) + 1
        out.append(line)

    for chain in sorted(kept):
        if capped.get(chain) != {"CA", "C", "O"}:
            sys.exit("chain %s: residue %d lacks one of CA, C, O -- cannot build "
                     "the ACE cap from it" % (chain, cap))
    if not kept:
        sys.exit("%s: no residues left after the cut" % src)

    with open(dst, "w") as fh:
        fh.writelines(out)
    print("%s: residues 1-%d removed, residue %d ACE-capped, chains %s"
          % (dst, cap, cap + 1, " ".join(sorted(kept))))


if __name__ == "__main__":
    main()
