#!/usr/bin/env python3
"""Append HEAD_A / HEAD_B index groups -- the two recognition helices -- to an
index file, so the trajectory can be measured for the quantity this project
actually cares about: how far apart the two DNA-reading heads sit.

Why this is not simply "resid 56 to 63": pdb2gmx numbers each chain from 1, so
both chains have a residue 56 and a selection by residue number alone would
silently return both of them. The chain boundary is found here by looking for
the point where the residue index stops increasing.

The range 56-63 is not arbitrary either. It is the EdcR segment that
structurally superposes onto TetR's alpha-3 recognition helix (residues 38-45),
established in validation/af3/06_convention.md -- the project's own
secondary-structure guess of 55-60 is about 1 A off and should not be used for
anything compared against the literature.

Usage:  head_groups.py <protein.gro> <index.ndx> [first_resid] [last_resid]
"""
import sys

FIRST, LAST = 56, 63


def main():
    gro, ndx = sys.argv[1], sys.argv[2]
    first = int(sys.argv[3]) if len(sys.argv) > 3 else FIRST
    last = int(sys.argv[4]) if len(sys.argv) > 4 else LAST

    lines = open(gro).read().splitlines()
    n = int(lines[1].strip())

    # (atom number, residue id, atom name), atom numbers 1-based as GROMACS counts
    atoms = []
    for i, line in enumerate(lines[2:2 + n], start=1):
        atoms.append((i, int(line[0:5]), line[10:15].strip()))

    # Chain boundary: the first atom whose residue id is lower than the one before.
    boundary = None
    for k in range(1, len(atoms)):
        if atoms[k][1] < atoms[k - 1][1]:
            boundary = atoms[k][0]
            break
    if boundary is None:
        sys.exit("only one chain found in %s -- expected a dimer" % gro)

    heads = {"HEAD_A": [], "HEAD_B": []}
    for num, resid, name in atoms:
        if name != "CA" or not (first <= resid <= last):
            continue
        heads["HEAD_A" if num < boundary else "HEAD_B"].append(num)

    expected = last - first + 1
    for key, members in heads.items():
        if len(members) != expected:
            sys.exit("%s has %d CA atoms, expected %d -- residue numbering in "
                     "%s is not what this script assumes"
                     % (key, len(members), expected, gro))

    with open(ndx, "a") as fh:
        for key, members in heads.items():
            fh.write("[ %s ]\n" % key)
            for i in range(0, len(members), 15):
                fh.write(" ".join("%d" % m for m in members[i:i + 15]) + "\n")

    print("chain boundary at atom %d" % boundary)
    print("HEAD_A: CA of residues %d-%d, atoms %s"
          % (first, last, " ".join(map(str, heads["HEAD_A"]))))
    print("HEAD_B: CA of residues %d-%d, atoms %s"
          % (first, last, " ".join(map(str, heads["HEAD_B"]))))
    print("appended 2 groups to %s" % ndx)


if __name__ == "__main__":
    main()
