#!/usr/bin/env python3
"""Concatenate GROMACS .gro files into one, keeping the last file's box.

Written out rather than done with shell text processing because the .gro format
is column-positional: the atom count on line 2 and the box vector on the last
line both have to be rewritten, and a stray space in the coordinate columns is
a silently wrong structure rather than an error.

Usage:  merge_gro.py <out.gro> <in1.gro> <in2.gro> [...]
"""
import sys


def read(path):
    lines = open(path).read().splitlines()
    if len(lines) < 3:
        sys.exit("%s: too short to be a .gro file" % path)
    n = int(lines[1].strip())
    atoms = lines[2:2 + n]
    if len(atoms) != n:
        sys.exit("%s: header says %d atoms, found %d" % (path, n, len(atoms)))
    return atoms, lines[2 + n]


def main():
    out_path, in_paths = sys.argv[1], sys.argv[2:]
    if not in_paths:
        sys.exit(__doc__)
    atoms, box = [], None
    for p in in_paths:
        a, box = read(p)
        atoms.extend(a)

    # Residue and atom numbering is rebuilt so the file is self-consistent;
    # GROMACS matches coordinates to the topology by ORDER, not by these
    # numbers, but a wrong count here is a hard error later and a confusing one.
    renumbered = []
    for i, line in enumerate(atoms, start=1):
        renumbered.append("%s%5d%s" % (line[:15], i % 100000, line[20:]))

    with open(out_path, "w") as fh:
        fh.write("EdcR system, merged by merge_gro.py\n")
        fh.write("%5d\n" % len(renumbered))
        fh.write("\n".join(renumbered) + "\n")
        fh.write(box + "\n")
    print("merged %d files -> %s (%d atoms)"
          % (len(in_paths), out_path, len(renumbered)))


if __name__ == "__main__":
    main()
