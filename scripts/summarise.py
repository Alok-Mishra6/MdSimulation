#!/usr/bin/env python3
"""Reduce the .xvg files from 05_analyse.sh to a few numbers and a verdict.

Deliberately plain: it reports what the trajectory did, including when the
answer is "nothing happened", and it does not dress up a short run as a result.

Usage:  summarise.py <analysis_dir> <system_name>
"""
import os
import sys

# The published reference values this project established, for context only.
AF_HEAD_MEAN = 43.97   # A, AlphaFold 3 prediction across 55 models
AF_HEAD_SD = 0.72
TETR_ON_DNA = 36.13    # A, 1QPI, a real TetR bound to its operator
TETR_INDUCED = 38.88   # A, 2TRT, induced (ligand-bound) TetR


def read_xvg(path):
    """[(x, y)] from an .xvg, ignoring comment and legend lines."""
    if not os.path.exists(path):
        return []
    out = []
    for line in open(path):
        if line.startswith(("#", "@")):
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                out.append((float(parts[0]), float(parts[1])))
            except ValueError:
                pass
    return out


def stats(values):
    if not values:
        return None
    n = len(values)
    mean = sum(values) / n
    var = sum((v - mean) ** 2 for v in values) / n if n > 1 else 0.0
    return mean, var ** 0.5, min(values), max(values)


def line(label, s, unit, fmt="%.2f"):
    if s is None:
        print("  %-26s --" % label)
        return
    mean, sd, lo, hi = s
    print(("  %-26s " + fmt + " +/- " + fmt + " %s   (range " + fmt
           + " to " + fmt + ")") % (label, mean, sd, unit, lo, hi))


def main():
    d, system = sys.argv[1], sys.argv[2]

    rmsd = read_xvg(os.path.join(d, "rmsd.xvg"))
    gyr = read_xvg(os.path.join(d, "gyrate.xvg"))
    head = read_xvg(os.path.join(d, "head_separation.xvg"))
    rmsf = read_xvg(os.path.join(d, "rmsf.xvg"))

    duration = head[-1][0] if head else (rmsd[-1][0] if rmsd else 0.0)
    print()
    print("  system                     %s" % system)
    print("  trajectory length          %.3f ns, %d frames"
          % (duration, len(head) or len(rmsd)))
    print()

    line("backbone RMSD", stats([y for _, y in rmsd]), "nm", "%.3f")
    line("radius of gyration", stats([y for _, y in gyr]), "nm", "%.3f")
    line("per-residue RMSF", stats([y for _, y in rmsf]), "nm", "%.3f")

    if head:
        # .xvg from gmx distance is in nm; this project quotes angstrom.
        ang = [y * 10.0 for _, y in head]
        s = stats(ang)
        print()
        line("head separation", s, "A")
        mean, sd, lo, hi = s
        print()
        print("  for comparison:")
        print("    AlphaFold 3 prediction   %.2f +/- %.2f A" % (AF_HEAD_MEAN, AF_HEAD_SD))
        print("    TetR induced  (2TRT)     %.2f A" % TETR_INDUCED)
        print("    TetR on DNA   (1QPI)     %.2f A" % TETR_ON_DNA)
        print()
        start = ang[0]
        drift = mean - start
        print("  started at %.2f A, mean %.2f A (%+.2f A)" % (start, mean, drift))
        closed_by = start - lo
        print("  closest approach %.2f A, i.e. %.2f A of closing" % (lo, closed_by))
        if lo <= TETR_ON_DNA:
            print("  -> reaches the operator-bound spacing at least once")
        elif lo <= TETR_INDUCED:
            print("  -> reaches the induced spacing, but not the operator-bound one")
        else:
            print("  -> stays wider than both TetR reference structures")
        if duration < 10:
            print()
            print("  CAUTION: %.3f ns is far too short to say anything about whether"
                  % duration)
            print("  the heads close. Domain motion on this scale needs 100 ns or more.")
            print("  A run this length only demonstrates that the pipeline works.")

    # Per-copy minimum distance to the protein: in contact this sits near
    # 1.7 A. A climb means that copy of the ligand is leaving its pocket.
    import glob
    ligand_files = sorted(glob.glob(os.path.join(d, "mindist_EST_*.xvg")))
    if ligand_files:
        print()
        print("  estradiol -- minimum distance to the protein:")
        for path in ligand_files:
            name = os.path.basename(path)[len("mindist_"):-len(".xvg")]
            vals = [y * 10.0 for _, y in read_xvg(path)]
            s_lig = stats(vals)
            if s_lig is None:
                continue
            mean, sd, lo, hi = s_lig
            verdict = ("in contact, has not moved" if hi < 4.0
                       else "LOOSENING -- check this" if hi < 8.0
                       else "LEFT THE POCKET")
            print("    %-8s %.2f +/- %.2f A  (max %.2f)  %s"
                  % (name, mean, sd, hi, verdict))
    print()


if __name__ == "__main__":
    main()
