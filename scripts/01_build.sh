#!/bin/sh
# Stage 1: turn a bare protein structure into a solvated, neutralised system.
#
#   pdb2gmx  -> topology + hydrogens for the chosen force field
#   editconf -> put it in a box
#   solvate  -> fill the box with water
#   genion   -> neutralise and add 0.15 M NaCl
#
# For the holo system the estradiol topology is spliced into topol.top between
# pdb2gmx and the box step, because the ligand has to be present before the
# box is sized and the water is added.
. "$(dirname "$0")/lib.sh"
find_gmx
system_inputs "${1:?usage: 01_build.sh <apo|holo>}"
SYS="$1"

# $WORK is derived from $RUNS and the system name, so check it really is the
# directory we think it is before deleting anything. A rm -rf built from
# variables deserves a guard, especially in a script meant to be handed around.
case "$WORK" in
    "$RUNS"/apo|"$RUNS"/holo) ;;
    *) die "refusing to delete '$WORK': not \$RUNS/apo or \$RUNS/holo" ;;
esac
if [ -d "$WORK" ]; then
    info "removing the previous $SYS run in $WORK"
    rm -rf "$WORK"
fi
mkdir -p "$WORK"
cd "$WORK"

case "$HMR" in yes|no) ;; *) die "HMR must be 'yes' or 'no', not '$HMR'" ;; esac
case "$NTERM_TRIM" in ''|*[!0-9]*) die "NTERM_TRIM must be a whole number" ;; esac
# Every later stage reads these back (load_settings in lib.sh).
printf 'NTERM_TRIM=%s\nHMR=%s\n' "$NTERM_TRIM" "$HMR" > settings.sh

if [ "$NTERM_TRIM" -gt 0 ]; then
    say "[$SYS] trimming residues 1-$NTERM_TRIM of each chain"
    python3 "$HERE/trim_nterm.py" "$PROTEIN" protein_trimmed.pdb "$NTERM_TRIM"
    PROTEIN=protein_trimmed.pdb
fi

say "[$SYS] 1/5  pdb2gmx: topology and hydrogens ($FORCEFIELD, $WATERMODEL)"
# AlphaFold gives heavy atoms only, so every hydrogen is built here. -ignh is
# set anyway: it makes the step behave identically if someone feeds in a
# structure that already has them.
$GMX pdb2gmx -f "$PROTEIN" -o protein.gro -p topol.top -i posre.itp \
             -ff "$FORCEFIELD" -water "$WATERMODEL" -ignh

if [ "$N_EST" -gt 0 ]; then
    say "[$SYS] 2/5  splicing in $N_EST copies of estradiol"
    cp "$LIGAND/est.itp" "$LIGAND/est_atomtypes.itp" "$LIGAND/posre_est.itp" .

    for ch in C D; do
        $GMX editconf -f "$LIGAND/est_pose_$ch.pdb" -o "est_$ch.gro"
    done

    # Coordinates: protein first, then the ligands, matching the order the
    # [ molecules ] section will declare.
    python3 "$HERE/merge_gro.py" complex.gro protein.gro est_C.gro est_D.gro
    START=complex.gro

    python3 - <<'PY'
import re

top = open("topol.top").read()

# [ atomtypes ] must precede every [ moleculetype ], so it goes directly after
# the force-field include; the molecule definition goes after the protein's.
ff = re.search(r'#include ".*forcefield\.itp"\n', top)
if not ff:
    raise SystemExit("could not find the forcefield include in topol.top")
top = (top[:ff.end()]
       + '\n; ---- estradiol (GAFF2 + AM1-BCC); see ligand/README.md ----\n'
         '#include "est_atomtypes.itp"\n'
       + top[ff.end():])

anchor = top.index("[ system ]")
top = (top[:anchor]
       + '; ---- estradiol molecule definition ----\n'
         '#include "est.itp"\n\n'
         '#ifdef POSRES\n'
         '#include "posre_est.itp"\n'
         '#endif\n\n'
       + top[anchor:])

# One EST entry per copy, appended to [ molecules ].
top = top.rstrip("\n") + "\nEST                 2\n"
open("topol.top", "w").write(top)
print("topol.top: estradiol atomtypes, molecule and restraints spliced in")
PY
else
    say "[$SYS] 2/5  no ligand (apo system)"
    START=protein.gro
fi

say "[$SYS] 3/5  box: $BOXTYPE, $BOXDIST nm to the edge"
# A dodecahedron holds ~71% of the water a cube needs for the same clearance,
# so the same physics costs ~29% less to simulate.
$GMX editconf -f "$START" -o boxed.gro -c -d "$BOXDIST" -bt "$BOXTYPE"

say "[$SYS] 4/5  solvate"
$GMX solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

say "[$SYS] 5/5  ions: neutralise + $SALT_CONC M NaCl"
# -maxwarn 1 here and nowhere else: at this point the system still has a net
# charge, so grompp warns about Ewald with a charged system. That warning is
# the entire reason the next command exists. Every other stage runs with
# -maxwarn 0 so that a real problem stops the pipeline instead of scrolling by.
$GMX grompp -f "$MDP/ions.mdp" -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1
printf 'SOL\n' | $GMX genion -s ions.tpr -o system.gro -p topol.top \
        -pname NA -nname CL -neutral -conc "$SALT_CONC"

say "[$SYS] index groups"
# SOLUTE / SOLVENT are defined by what they are NOT, so the identical
# definition works for the apo and the holo system.
$GMX select -s system.gro -on index.ndx \
    -select '"SOLUTE" not resname SOL NA CL' '"SOLVENT" resname SOL NA CL'

say "[$SYS] recognition-helix groups"
# Adds HEAD_A / HEAD_B so the trajectory can be measured for the separation of
# the two DNA-reading heads -- the observable this project is actually about.
python3 "$HERE/head_groups.py" protein.gro index.ndx

say "[$SYS] built"
$GMX check -f system.gro 2>&1 | grep -iE "atoms|box" | head -4 || true
info "system:   $WORK/system.gro"
info "topology: $WORK/topol.top"
