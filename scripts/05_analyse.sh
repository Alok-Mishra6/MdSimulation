#!/bin/sh
# Stage 5: measure the trajectory.
#
# Everything here is a standard GROMACS tool except the last measurement, which
# is the one specific to this project: the distance between the two
# recognition helices. That is the number the whole EdcR analysis turns on --
# AlphaFold put the heads 43.97 +/- 0.72 A apart, a real TetR on DNA measures
# 36.13 A, and the open question is whether they close when the structure is
# allowed to move.
. "$(dirname "$0")/lib.sh"
find_gmx
system_inputs "${1:?usage: 05_analyse.sh <apo|holo>}"
cd "$WORK"

[ -f md.xtc ] || die "no md.xtc -- run 04_production.sh first"
mkdir -p analysis

# md.mdp saves only the solute, so the trajectory has ~6,600 atoms while md.tpr
# describes ~300,000. Every tool below needs a reference that matches the
# trajectory, so a solute-only .tpr is made once here.
#
# The HEAD_A / HEAD_B atom numbers stay valid in this subset without
# renumbering, because the solute is the first block of atoms in the system.
say "[$1] matching reference for the solute-only trajectory"
printf 'SOLUTE\n' | $GMX convert-tpr -s md.tpr -n index.ndx -o analysis/solute.tpr
REF=analysis/solute.tpr

say "[$1] removing periodic jumps"
# Without this, a molecule stepping across the box boundary looks like a huge
# sudden motion, and every measurement below inherits the artefact.
#
# The dimer is two molecules, so keeping each one whole is not enough: -pbc mol
# can wrap the two chains into different periodic images, and every whole-dimer
# measurement then sees two halves several nm apart. (It did: Rg jumped from
# 2.3 to 4.5 nm within 50 ps.) -pbc cluster does not fix it either: the dimer
# spans most of the box, and the clustering settled on the wrong image.
#
# -pbc nojump does, because it starts from the structure given with -s: the
# equilibrated start of this run, in which the chains and ligands are still
# assembled as built, and it keeps every frame continuous with that. The dimer
# may drift out of the unit cell over time; nothing below cares.
printf 'System\n' | $GMX trjconv -s "$REF" -f md.xtc \
    -o analysis/md_fit.xtc -pbc nojump
printf 'System\n' | $GMX trjconv -s "$REF" -f analysis/md_fit.xtc \
    -o analysis/first_frame.pdb -dump 0

say "[$1] backbone RMSD against the starting structure"
printf 'Backbone\nBackbone\n' | $GMX rms -s "$REF" -f analysis/md_fit.xtc \
    -o analysis/rmsd.xvg -tu ns

say "[$1] per-residue fluctuation (RMSF)"
printf 'Backbone\n' | $GMX rmsf -s "$REF" -f analysis/md_fit.xtc \
    -o analysis/rmsf.xvg -res

say "[$1] radius of gyration"
printf 'Protein\n' | $GMX gyrate -s "$REF" -f analysis/md_fit.xtc \
    -o analysis/gyrate.xvg

say "[$1] recognition-helix separation (the project observable)"
$GMX distance -s "$REF" -f analysis/md_fit.xtc -n index.ndx \
    -select 'com of group "HEAD_A" plus com of group "HEAD_B"' \
    -oall analysis/head_separation.xvg -tu ns

if [ "$N_EST" -gt 0 ]; then
    say "[$1] has estradiol stayed in its pocket?"
    # The MINIMUM distance from each ligand to the protein, not the distance
    # between centres of mass. Estradiol sits in a buried pocket inside one
    # subunit, so its centre of mass is tens of angstrom from the centre of
    # mass of the whole dimer even when it has not moved at all -- a
    # centre-to-centre number says nothing about whether it is still bound.
    # Minimum distance does: in contact it reads ~0.17 nm, and it climbs the
    # moment the ligand leaves.
    #
    # Each copy is measured separately, because an average over both would hide
    # one of them escaping.
    $GMX select -s "$REF" -on analysis/ligand.ndx \
        -select '"PROT" not resname EST' "$(
            awk '/^(ATOM|HETATM)/ && substr($0,18,3)=="EST" {
                     r = substr($0,23,4)+0; if (!(r in seen)) { seen[r]; print r }
                 }' analysis/first_frame.pdb |
            while read -r r; do printf '"EST_%s" resname EST and resid %s; ' "$r" "$r"; done
        )"

    for grp in $(grep -o '^\[ EST_[0-9]* \]' analysis/ligand.ndx | tr -d '[] '); do
        printf '%s\nPROT\n' "$grp" | $GMX mindist -s "$REF" -f analysis/md_fit.xtc \
            -n analysis/ligand.ndx -od "analysis/mindist_$grp.xvg" -tu ns \
            >/dev/null 2>&1 || info "mindist failed for $grp"
    done
fi

say "[$1] summary"
python3 "$HERE/summarise.py" analysis "$1"
