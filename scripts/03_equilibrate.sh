#!/bin/sh
# Stage 3: equilibration, in two steps.
#   NVT  heats the system to 300 K at fixed volume
#   NPT  lets the box relax so the water reaches the right density
# The solute is position-restrained throughout, so the water adapts to the
# structure rather than the structure being pushed around by cold solvent.
. "$(dirname "$0")/lib.sh"
find_gmx
system_inputs "${1:?usage: 03_equilibrate.sh <apo|holo>}"
cd "$WORK"

[ -f em.gro ] || die "no em.gro -- run 02_minimise.sh first"

say "[$1] NVT, 100 ps at ${TEMPERATURE} K"
$GMX grompp -f "$MDP/nvt.mdp" -c em.gro -r em.gro -p topol.top -n index.ndx \
            -o nvt.tpr -maxwarn 0
$GMX mdrun -deffnm nvt -nb gpu -pme gpu -bonded gpu -update gpu ${GMX_MDRUN_EXTRA:-}

say "[$1] NPT, 100 ps at 1 bar"
$GMX grompp -f "$MDP/npt.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top \
            -n index.ndx -o npt.tpr -maxwarn 0
$GMX mdrun -deffnm npt -nb gpu -pme gpu -bonded gpu -update gpu ${GMX_MDRUN_EXTRA:-}

say "[$1] equilibration checks"
printf 'Temperature\n' | $GMX energy -f nvt.edr -o nvt_temperature.xvg 2>/dev/null \
    | grep -i "Temperature" || true
printf 'Density\n' | $GMX energy -f npt.edr -o npt_density.xvg 2>/dev/null \
    | grep -i "Density" || true
