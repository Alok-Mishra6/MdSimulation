#!/bin/sh
# Stage 2: energy minimisation. Fixes the overlaps that solvate and genion
# leave behind. If Fmax does not come down here, nothing later will work.
. "$(dirname "$0")/lib.sh"
find_gmx
system_inputs "${1:?usage: 02_minimise.sh <apo|holo>}"
cd "$WORK"

say "[$1] minimising"
$GMX grompp -f "$MDP/em.mdp" -c system.gro -r system.gro -p topol.top -o em.tpr -maxwarn 0
# shellcheck disable=SC2086
$GMX mdrun -deffnm em $EM_OFFLOAD ${GMX_MDRUN_EXTRA:-}

say "[$1] result"
printf 'Potential\nMax-Force\n' | $GMX energy -f em.edr -o em_potential.xvg 2>/dev/null \
    | grep -iE "Potential|Max-Force" || true
