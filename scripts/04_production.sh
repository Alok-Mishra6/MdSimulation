#!/bin/sh
# Stage 4: production MD. This is the run that produces the trajectory.
#
# Restraints are off and the barostat is c-rescale, which is what lets mdrun
# keep the integration on the GPU (-update gpu). On one GPU that is typically
# worth a factor of two or more.
#
#   NS=<nanoseconds>   override the length in md.mdp
#   GMX_MDRUN_EXTRA    extra mdrun flags, e.g. "-ntomp 8 -gpu_id 0"
. "$(dirname "$0")/lib.sh"
find_gmx
system_inputs "${1:?usage: 04_production.sh <apo|holo> [NS=<n>]}"
cd "$WORK"

[ -f npt.gro ] || die "no npt.gro -- run 03_equilibrate.sh first"

MDPFILE="$MDP/md.mdp"
if [ -n "${NS:-}" ]; then
    STEPS=$(awk -v ns="$NS" 'BEGIN{printf "%d", ns*500000}')   # dt = 2 fs
    sed "s/^nsteps.*/nsteps                  = $STEPS   ; ${NS} ns (set by run.sh)/" \
        "$MDP/md.mdp" > md_patched.mdp
    MDPFILE=md_patched.mdp
    info "production length overridden: ${NS} ns ($STEPS steps)"
fi

say "[$1] production"
$GMX grompp -f "$MDPFILE" -c npt.gro -t npt.cpt -p topol.top -n index.ndx \
            -o md.tpr -maxwarn 0

# If an earlier run was interrupted, resume from its checkpoint instead of
# starting over. These flags are only passed when a checkpoint actually exists:
# GROMACS refuses -append when there is nothing to append to, which would turn
# a fresh run into a hard error.
RESUME=""
if [ -f md.cpt ]; then
    info "found md.cpt -- resuming from the last checkpoint"
    RESUME="-cpi md.cpt -append"
fi
# shellcheck disable=SC2086
$GMX mdrun -deffnm md -nb gpu -pme gpu -bonded gpu -update gpu \
           $RESUME ${GMX_MDRUN_EXTRA:-}

say "[$1] trajectory written"
info "$WORK/md.xtc"
$GMX check -f md.xtc 2>&1 | grep -iE "^Reading|Coords|Last frame" | head -4 || true
