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
load_settings

[ -f npt.gro ] || die "no npt.gro -- run 03_equilibrate.sh first"

MDPFILE=$(stage_mdp md)
if [ -n "${NS:-}" ]; then
    DT=$(awk '$1 == "dt" {print $3}' "$MDPFILE")             # ps
    STEPS=$(awk -v ns="$NS" -v dt="$DT" 'BEGIN{printf "%d", ns * 1000 / dt + 0.5}')
    sed "s/^nsteps.*/nsteps                  = $STEPS   ; ${NS} ns (set by run.sh)/" \
        "$MDPFILE" > "$MDPFILE.tmp" && mv "$MDPFILE.tmp" "$MDPFILE"
    info "production length overridden: ${NS} ns ($STEPS steps of $DT ps)"
fi

if [ "$HMR" = yes ]; then
    # At 4 fs the methionine S-CH3 methyls (atom CE) now and then rotate more
    # than 30 degrees in one step, which LINCS reports as a warning: ~740 an
    # hour on the EdcR dimer, every one of them on a Met CE hydrogen, with
    # constraint deviations under 1.5 %. mdrun aborts after 1000 warnings per
    # process, which here is every ~80 minutes. A real blow-up produces
    # thousands a second, so a limit of 100000 still catches one.
    GMX_MAXCONSTRWARN="${GMX_MAXCONSTRWARN:-100000}"
    export GMX_MAXCONSTRWARN
    info "HMR: LINCS warning limit raised to $GMX_MAXCONSTRWARN (Met CE methyls)"
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
$GMX mdrun -deffnm md $MD_OFFLOAD \
           $RESUME ${GMX_MDRUN_EXTRA:-}

say "[$1] trajectory written"
info "$WORK/md.xtc"
$GMX check -f md.xtc 2>&1 | grep -iE "^Reading|Coords|Last frame" | head -4 || true
