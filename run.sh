#!/bin/sh
# Quiver / EdcR molecular dynamics -- one entry point.
#
#   ./run.sh apo            build, equilibrate and run the ligand-free dimer
#   ./run.sh holo           the same with estradiol bound
#   ./run.sh both           one after the other
#   NS=10 ./run.sh apo      a 10 ns production run instead of the default 100
#   ./run.sh apo --build    stop after building the system
#   ./run.sh apo --analyse  only re-run the analysis on an existing trajectory
#
# Set GROMACS_ROOT=/path/to/gromacs if 'gmx' is not already on PATH.
set -eu
HERE="$(cd "$(dirname "$0")" && pwd)"

SYSTEMS="${1:?usage: run.sh <apo|holo|both> [--build|--analyse]}"
STAGE="${2:-all}"
[ "$SYSTEMS" = "both" ] && SYSTEMS="apo holo"

for sys in $SYSTEMS; do
    case "$STAGE" in
        --build)
            "$HERE/scripts/01_build.sh"       "$sys"
            ;;
        --analyse)
            "$HERE/scripts/05_analyse.sh"     "$sys"
            ;;
        all)
            "$HERE/scripts/01_build.sh"       "$sys"
            "$HERE/scripts/02_minimise.sh"    "$sys"
            "$HERE/scripts/03_equilibrate.sh" "$sys"
            "$HERE/scripts/04_production.sh"  "$sys"
            "$HERE/scripts/05_analyse.sh"     "$sys"
            ;;
        *)
            echo "unknown stage '$STAGE'" >&2; exit 1 ;;
    esac
done

printf '\n\033[1m==> done\033[0m\n'
for sys in $SYSTEMS; do
    printf '    %s: %s/runs/%s/md.xtc\n' "$sys" "$HERE" "$sys"
done
