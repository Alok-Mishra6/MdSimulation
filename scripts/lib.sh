# Shared settings and helpers. Sourced by every stage script; not run directly.
set -eu

GROMACS_ROOT="${GROMACS_ROOT:-}"      # optional: a GROMACS install to source
FORCEFIELD="${FORCEFIELD:-amber99sb-ildn}"
WATERMODEL="${WATERMODEL:-tip3p}"
BOXTYPE="${BOXTYPE:-dodecahedron}"
BOXDIST="${BOXDIST:-1.2}"             # nm from solute to box edge
SALT_CONC="${SALT_CONC:-0.15}"        # mol/L NaCl, physiological
TEMPERATURE="${TEMPERATURE:-300}"

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
MDP="$ROOT/mdp"
INPUT="$ROOT/input"
LIGAND="$ROOT/ligand"
RUNS="${RUNS:-$ROOT/runs}"

say()  { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
info() { printf '    %s\n' "$*"; }
die()  { printf '\n\033[1;31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }

# GROMACS is not bundled. Find it, or tell the user exactly what to do.
find_gmx() {
    if [ -n "$GROMACS_ROOT" ] && [ -f "$GROMACS_ROOT/bin/GMXRC" ]; then
        # shellcheck disable=SC1091
        . "$GROMACS_ROOT/bin/GMXRC"
    fi
    if command -v gmx >/dev/null 2>&1; then
        GMX=gmx
    elif command -v gmx_mpi >/dev/null 2>&1; then
        GMX=gmx_mpi
    else
        die "No 'gmx' on PATH.

GROMACS is not included in this folder -- it is a 45 MB source tree and a
20-minute build, and your cluster almost certainly has it already. Try:

    module load gromacs            # on most HPC systems
    source /path/to/gromacs/bin/GMXRC
    GROMACS_ROOT=/path/to/gromacs ./run.sh ...

See README.md for building it yourself."
    fi
    export GMX
}

# One system name -> the structure files it is built from.
system_inputs() {
    case "$1" in
        apo)  PROTEIN="$INPUT/edcr_apo.pdb";           N_EST=0 ;;
        holo) PROTEIN="$INPUT/edcr_holo_protein.pdb";  N_EST=2 ;;
        *)    die "unknown system '$1' (expected 'apo' or 'holo')" ;;
    esac
    [ -f "$PROTEIN" ] || die "missing input structure: $PROTEIN"
    WORK="$RUNS/$1"
}
