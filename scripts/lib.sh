# Shared settings and helpers. Sourced by every stage script; not run directly.
set -eu

GROMACS_ROOT="${GROMACS_ROOT:-}"      # optional: a GROMACS install to source
FORCEFIELD="${FORCEFIELD:-amber99sb-ildn}"
WATERMODEL="${WATERMODEL:-tip3p}"
BOXTYPE="${BOXTYPE:-dodecahedron}"
BOXDIST="${BOXDIST:-1.2}"             # nm from solute to box edge
SALT_CONC="${SALT_CONC:-0.15}"        # mol/L NaCl, physiological
TEMPERATURE="${TEMPERATURE:-300}"

# Two ways to make a run several times cheaper; both off by default. They are
# chosen at build time and recorded in the run directory (see load_settings),
# so every later stage of that run uses the same values.
NTERM_TRIM="${NTERM_TRIM:-0}"         # drop residues 1..N of each chain, ACE-cap the rest
HMR="${HMR:-no}"                      # yes: hydrogen mass repartitioning, 4 fs production

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

    # Where mdrun does its work. A GPU build gets the GPU-resident path; a
    # CPU-only build gets no offload flags at all, because asking it for
    # '-nb gpu' is a hard error rather than a fallback. USE_GPU=no forces the
    # CPU path on a GPU build, e.g. on a node without a card.
    USE_GPU="${USE_GPU:-auto}"
    if [ "$USE_GPU" = auto ]; then
        if $GMX --version 2>/dev/null | grep -qiE '^GPU support:[[:space:]]+disabled'; then
            USE_GPU=no
        else
            USE_GPU=yes
        fi
    fi
    if [ "$USE_GPU" = yes ]; then
        EM_OFFLOAD="-nb gpu"
        MD_OFFLOAD="-nb gpu -pme gpu -bonded gpu -update gpu"
    else
        EM_OFFLOAD=""
        MD_OFFLOAD=""
    fi
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

# The settings a run was built with, written by 01_build.sh. They override the
# environment: a run built without HMR must not have it switched on halfway.
load_settings() {
    [ -f "$WORK/settings.sh" ] || die "no $WORK/settings.sh -- run 01_build.sh first"
    # shellcheck disable=SC1091
    . "$WORK/settings.sh"
    info "built with NTERM_TRIM=$NTERM_TRIM HMR=$HMR"
}

# The .mdp a stage actually runs: the shipped one with the run's settings
# applied, written into the run directory so what was run stays on record.
#   stage_mdp <nvt|npt|md>    prints the path of the file it wrote
stage_mdp() {
    out="$WORK/$1_run.mdp"
    cp "$MDP/$1.mdp" "$out"
    if [ "$HMR" = yes ]; then
        # Each hydrogen is made 3x heavier, the mass taken from the atom it is
        # bonded to. That slows the fastest motions enough for a 4 fs step;
        # equilibrium properties do not depend on the masses. Equilibration
        # keeps 2 fs, so only production takes the larger step.
        printf '\n; ---- HMR=yes (scripts/lib.sh) ----\nmass-repartition-factor = 3\n' >> "$out"
        if [ "$1" = md ]; then
            # Step counts halve so the run length and the 10 ps output
            # interval stay what md.mdp says they are.
            awk '$1 == "dt" { printf "%-24s= 0.004\n", $1; next }
                 $1 ~ /^(nsteps|nstxout-compressed|nstlog|nstenergy)$/ {
                     printf "%-24s= %d\n", $1, $3 / 2; next }
                 { print }' "$out" > "$out.tmp" && mv "$out.tmp" "$out"
        fi
    fi
    printf '%s\n' "$out"
}
