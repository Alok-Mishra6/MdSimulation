#!/bin/sh
# Check the machine before committing to a multi-day run. Prints what it finds
# and says plainly whether the GPU path is available.
. "$(dirname "$0")/lib.sh"

say "GROMACS"
find_gmx
$GMX --version 2>/dev/null | grep -iE "^(GROMACS version|Precision|GPU support|SIMD instructions|CUDA runtime)" \
    | sed 's/^/    /'

say "GPU"
if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=name,memory.total,driver_version,compute_cap \
               --format=csv,noheader | sed 's/^/    /'
else
    info "no nvidia-smi found -- mdrun will fall back to the CPU, which for a"
    info "300,000-atom system means weeks rather than days"
fi

say "CPU and memory"
info "logical cores: $( (nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo '?') )"
info "memory:        $(awk '/MemTotal/{printf "%.1f GB", $2/1048576}' /proc/meminfo 2>/dev/null || echo '?')"
info "free disk:     $(df -h "$ROOT" | awk 'NR==2{print $4}')"

say "inputs"
for f in "$INPUT/edcr_apo.pdb" "$INPUT/edcr_holo_protein.pdb" \
         "$LIGAND/est.itp" "$LIGAND/est_atomtypes.itp" \
         "$LIGAND/est_pose_C.pdb" "$LIGAND/est_pose_D.pdb"; do
    if [ -f "$f" ]; then info "ok      $(basename "$f")"
    else                 info "MISSING $(basename "$f")"; fi
done

say "force field"
if $GMX pdb2gmx -h 2>&1 | grep -q .; then
    if [ -d "$($GMX --version 2>/dev/null | awk '/Data prefix/{print $3}')/share/gromacs/top/$FORCEFIELD.ff" ]; then
        info "ok      $FORCEFIELD ships with this GROMACS"
    else
        info "check   $FORCEFIELD -- not found in the GROMACS share directory"
    fi
fi

say "ready"
info "disk needed: roughly 400 MB per 100 ns per system (solute-only trajectory)"
