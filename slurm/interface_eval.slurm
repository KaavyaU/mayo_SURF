#!/bin/bash
#SBATCH --job-name=chai_homodimer_eval
#SBATCH --partition=gen-cpu.p 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=45
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/rodata/research/brownfieldd/m355715/slurm/slurmout/%N.%j.stdout
#SBATCH --error=/rodata/research/brownfieldd/m355715/slurm/slurmout/%N.%j.stderr

# ── Step 0: Activate Python virtualenv ─────────────────────────────────────
# Activate virtual environment
source /rodata/research/brownfieldd/m355715/interface_eval/bin/activate

echo "Using python: $(python --version) from $(which python)"

# from_chai_homodimer.slurm — Homodimer evaluation pipeline via SLURM
#
# Usage (from login node):
#   sbatch from_chai_homodimer.slurm <prediction_dir> [--resume]
#
# Expects:
#   <prediction_dir>/     *.cif  *.npz   (Chai-1 predictions + scores)
#   <prediction_dir>/../pdb_files/       (design PDBs, same stem without _0 suffix)

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: sbatch $0 <prediction_dir> [--resume]"
    echo ""
    echo "  prediction_dir   Directory containing .cif and .npz files from Chai-1"
    echo "  --resume         Skip designs already in the output CSV"
    exit 1
fi

PRED_DIR=$(realpath "$1")
shift

OUTPUT_CSV="${PRED_DIR}/preds_evaluation.csv"

# ── Step 1: Organize files ─────────────────────────────────────────────────

CIF_DIR="${PRED_DIR}/cif_files"
NPZ_DIR="${PRED_DIR}/npz_files"
mkdir -p "${CIF_DIR}" "${NPZ_DIR}"

n_cif=0
n_npz=0

for f in "${PRED_DIR}"/*.cif; do
    [[ -e "$f" ]] || continue
    mv "$f" "${CIF_DIR}/"
    ((n_cif++)) || true
done

for f in "${PRED_DIR}"/*.npz; do
    [[ -e "$f" ]] || continue
    mv "$f" "${NPZ_DIR}/"
    ((n_npz++)) || true
done

# Make directory for pdb files
PDB_DIR="$(realpath "${PRED_DIR}/../pdb_files")"
mkdir -p "${PDB_DIR}"

python "/rodata/research/brownfieldd/m355715/scripts/converter.py" \
    "${CIF_DIR}" \
    --output "${PDB_DIR}"

# Count files already in subfolders (if run_eval.sh was already run once)
n_cif_existing=$(find "${CIF_DIR}" -maxdepth 1 -name '*.cif' | wc -l)
n_npz_existing=$(find "${NPZ_DIR}" -maxdepth 1 -name '*.npz' | wc -l)

echo ""
echo "============================================================"
echo " Homodimer Evaluation Pipeline"
echo "============================================================"
echo "  Prediction dir : ${PRED_DIR}"
echo "  CIF files      : ${n_cif_existing}"
echo "  NPZ files      : ${n_npz_existing}"
echo "  Output CSV     : ${OUTPUT_CSV}"
echo "============================================================"
echo ""

if [[ "${n_cif_existing}" -eq 0 ]]; then
    echo "ERROR: No .cif files found in ${CIF_DIR}"
    exit 1
fi

# ── Step 2–3: Filter + Evaluate ────────────────────────────────────────────

python "/rodata/research/brownfieldd/m355715/scripts/evaluate_homodimer.py" \
    --pred_dir "${PRED_DIR}" \
    --output_csv "${OUTPUT_CSV}" \
    "$@"

echo ""
echo "Done. Results: ${OUTPUT_CSV}"
