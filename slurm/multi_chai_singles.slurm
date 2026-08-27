#!/bin/bash
#SBATCH --job-name=chai_batch
#SBATCH --partition=gen-a100.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --output=/rodata/research/brownfieldd/m355715/slurm/slurmout/%N.%j.stdout
#SBATCH --error=/rodata/research/brownfieldd/m355715/slurm/slurmout/%N.%j.stderr

unset PYTHONHOME
export CHAI_DOWNLOADS_DIR=/rodata/research/brownfieldd/tools/chai_downloads

INPUT_DIR=$(realpath "$1")
OUTPUT_DIR=$(realpath -m "$2")
USE_MSA="${3:-true}"
NUM_SAMPLES="${4:-1}"

if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: sbatch chai_batch.slurm <input_dir> <output_dir> [use_msa true|false] [num_samples N]"
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

ARGS=()
if [ "$USE_MSA" = "true" ]; then
    ARGS+=(--use-msa-server)
fi
ARGS+=(--no-use-templates-server)
ARGS+=(--num-diffn-samples "$NUM_SAMPLES")

SUCCEEDED=0
FAILED=0
FAILED_NAMES=()

for FASTA in "$INPUT_DIR"/*.fasta; do
    [ -f "$FASTA" ] || { echo "WARNING: No .fasta files found in $INPUT_DIR"; break; }

    STEM=$(basename "$FASTA" .fasta)
    TMP_DIR="$OUTPUT_DIR/tmp_${STEM}"
    mkdir -p "$TMP_DIR"

    echo "----------------------------------------"
    echo "Folding: $STEM"
    echo "  MSA:     $USE_MSA"
    echo "  SAMPLES: $NUM_SAMPLES"

    if chai-lab fold "${ARGS[@]}" "$FASTA" "$TMP_DIR"; then

        for SRC in "$TMP_DIR"/pred.model_idx_*.cif; do
            [ -f "$SRC" ] || continue
            IDX=$(echo "$SRC" | grep -o 'idx_[0-9]*' | grep -o '[0-9]*')
            mv "$SRC" "$OUTPUT_DIR/${STEM}_${IDX}.cif"
        done

        for SRC in "$TMP_DIR"/scores.model_idx_*.npz; do
            [ -f "$SRC" ] || continue
            IDX=$(echo "$SRC" | grep -o 'idx_[0-9]*' | grep -o '[0-9]*')
            mv "$SRC" "$OUTPUT_DIR/${STEM}_${IDX}.npz"
        done

        if [ -d "$TMP_DIR/msa" ]; then
            mv "$TMP_DIR/msa" "$OUTPUT_DIR/${STEM}_msa"
        fi

        rm -rf "$TMP_DIR"
        echo "  SUCCESS: outputs saved as ${STEM}_*.cif / *.npz"
        ((SUCCEEDED++))

    else
        echo "  ERROR: chai-lab failed on $STEM — skipping"
        rm -rf "$TMP_DIR"
        FAILED_NAMES+=("$STEM")
        ((FAILED++))
    fi
done

echo "========================================"
echo "BATCH COMPLETE"
echo "  Succeeded: $SUCCEEDED"
echo "  Failed:    $FAILED"
if [ ${#FAILED_NAMES[@]} -gt 0 ]; then
    echo "  Failed files:"
    for NAME in "${FAILED_NAMES[@]}"; do
        echo "    - $NAME"
    done
fi
echo "========================================"
