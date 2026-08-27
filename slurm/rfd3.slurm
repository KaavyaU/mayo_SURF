#!/bin/bash
#SBATCH --job-name=rfd3_design
#SBATCH --partition=gen-a100.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/rodata/research/brownfieldd/m355715/slurm/slurmout/%N.%j.stdout
#SBATCH --error=/rodata/research/brownfieldd/m355715/slurm/slurmout/%N.%j.stderr

export FOUNDRY_CHECKPOINT_DIRS=/rodata/research/brownfieldd/tools/foundry_checkpoints

JSON=$(realpath "$1")
OUT_DIR=$(realpath -m "$2")
BATCH_SIZE="${3:-2}"
N_BATCHES="${4:-10}"

if [ -z "$JSON" ] || [ -z "$OUT_DIR" ]; then
    echo "Usage: sbatch rfd3.slurm <json_file> <out_dir> [diffusion_batch_size N] [n_batches N]"
    exit 1
fi

if [ ! -f "$JSON" ]; then
    echo "ERROR: JSON file not found: $JSON"
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "Running RFD3"
echo "  JSON:        $JSON"
echo "  OUTPUT:      $OUT_DIR"
echo "  BATCH_SIZE:  $BATCH_SIZE"
echo "  N_BATCHES:   $N_BATCHES"

if rfd3 design \
    out_dir="$OUT_DIR" \
    inputs="$JSON" \
    diffusion_batch_size="$BATCH_SIZE" \
    n_batches="$N_BATCHES"; then
    echo "SUCCESS: outputs saved to $OUT_DIR"
else
    echo "ERROR: rfd3 design failed"
    exit 1
fi
