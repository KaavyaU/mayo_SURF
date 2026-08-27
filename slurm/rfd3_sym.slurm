#!/bin/bash
#SBATCH --job-name=rfd3_sym
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
BATCH_SIZE="${3:-1}"
N_BATCHES="${4:-20}"

if [ -z "$JSON" ] || [ -z "$OUT_DIR" ]; then
    echo "Usage: sbatch rfd3_sym.slurm <json_file> <out_dir> [diffusion_batch_size 1] [n_batches 20]"
    exit 1
fi

if [ ! -f "$JSON" ]; then
    echo "ERROR: JSON file not found: $JSON"
    exit 1
fi

mkdir -p "$OUT_DIR"

if rfd3 design \
    out_dir="$OUT_DIR" \
    inputs="$JSON" \
    inference_sampler.kind=symmetry \
    diffusion_batch_size="$BATCH_SIZE" \
    n_batches="$N_BATCHES"; then
    echo "SUCCESS: outputs saved to $OUT_DIR"
else
    echo "ERROR: rfd3 design failed"
    exit 1
fi

