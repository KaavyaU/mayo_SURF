#!/bin/bash
#SBATCH --job-name=mpnn
#SBATCH --partition=gen-a100.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/rodata/research/brownfieldd/m355715/slurm/slurmout/%N.%j.stdout
#SBATCH --error=/rodata/research/brownfieldd/m355715/slurm/slurmout/%N.%j.stderr

export FOUNDRY_CHECKPOINT_DIRS=/rodata/research/brownfieldd/tools/foundry_checkpoints/ligandmpnn_v_32_010_25.pt

PDB=$(realpath "$1")
OUT_DIR=$(realpath -m "$2")
#FIXED_RESIDES="${3:""}"
BATCH_SIZE="${3:-2}"
N_BATCHES="${4:-10}"
DESIGN_MODE="${5:-binder_target}"
OMIT='["CYS", "MET", "LYS"]'


if [ -z "$PDB" ] || [ -z "$OUT_DIR" ]; then
    echo "Usage: sbatch mpnn.slurm <pdb_file> <out_dir> [conserved_residues] [diffusion_batch_size N] [n_batches N] [design_mode]"
    exit 1
fi

if [ ! -f "$PDB" ]; then
    echo "ERROR: PDB file not found: $PDB"
    exit 1
fi


mkdir -p "$OUT_DIR"

echo "Running MPNN"
echo "  PDB:         $PDB"
echo "  OUTPUT:      $OUT_DIR"
echo "  BATCH_SIZE:  $BATCH_SIZE"
echo "  N_BATCHES:   $N_BATCHES"

# --fixed_residues "$FIXED_RESIDUES"
MPNN_ARGS=(
    --structure_path "$PDB"
    --out_directory "$OUT_DIR"
    --checkpoint_path "$FOUNDRY_CHECKPOINT_DIRS"
    --model_type ligand_mpnn
    --is_legacy_weights True
    --batch_size "$BATCH_SIZE"
    --number_of_batches "$N_BATCHES"
    --write_structures False
    --omit "$OMIT"
    --temperature 0.2
    --seed 111
)

if [ "$DESIGN_MODE" = "homodimer" ]; then
    MPNN_ARGS+=(--homo_oligomer '[["A","B"]]')
elif [ "$DESIGN_MODE" = "binder_target" ] || [ "$DESIGN_MODE" = "binder_dna" ]; then
    MPNN_ARGS+=(--designed_chains '["B"]')
else
    echo "ERROR: Unknown DESIGN_MODE: $DESIGN_MODE (use 'homodimer', 'binder_target', or 'binder_dna')"
    exit 1
fi

if mpnn "${MPNN_ARGS[@]}" ; then
    echo "SUCCESS: outputs saved to $OUT_DIR"
else
    echo "ERROR: MPNN design failed"
    exit 1
fi


