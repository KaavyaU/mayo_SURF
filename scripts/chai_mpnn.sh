#!/usr/bin/env bash
set -euo pipefail

# USAGE: bash chai_mpnn.sh [PDB Filepath] [Output Directory]

PDB=$(realpath "$1")
OUT_DIR=$(realpath -m "$2")
#FIXED_RESIDUE_FILE=$(realpath -m "$3")
WORKING_DIR=$(dirname "$PDB")
SLURM_DIR="/rodata/research/brownfieldd/m355715/slurm"
PY_UTILS="/rodata/research/brownfieldd/m355715/scripts/chai_mpnn_utils.py"

MAX_ITERS=5
PATIENCE=2
IPTM_MARGIN=0.01
PLDDT_MARGIN=0.005
USE_MSA=true
NUM_SAMPLES=1
DIFF_BATCH_SIZE=5
N_BATCHES=1
USE_INTERFACE=true
IGNORE_STALL=true

DESIGN_MODE="binder_target" # "homodimer" or "binder_target" or "binder_dna"
BINDER_CHAIN="B"
TARGET_CHAIN="A"
TARGET_SEQ=$(python "$PY_UTILS" get-chain-seq "$PDB" "$TARGET_CHAIN")

# Thresholds
#MIN_SC=0.65
#MAX_UNSAT=1
#MAX_SURF_H=0.35
#MIN_HBONDS=2

#PDB_CLEAN="${PDB%.pdb}_clean.pdb"
#python "$PY_UTILS" strip-non-residues "$PDB" "$PDB_CLEAN" "$BINDER_CHAIN"
#PDB="$PDB_CLEAN"

echo "Cleaning previous iter_* and final_output directories"
rm -rf "$WORKING_DIR"/iter_*
rm -rf "$OUT_DIR"/final_output

mkdir -p "$OUT_DIR" "$WORKING_DIR"

stall=0
#current_struct="$PDB"
#best_struct="$PDB"
history_csv="$OUT_DIR/history.csv"
echo "iter,cif,seq,iptm,plddt,ptm,aggregate_score,length,improved,stall,interface_residues,overall_hamming,parent_hamming,overall_complex_rmsd,parent_complex_rmsd,dG,dSASA,dG_dSASA,interface_SASA_pct,interface_hydrophobicity,surface_hydrophobicity,n_unsat_hbond,unsat_hbond_pct,shape_complimentarity,helix_pct,sheet_pct,loop_pct,interface_helix_pct,interface_sheet_pct,interface_loop_pct" > "$history_csv"
WARMUP=1

echo "Creating input files"
INPUT_FASTA_NAME="${PDB%.pdb}.fasta"
INPUT_FASTA=$(python "$PY_UTILS" pdb-to-fasta "$PDB" "$INPUT_FASTA_NAME")
sbatch --wait "$SLURM_DIR/multi_chai_singles.slurm" "$WORKING_DIR" "$WORKING_DIR" "$USE_MSA" 1 >/dev/null
INPUT_NPZ="${PDB%.pdb}_0.npz"
INPUT_CIF="${PDB%.pdb}_0.cif"
INPUT_RELAXED_PDB="${INPUT_CIF%.cif}_relaxed.pdb"
INPUT_ITER=0
NA=""
read -r INPUT_IPTM INPUT_PLDDT INPUT_PTM INPUT_AGG<<< "$(python "$PY_UTILS" metrics "$INPUT_CIF" "$INPUT_NPZ")"
METRICS_RAW=$(python "$PY_UTILS" interface-eval "$INPUT_CIF" "$INPUT_RELAXED_PDB" "$BINDER_CHAIN" "$TARGET_CHAIN")
METRICS=$(printf '%s\n' "$METRICS_RAW" | tail -n 1)
read -r INPUT_D_G INPUT_D_SASA INPUT_D_G_SASA INPUT_IF_HYDRO INPUT_N_UNSAT \
        INPUT_SEQ INPUT_LEN INPUT_SC INPUT_UNSAT_HBOND_PCT INPUT_N_IFACE_RES INPUT_SURF_HYDRO INPUT_IFACE_SASA_PCT \
        INPUT_HELIX INPUT_SHEET INPUT_LOOP \
        INPUT_IFACE_HELIX INPUT_IFACE_SHEET INPUT_IFACE_LOOP <<< "$METRICS"
echo "$INPUT_ITER,$PDB,$INPUT_SEQ,$INPUT_IPTM,$INPUT_PLDDT,$INPUT_PTM,$INPUT_AGG,$INPUT_LEN,$NA,$NA,$INPUT_N_IFACE_RES,$NA,$NA,$NA,$NA,$INPUT_D_G,$INPUT_D_SASA,$INPUT_D_G_SASA,$INPUT_IFACE_SASA_PCT,$INPUT_IF_HYDRO,$INPUT_SURF_HYDRO,$INPUT_N_UNSAT,$INPUT_UNSAT_HBOND_PCT,$INPUT_SC,$INPUT_HELIX,$INPUT_SHEET,$INPUT_LOOP,$INPUT_IFACE_HELIX,$INPUT_IFACE_SHEET,$INPUT_IFACE_LOOP" >> "$history_csv"


best_iptm="${INPUT_IPTM}"
best_plddt="${INPUT_PLDDT}"
current_struct="$INPUT_CIF"
best_struct="$INPUT_CIF"

get_metric_from_csv() {
    local csv="$1"
    local cif="$2"
    local col="$3"
    awk -F',' -v cif="$cif" -v col="$col" '
    NR==1 {
        for (i=1; i<=NF; i++) if ($i == col) { idx=i; break }
        next
    }
    $1 == cif { print $idx; exit }
    ' "$csv"
}

for iter in $(seq 1 "$MAX_ITERS"); do
    ITER_DIR="$WORKING_DIR/iter_$iter"
    MPNN_DIR="$ITER_DIR/mpnn"
    CHAI_IN="$ITER_DIR/chai_in"
    CHAI_OUT="$ITER_DIR/chai_out"
    mkdir -p "$MPNN_DIR" "$CHAI_IN" "$CHAI_OUT"

    rm -f "$MPNN_DIR"/*.fa "$MPNN_DIR"/*.fasta

    echo "=== Iteration $iter: MPNN ==="
    #echo "Parsing conserved residues text file"
    #FIXED_RESIDUES=$(python "$PY_UTILS" fixed-residues "$FIXED_RESIDUE_FILE")

    echo "Started running mpnn"
    #sbatch --wait "$SLURM_DIR/mpnn.slurm" "$current_struct" "$MPNN_DIR" "$FIXED_RESIDUES" "$DIFF_BATCH_SIZE" "$N_BATCHES" >/dev/null
    sbatch --wait "$SLURM_DIR/mpnn.slurm" "$current_struct" "$MPNN_DIR" "$DIFF_BATCH_SIZE" "$N_BATCHES" "$DESIGN_MODE" >/dev/null
    echo "Finished running MPNN"

    shopt -s nullglob
    FASTA_FILES=("$MPNN_DIR"/*.fa "$MPNN_DIR"/*.fasta)
    shopt -u nullglob
    
    if [ ${#FASTA_FILES[@]} -eq 0 ]; then
	    echo "ERROR: No FASTA found in $MPNN_DIR"
	    exit 1
    fi
    
    FASTA="${FASTA_FILES[0]}"
    rm -f "$CHAI_IN"/*.fasta
    rm -f "$CHAI_OUT"/*.cif "$CHAI_OUT"/*.npz

    echo "Splitting MPNN fastas for Chai"
    python "$PY_UTILS" split-fasta "$FASTA" "$CHAI_IN" "$DESIGN_MODE" "$TARGET_SEQ" "$TARGET_CHAIN"
   
 
    echo "=== Iteration $iter: Chai ==="
    echo "Started running Chai"
    sbatch --wait "$SLURM_DIR/multi_chai_singles.slurm" "$CHAI_IN" "$CHAI_OUT" "$USE_MSA" "$NUM_SAMPLES" >/dev/null
    echo "Finished running Chai"

    
    echo "Calculating interface metrics for structures"
    METRICS_CSV="$CHAI_OUT/metrics.csv"
    FILTERED_CSV="$CHAI_OUT/filtered_ranked_metrics.csv"
    
    echo "cif,npz,relaxed_pdb,iptm,plddt,ptm,aggregate_score,dG,dSASA,dG_dSASA,interface_hydrophobicity,n_unsat_hbonds,binder_seq,binder_length,shape_complementarity,unsat_hbond_pct,n_interface_residues,surface_hydrophobicity,interface_sasa_pct,binder_helix_pct,binder_sheet_pct,binder_loop_pct,interface_helix_pct,interface_sheet_pct,interface_loop_pct" > "$METRICS_CSV"
    
    shopt -s nullglob
    CIF_FILES=("$CHAI_OUT"/*.cif)
    shopt -u nullglob
    
    for CIF in "${CIF_FILES[@]}"; do
	    NPZ="${CIF%.cif}.npz"
	    RELAXED_PDB="${CIF%.cif}_relaxed.pdb"
	    
	    # Skip if npz is missing
	    if [ ! -f "$NPZ" ]; then
		    continue
	    fi
	    
	    # Confidence metrics
	    read -r IPTM PLDDT PTM AGG <<< "$(python "$PY_UTILS" metrics "$CIF" "$NPZ")"
	    
	    # Interface metrics (your existing interface-eval)
	    METRICS_RAW=$(python "$PY_UTILS" interface-eval "$CIF" "$RELAXED_PDB" "$BINDER_CHAIN" "$TARGET_CHAIN")
	    METRICS_LINE=$(printf '%s\n' "$METRICS_RAW" | tail -n 1)
	    
	    read -r D_G D_SASA D_G_SASA IF_HYDRO N_UNSAT \
		    BINDER_SEQ BINDER_LEN SC UNSAT_HBOND_PCT N_IFACE_RES SURF_HYDRO IFACE_SASA_PCT \
		    BINDER_HELIX BINDER_SHEET BINDER_LOOP \
		    IFACE_HELIX IFACE_SHEET IFACE_LOOP <<< "$METRICS_LINE"
	    
	    # Write one row per CIF to metrics.csv
	    echo "$CIF,$NPZ,$RELAXED_PDB,$IPTM,$PLDDT,$PTM,$AGG,$D_G,$D_SASA,$D_G_SASA,$IF_HYDRO,$N_UNSAT,$BINDER_SEQ,$BINDER_LEN,$SC,$UNSAT_HBOND_PCT,$N_IFACE_RES,$SURF_HYDRO,$IFACE_SASA_PCT,$BINDER_HELIX,$BINDER_SHEET,$BINDER_LOOP,$IFACE_HELIX,$IFACE_SHEET,$IFACE_LOOP" >> "$METRICS_CSV"
    done

    echo "Ranking candidate structures"
    
    python "$PY_UTILS" rank-filtered-csv \
	    "$METRICS_CSV" "$FILTERED_CSV" "$USE_INTERFACE"
    
    # Now metrics.csv itself is filtered, ranked, and has a 'rank' column
    # Get best and second-best CIFs (assuming 'cif' column)
    BEST_CIF=$(sed -n '2p' "$FILTERED_CSV" | cut -d',' -f1)
    SECOND_CIF=$(sed -n '3p' "$FILTERED_CSV" | cut -d',' -f1)


    # If no candidates passed filters
    if [ -z "$BEST_CIF" ]; then
	    echo "No candidates passed thresholds this iteration; keeping previous best."
	    current_struct="$best_struct"
	    stall=$((stall + 1))
	    continue
    fi
    
    # Use best as current candidate
    CUR_CIF="$BEST_CIF"
    CUR_NPZ="${CUR_CIF%.cif}.npz"
    RELAXED_PDB="${CUR_CIF%.cif}_relaxed.pdb"
    
    echo "Selected candidate: $CUR_CIF"

    CUR_IPTM=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "iptm")
    CUR_PLDDT=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "plddt")
    CUR_PTM=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "ptm")
    CUR_AGG=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "aggregate_score")
    
    D_G=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "dG")
    D_SASA=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "dSASA")
    D_G_SASA=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "dG_dSASA")
    IF_HYDRO=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "interface_hydrophobicity")
    N_UNSAT=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "n_unsat_hbonds")

    BINDER_SEQ=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "binder_seq")
    BINDER_LEN=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "binder_length")
    SC=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "shape_complementarity")
    UNSAT_HBOND_PCT=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "unsat_hbond_pct")
    N_IFACE_RES=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "n_interface_residues")
    SURF_HYDRO=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "surface_hydrophobicity")
    IFACE_SASA_PCT=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "interface_sasa_pct")
    
    BINDER_HELIX=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "binder_helix_pct")
    BINDER_SHEET=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "binder_sheet_pct")
    BINDER_LOOP=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "binder_loop_pct")
    IFACE_HELIX=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "interface_helix_pct")
    IFACE_SHEET=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "interface_sheet_pct")
    IFACE_LOOP=$(get_metric_from_csv "$METRICS_CSV" "$CUR_CIF" "interface_loop_pct")

    echo "Calculating hamming distances"
    HAMMING_OVERALL=$(python "$PY_UTILS" hamming-dist "$PDB" "$CUR_CIF" "$BINDER_CHAIN")
    HAMMING_PARENT=$(python "$PY_UTILS" hamming-dist "$current_struct" "$CUR_CIF" "$BINDER_CHAIN")
    
    echo "Calculating complex RMSD"
    RMSD_OVERALL=$(python "$PY_UTILS" complex-rmsd "$PDB" "$CUR_CIF")
    RMSD_PARENT=$(python "$PY_UTILS" complex-rmsd "$current_struct" "$CUR_CIF")

    if [ -z "$best_iptm" ]; then
	    best_iptm="$CUR_IPTM"
	    best_plddt="$CUR_PLDDT"
	    best_struct="$CUR_CIF"
	    improved=true
	    stall=0
    else
	    read -r IMPROVED BEST_IPTM BEST_PLDDT <<< "$(python "$PY_UTILS" update-convergence \
		    "$best_iptm" "$best_plddt" "$CUR_CIF" "$CUR_NPZ" "$USE_INTERFACE" "$IPTM_MARGIN" "$PLDDT_MARGIN")"
	    improved=$([ "$IMPROVED" -eq 1 ] && echo true || echo false)
	    best_iptm="$BEST_IPTM"
	    best_plddt="$BEST_PLDDT"
	    
	    if [ "$improved" = true ]; then
		    best_struct="$CUR_CIF"
		    stall=0
	    else
		    stall=$((stall + 1))
		    
		    # If stalled for PATIENCE iterations, switch to second-best from this iteration
		    if [ "$IGNORE_STALL" = false ] && [ "$stall" -ge "$PATIENCE" ]; then
			    if [ -n "$SECOND_CIF" ] && [ -f "$SECOND_CIF" ]; then
				    echo "Stalled for $PATIENCE iterations; switching to 2nd-best candidate: $SECOND_CIF"
				    current_struct="$SECOND_CIF"
				    stall=0
			    else
				    echo "Stalled for $PATIENCE iterations; no second-best available; reverting to best_struct."
				    current_struct="$best_struct"
				    stall=0
			    fi
		    else
			    current_struct="$CUR_CIF"
		    fi
	    fi
    fi
    
    # Append the best candidate’s row to history.csv
    echo "Updating history.csv"
    echo "$iter,$CUR_CIF,$BINDER_SEQ,$CUR_IPTM,$CUR_PLDDT,$CUR_PTM,$CUR_AGG,$BINDER_LEN,$improved,$stall,$N_IFACE_RES,$HAMMING_OVERALL,$HAMMING_PARENT,$RMSD_OVERALL,$RMSD_PARENT,$D_G,$D_SASA,$D_G_SASA,$IFACE_SASA_PCT,$IF_HYDRO,$SURF_HYDRO,$N_UNSAT,$UNSAT_HBOND_PCT,$SC,$BINDER_HELIX,$BINDER_SHEET,$BINDER_LOOP,$IFACE_HELIX,$IFACE_SHEET,$IFACE_LOOP" >> "$history_csv"

done

FINAL_DIR="$OUT_DIR/final_output"
mkdir -p "$FINAL_DIR"

if [ ! -f "$current_struct" ]; then
    echo "ERROR: final structure not found at $current_struct"
    exit 1
fi

rm -f "$FINAL_DIR"/*.cif

cp "$current_struct" "$FINAL_DIR/"


