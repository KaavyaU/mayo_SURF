#!/usr/bin/env python3
"""
evaluate_homodimers.py — Homodimer design evaluation pipeline

Evaluates Chai-1 homodimer predictions: confidence filtering, FastRelax,
interface scoring, surface hydrophobicity, and RMSD to the design PDB.

Pipeline:
    1. Parse confidence scores (pLDDT from CIF B-factors, ipTM from .npz)
    2. Filter: mean pLDDT > 0.8  AND  ipTM > 0.7
    3. Convert passing CIFs → PDB
    4. For each passing design:
       a. Pre-relax clash check
       b. FastRelax (200 iter, constrained)
       c. Post-relax clash check
       d. Interface scoring (dG, SC, packstat, BUNs, surface hydrophobicity, ...)
       e. Secondary structure + pLDDT (on unrelaxed PDB)
       f. RMSD: predicted complex → design PDB  (CA-only)

Directory layout:
    <pred_dir>/
    ├── cif_files/       *.cif    (Chai-1 predictions, organized by run_eval.sh)
    ├── npz_files/       *.npz    (Chai-1 confidence scores)
    └── ../pdb_files/    *.pdb    (design PDBs — stem matches after stripping _0)

Usage:
    python evaluate_homodimers.py --pred_dir ./preds --output_csv ./preds/evaluation.csv
    python evaluate_homodimers.py --pred_dir ./preds --output_csv ./preds/evaluation.csv --resume
"""

import os
import sys
import re
import csv
import argparse
import traceback
import warnings
import time
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Superimposer
from Bio.PDB.Polypeptide import is_aa

from evaluate_homodimer_functions import (
    init_pyrosetta,
    pr_relax,
    score_interface,
    calc_ss_percentage,
    calculate_clash_score,
    get_binder_sequence,
)

# ── Constants ─────────────────────────────────────────────────────────────

BINDER_CHAIN  = 'A'
TARGET_CHAIN  = 'B'
N_WORKERS     = 8
PLDDT_MIN     = 0.8
IPTM_MIN      = 0.7
RELAX_ITER    = 200
CLASH_CUTOFF  = 2.4

# ── CSV schema ────────────────────────────────────────────────────────────

COLUMNS = [
    # Identity
    'pred_name', 'design_name', 'input_pdb', 'design_pdb',

    # Confidence (from npz + CIF)
    'mean_plddt', 'aggregate_score', 'ptm', 'iptm',
    'ptm_chain_A', 'ptm_chain_B',
    'iptm_A_to_B', 'iptm_B_to_A',
    'has_inter_chain_clashes',

    # RMSD
    'complex_rmsd',

    # Sequence
    'sequence', 'binder_length',

    # Clashes
    'clashes_prerelax_allatom', 'clashes_prerelax_ca_only',
    'clashes_postrelax_allatom', 'clashes_postrelax_ca_only',

    # Interface energetics (relaxed)
    'dG', 'dSASA', 'dG_dSASA_ratio', 'binder_score',
    'shape_complementarity', 'packstat',
    'n_hbonds', 'hbond_pct', 'n_unsat_hbonds', 'unsat_hbond_pct',
    'n_interface_residues', 'interface_hydrophobicity',
    'surface_hydrophobicity', 'interface_sasa_pct',

    # Secondary structure (unrelaxed)
    'binder_helix_pct', 'binder_sheet_pct', 'binder_loop_pct',
    'interface_helix_pct', 'interface_sheet_pct', 'interface_loop_pct',
    'i_plddt', 'ss_plddt',

    # Status
    'status', 'error_message',
]


# ── Naming helpers ────────────────────────────────────────────────────────

def strip_model_suffix(stem):
    """Remove trailing Chai model index: 'foo_0' → 'foo'."""
    return re.sub(r'_\d+$', '', stem)

def find_design_stem(pred_stem, design_stems):
    """Find the longest design stem that is a prefix of pred_stem."""
    matches = [s for s in design_stems if pred_stem.startswith(s + '_')]
    if not matches:
        return None
    return max(matches, key=len)

# ── Confidence extraction ─────────────────────────────────────────────────

def extract_plddt_from_cif(cif_path):
    """Return mean pLDDT (0–1) from CA B-factors in a CIF file."""
    parser = MMCIFParser(QUIET=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        structure = parser.get_structure('s', cif_path)

    bfactors = []
    for atom in structure.get_atoms():
        if atom.get_name() == 'CA':
            bfactors.append(atom.get_bfactor())

    if not bfactors:
        return None
    return round(np.mean(bfactors) / 100.0, 4)  # Chai stores pLDDT as 0–100


def parse_npz(path):
    """Return flat dict of Chai-1 confidence metrics from .npz."""
    data   = np.load(path, allow_pickle=True)
    result = {}

    for key in ['aggregate_score', 'ptm', 'iptm']:
        if key in data:
            result[key] = round(float(data[key][0]), 4)

    if 'per_chain_ptm' in data:
        for i, val in enumerate(data['per_chain_ptm'][0]):
            result[f'ptm_chain_{chr(ord("A") + i)}'] = round(float(val), 4)

    if 'per_chain_pair_iptm' in data:
        mat = data['per_chain_pair_iptm'][0]
        for i in range(mat.shape[0]):
            for j in range(mat.shape[0]):
                if i != j:
                    ci = chr(ord('A') + i)
                    cj = chr(ord('A') + j)
                    result[f'iptm_{ci}_to_{cj}'] = round(float(mat[i, j]), 4)

    if 'has_inter_chain_clashes' in data:
        result['has_inter_chain_clashes'] = bool(data['has_inter_chain_clashes'][0])

    return result


# ── CIF → PDB conversion ─────────────────────────────────────────────────

def convert_cif_to_pdb(cif_path, pdb_path):
    parser = MMCIFParser(QUIET=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        structure = parser.get_structure(Path(cif_path).stem, cif_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path)


# ── RMSD (complex CA-only) ────────────────────────────────────────────────

def compute_complex_rmsd(pred_pdb, design_pdb, chain='A'):
    """
    CA-only RMSD: align predicted chain against design chain.
    For homodimers chain A is arbitrary; both chains are the same protein.
    """
    parser = PDBParser(QUIET=True)
    ref_model = parser.get_structure('ref', design_pdb)[0]
    mob_model = parser.get_structure('mob', pred_pdb)[0]

    if chain not in ref_model or chain not in mob_model:
        return None

    ref_res = [r for r in ref_model[chain] if is_aa(r, standard=True)]
    mob_res = [r for r in mob_model[chain] if is_aa(r, standard=True)]

    n = min(len(ref_res), len(mob_res))
    if n < 3:
        return None

    ref_atoms = [r['CA'] for r in ref_res[:n] if r.has_id('CA')]
    mob_atoms = [r['CA'] for r in mob_res[:n] if r.has_id('CA')]
    n = min(len(ref_atoms), len(mob_atoms))
    if n < 3:
        return None

    sup = Superimposer()
    sup.set_atoms(ref_atoms, mob_atoms)
    return round(sup.rms, 3)


# ── CSV helpers ───────────────────────────────────────────────────────────

def init_csv(path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    if not os.path.exists(path):
        with open(path, 'w', newline='') as f:
            csv.DictWriter(f, fieldnames=COLUMNS).writeheader()


def append_row(path, row_dict):
    row_dict.pop('interface_aa_counts', None)
    row = {col: row_dict.get(col, '') for col in COLUMNS}
    with open(path, 'a', newline='') as f:
        csv.DictWriter(f, fieldnames=COLUMNS).writerow(row)


def get_completed(csv_path):
    if not os.path.exists(csv_path):
        return set()
    done = set()
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            name = row.get('pred_name', '').strip()
            if name:
                done.add(name)
    return done


def build_error_row(pred_name, design_name, input_pdb, design_pdb, error_msg):
    row = {col: '' for col in COLUMNS}
    row.update({
        'pred_name':     pred_name,
        'design_name':   design_name,
        'input_pdb':     input_pdb,
        'design_pdb':    design_pdb,
        'status':        'FAILED',
        'error_message': str(error_msg),
    })
    return row


# ── Per-design evaluation ─────────────────────────────────────────────────

def evaluate_one(
    pred_name, design_name, pred_pdb, design_pdb,
    relax_dir, chai_scores, mean_plddt,
):
    """Full evaluation pipeline for a single passing prediction."""

    relaxed_pdb = os.path.join(relax_dir, f'{pred_name}_relaxed.pdb')

    row = {
        'pred_name':    pred_name,
        'design_name':  design_name,
        'input_pdb':    pred_pdb,
        'design_pdb':   design_pdb,
        'mean_plddt':   mean_plddt,
        'status':       'OK',
        'error_message': '',
    }
    row.update(chai_scores)

    # ── RMSD ──
    try:
        row['complex_rmsd'] = compute_complex_rmsd(
            pred_pdb, design_pdb, chain=BINDER_CHAIN
        )
    except Exception as e:
        row['complex_rmsd'] = ''

    # ── Pre-relax clashes ──
    clashes_pre_aa = calculate_clash_score(pred_pdb, threshold=CLASH_CUTOFF, only_ca=False)
    clashes_pre_ca = calculate_clash_score(pred_pdb, threshold=2.5, only_ca=True)
    row['clashes_prerelax_allatom'] = clashes_pre_aa
    row['clashes_prerelax_ca_only'] = clashes_pre_ca

    if clashes_pre_ca > 0:
        row['status']        = 'FAILED_CA_CLASHES'
        row['error_message'] = f'Cα clashes pre-relax: {clashes_pre_ca}'
        return row

    if clashes_pre_aa >= 6:
        row['status']        = 'FAILED_AA_CLASHES'
        row['error_message'] = f'Sidechain clashes pre-relax: {clashes_pre_aa}'
        return row

    # ── FastRelax ──
    pr_relax(pred_pdb, relaxed_pdb, max_iter=RELAX_ITER)
    if not os.path.exists(relaxed_pdb):
        raise RuntimeError(f'FastRelax produced no output at {relaxed_pdb}')

    # ── Post-relax clashes ──
    clashes_post_aa = calculate_clash_score(relaxed_pdb, threshold=CLASH_CUTOFF, only_ca=False)
    clashes_post_ca = calculate_clash_score(relaxed_pdb, threshold=2.5, only_ca=True)
    row['clashes_postrelax_allatom'] = clashes_post_aa
    row['clashes_postrelax_ca_only'] = clashes_post_ca

    # ── Interface scoring (relaxed) ──
    iface = score_interface(relaxed_pdb, binder_chain=BINDER_CHAIN, target_chain=TARGET_CHAIN)
    row.update(iface)

    # ── Secondary structure + pLDDT (unrelaxed) ──
    ss = calc_ss_percentage(pred_pdb, binder_chain=BINDER_CHAIN)
    row.update(ss)

    # ── Sequence ──
    sequence = get_binder_sequence(pred_pdb, binder_chain=BINDER_CHAIN)
    row['sequence']      = sequence
    row['binder_length'] = len(sequence)

    return row


# ── Parallel worker plumbing ──────────────────────────────────────────────

def _worker_init(silent):
    init_pyrosetta(silent=silent)


def _worker_task(kwargs):
    try:
        row = evaluate_one(**kwargs)
        return row, None
    except Exception as e:
        return None, (kwargs['pred_name'], traceback.format_exc(), str(e))


# ── Confidence CSV (all predictions, pass + fail) ─────────────────────────

CONFIDENCE_COLUMNS = [
    'pred_name', 'design_name',
    'mean_plddt', 'aggregate_score', 'ptm', 'iptm',
    'ptm_chain_A', 'ptm_chain_B',
    'iptm_A_to_B', 'iptm_B_to_A',
    'has_inter_chain_clashes',
    'passed',
]


def write_confidence_csv(path, records):
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=CONFIDENCE_COLUMNS)
        writer.writeheader()
        for rec in records:
            writer.writerow({col: rec.get(col, '') for col in CONFIDENCE_COLUMNS})


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        prog='evaluate_homodimers.py',
        description='Evaluate Chai-1 homodimer predictions.',
    )
    parser.add_argument('--pred_dir', '-d', required=True,
                        help='Prediction directory (contains cif_files/ and npz_files/).')
    parser.add_argument('--output_csv', '-o', required=True,
                        help='Output CSV path.')
    parser.add_argument('--resume', action='store_true', default=False,
                        help='Skip predictions already in the output CSV.')
    args = parser.parse_args()

    pred_dir = os.path.realpath(args.pred_dir)
    cif_dir  = os.path.join(pred_dir, 'cif_files')
    npz_dir  = os.path.join(pred_dir, 'npz_files')
    pdb_dir  = os.path.join(pred_dir, '..', 'pdb_files')
    pdb_out  = os.path.join(pred_dir, 'pdb_files')
    relax_dir = os.path.join(pred_dir, 'relaxed')

    for d, label in [(cif_dir, 'cif_files'), (npz_dir, 'npz_files')]:
        if not os.path.isdir(d):
            print(f'ERROR: Expected {label}/ not found: {d}')
            sys.exit(1)

    if not os.path.isdir(pdb_dir):
        print(f'ERROR: Design PDB directory not found: {os.path.realpath(pdb_dir)}')
        print(f'       Expected at ../pdb_files/ relative to --pred_dir')
        sys.exit(1)

    os.makedirs(pdb_out, exist_ok=True)
    os.makedirs(relax_dir, exist_ok=True)

    # ── Collect design PDBs ───────────────────────────────────────────────

    design_pdbs = {}
    for f in os.listdir(pdb_dir):
        if f.lower().endswith('.pdb'):
            stem = Path(f).stem
            design_pdbs[stem] = os.path.join(pdb_dir, f)

    print(f'\nFound {len(design_pdbs)} design PDB(s) in {os.path.realpath(pdb_dir)}')

    # ── Step 1: Parse confidence scores ───────────────────────────────────

    print(f'\nParsing confidence scores...')

    cif_files = sorted(Path(cif_dir).glob('*.cif'))
    npz_lookup = {Path(f).stem: str(f) for f in Path(npz_dir).glob('*.npz')}

    confidence_records = []
    passing = []
    n_total = 0
    n_pass  = 0

    for cif_path in cif_files:
        raw_stem     = cif_path.stem                    # e.g. foo_0
        stripped_stem = strip_model_suffix(raw_stem)     # e.g. foo
        n_total += 1

        # Extract pLDDT from CIF
        mean_plddt = extract_plddt_from_cif(str(cif_path))

        # Extract Chai scores from NPZ
        chai_scores = {}
        if raw_stem in npz_lookup:
            try:
                chai_scores = parse_npz(npz_lookup[raw_stem])
            except Exception as e:
                print(f'  [WARN] Could not parse npz for {raw_stem}: {e}')

        iptm = chai_scores.get('iptm')

        # Filter
        plddt_ok = mean_plddt is not None and mean_plddt >= PLDDT_MIN
        iptm_ok  = iptm is not None and iptm >= IPTM_MIN
        passed   = plddt_ok and iptm_ok

        plddt_str = f'{mean_plddt:.4f}' if mean_plddt is not None else 'N/A'
        iptm_str  = f'{iptm:.4f}' if iptm is not None else 'N/A'

        record = {
            'pred_name':    stripped_stem,
            'design_name':  stripped_stem,
            'mean_plddt':   plddt_str,
            'passed':       int(passed),
        }
        record.update(chai_scores)
        confidence_records.append(record)

        if not passed:
            tag = []
            if not plddt_ok:
                tag.append(f'pLDDT={plddt_str}')
            if not iptm_ok:
                tag.append(f'ipTM={iptm_str}')
            print(f'  [FAIL] {stripped_stem}: {", ".join(tag)}')
            continue

        design_stem = find_design_stem(stripped_stem, design_pdbs.keys())
        if design_stem is None:
            print(f'  [SKIP] {stripped_stem}: no matching design PDB in ../pdb_files/')
            continue

        n_pass += 1
        print(f'  [PASS] {stripped_stem}: pLDDT={plddt_str}  ipTM={iptm_str}')

        passing.append({
            'pred_name':    stripped_stem,
            'design_name':  stripped_stem,
            'cif_path':     str(cif_path),
            'design_pdb':   design_pdbs[design_stem],
            'chai_scores':  chai_scores,
            'mean_plddt':   mean_plddt,
        })

    # Write confidence CSV (all predictions)
    confidence_csv = os.path.join(pred_dir, 'confidence_scores.csv')
    write_confidence_csv(confidence_csv, confidence_records)
    print(f'\n  {n_pass}/{n_total} passed confidence filter (pLDDT>{PLDDT_MIN}, ipTM>{IPTM_MIN})')
    print(f'  Confidence CSV: {confidence_csv}\n')

    if not passing:
        print('Nothing to evaluate. Exiting.')
        sys.exit(0)

    # ── Convert passing CIFs → PDB ───────────────────────────────────────

    print('Converting passing CIFs to PDB...')
    for entry in passing:
        pdb_path = os.path.join(pdb_out, entry['pred_name'] + '.pdb')
        if not os.path.exists(pdb_path):
            convert_cif_to_pdb(entry['cif_path'], pdb_path)
        entry['pred_pdb'] = pdb_path
    print(f'  {len(passing)} PDB files ready.\n')

    # ── Resume ────────────────────────────────────────────────────────────

    completed = get_completed(args.output_csv) if args.resume else set()
    todo = [p for p in passing if p['pred_name'] not in completed]

    if args.resume and completed:
        print(f'  Resume: {len(completed)} already done, {len(todo)} remaining.\n')

    if not todo:
        print('All passing predictions already evaluated.')
        sys.exit(0)

    # ── Worker count ──────────────────────────────────────────────────────

    n_workers = min(N_WORKERS, len(todo))

    # ── Banner ────────────────────────────────────────────────────────────

    print(f'Evaluation')
    print(f'{"─" * 56}')
    print(f'  pred_dir     : {pred_dir}')
    print(f'  Output CSV   : {args.output_csv}')
    print(f'  Chains       : binder={BINDER_CHAIN}  target={TARGET_CHAIN}')
    print(f'  Relax iters  : {RELAX_ITER}')
    print(f'  Workers      : {n_workers}')
    print(f'  To evaluate  : {len(todo)}')
    print(f'{"─" * 56}\n')

    # ── Init CSV ──────────────────────────────────────────────────────────

    init_csv(args.output_csv)

    # ── Build job kwargs ──────────────────────────────────────────────────

    job_list = [
        dict(
            pred_name    = e['pred_name'],
            design_name  = e['design_name'],
            pred_pdb     = e['pred_pdb'],
            design_pdb   = e['design_pdb'],
            relax_dir    = relax_dir,
            chai_scores  = e['chai_scores'],
            mean_plddt   = e['mean_plddt'],
        )
        for e in todo
    ]

    # ── Run ────────────────────────────────────────────────────────────────

    run_start = time.perf_counter()
    n_ok = n_failed = 0
    n_total_eval = len(todo)

    if n_workers == 1:
        print('Initializing PyRosetta...')
        init_pyrosetta(silent=True)
        print('PyRosetta ready.\n')

        for idx, entry in enumerate(job_list, start=1):
            print(f'[{idx:>4}/{n_total_eval}] {entry["pred_name"]}')
            try:
                row = evaluate_one(**entry)
                append_row(args.output_csv, row)
                status = row.get('status', 'OK')
                if status == 'OK':
                    n_ok += 1
                    print(
                        f'         ✓  dG={row.get("dG", "?")}  '
                        f'SC={row.get("shape_complementarity", "?")}  '
                        f'RMSD={row.get("complex_rmsd", "?")}'
                    )
                else:
                    n_failed += 1
                    print(f'         ✗  {status}: {row.get("error_message", "")}')
            except KeyboardInterrupt:
                print('\nInterrupted — progress saved.')
                break
            except Exception as e:
                n_failed += 1
                print(f'         ✗  EXCEPTION: {e}')
                append_row(args.output_csv, build_error_row(
                    entry['pred_name'], entry['design_name'],
                    entry['pred_pdb'], entry['design_pdb'], str(e)
                ))
    else:
        print(f'Parallel mode: {n_workers} workers\n')
        try:
            with ProcessPoolExecutor(
                max_workers=n_workers,
                initializer=_worker_init,
                initargs=(True,),
            ) as pool:
                future_to_entry = {
                    pool.submit(_worker_task, entry): entry
                    for entry in job_list
                }
                done_count = 0
                for future in as_completed(future_to_entry):
                    entry = future_to_entry[future]
                    done_count += 1
                    try:
                        row, err = future.result()
                        if err is not None:
                            _, tb, msg = err
                            n_failed += 1
                            print(f'[{done_count:>4}/{n_total_eval}] {entry["pred_name"]}  ✗  {msg}')
                            append_row(args.output_csv, build_error_row(
                                entry['pred_name'], entry['design_name'],
                                entry['pred_pdb'], entry['design_pdb'], msg
                            ))
                        else:
                            append_row(args.output_csv, row)
                            status = row.get('status', 'OK')
                            if status == 'OK':
                                n_ok += 1
                                print(
                                    f'[{done_count:>4}/{n_total_eval}] {entry["pred_name"]}  ✓  '
                                    f'dG={row.get("dG", "?")}  '
                                    f'SC={row.get("shape_complementarity", "?")}  '
                                    f'RMSD={row.get("complex_rmsd", "?")}'
                                )
                            else:
                                n_failed += 1
                                print(
                                    f'[{done_count:>4}/{n_total_eval}] {entry["pred_name"]}  ✗  '
                                    f'{status}: {row.get("error_message", "")}'
                                )
                    except Exception as e:
                        n_failed += 1
                        print(f'[{done_count:>4}/{n_total_eval}] {entry["pred_name"]}  ✗  EXCEPTION: {e}')
        except KeyboardInterrupt:
            print('\nInterrupted — results written so far are in the CSV.')

    # ── Summary ───────────────────────────────────────────────────────────

    elapsed = time.perf_counter() - run_start
    print(f'\n{"─" * 56}')
    print(f'  Done: {n_ok} passed  {n_failed} failed  ({elapsed:.0f}s)')
    print(f'  Output CSV : {args.output_csv}')
    print(f'  Relaxed    : {relax_dir}')


if __name__ == '__main__':
    main()
