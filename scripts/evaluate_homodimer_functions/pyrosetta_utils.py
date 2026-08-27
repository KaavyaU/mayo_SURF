"""
pyrosetta_utils.py — PyRosetta-based evaluation for homodimer designs.

Provides PyRosetta initialization, FastRelax, interface scoring (including
surface hydrophobicity via LayerSelector), and secondary structure analysis.
"""

import os
from collections import defaultdict

import pyrosetta as pr
from pyrosetta import rosetta
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.protocols.simple_moves import AlignChainMover
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta.rosetta.core.scoring.dssp import Dssp

from .biopython_utils import get_interface_residues


# ── Initialization ────────────────────────────────────────────────────────

def init_pyrosetta(silent=True):
    mute_flag = '-mute all' if silent else ''
    pr.init(
        f'-ignore_unrecognized_res '
        f'-ignore_zero_occupancy '
        f'{mute_flag} '
        f'-corrections::beta_nov16 true '
        f'-relax:default_repeats 1'
    )


# ── FastRelax ─────────────────────────────────────────────────────────────

def pr_relax(input_pdb, output_pdb, max_iter=200):
    """Constrained FastRelax with B-factor preservation and realignment."""

    if os.path.exists(output_pdb):
        return output_pdb

    pose       = pr.pose_from_pdb(input_pdb)
    start_pose = pose.clone()

    mm = MoveMap()
    mm.set_chi(True)
    mm.set_bb(True)
    mm.set_jump(False)

    scorefxn  = pr.get_fa_scorefxn()
    fastrelax = FastRelax()
    fastrelax.set_scorefxn(scorefxn)
    fastrelax.set_movemap(mm)
    fastrelax.max_iter(max_iter)
    fastrelax.min_type('lbfgs_armijo_nonmonotone')
    fastrelax.constrain_relax_to_start_coords(True)
    fastrelax.apply(pose)

    # Re-align to original so coordinates don't drift
    align = AlignChainMover()
    align.source_chain(0)
    align.target_chain(0)
    align.pose(start_pose)
    align.apply(pose)

    # Preserve B-factors (pLDDT) from the original prediction
    for resid in range(1, pose.total_residue() + 1):
        if pose.residue(resid).is_protein():
            bfactor = start_pose.pdb_info().bfactor(resid, 1)
            for atom_id in range(1, pose.residue(resid).natoms() + 1):
                pose.pdb_info().bfactor(resid, atom_id, bfactor)

    pose.dump_pdb(output_pdb)
    _clean_pdb(output_pdb)

    return output_pdb


# ── Interface scoring ─────────────────────────────────────────────────────

def score_interface(relaxed_pdb, binder_chain='A', target_chain='B'):
    """
    Full interface analysis on a relaxed complex PDB.

    Returns dict with: dG, dSASA, shape complementarity, packstat, H-bonds,
    buried unsatisfied H-bonds, interface/surface hydrophobicity, and more.
    """
    pose     = pr.pose_from_pdb(relaxed_pdb)
    scorefxn = pr.get_fa_scorefxn()

    # ── InterfaceAnalyzerMover ──
    iam = rosetta.protocols.analysis.InterfaceAnalyzerMover()
    iam.set_interface_jump(1) 
    iam.set_scorefunction(scorefxn)
    iam.set_compute_packstat(True)
    iam.set_compute_interface_energy(True)
    iam.set_calc_dSASA(True)
    iam.set_calc_hbond_sasaE(True)
    iam.set_compute_interface_sc(True)
    iam.set_pack_separated(True)
    iam.apply(pose)

    idata    = iam.get_all_data()
    dG       = iam.get_interface_dG()
    dSASA    = iam.get_interface_delta_sasa()
    sc       = idata.sc_value
    packstat = iam.get_interface_packstat()
    n_hbonds = idata.interface_hbonds
    dG_dSASA = idata.dG_dSASA_ratio * 100

    # ── Buried unsatisfied H-bonds ──
    try:
        buns_filter = XmlObjects.static_get_filter(
            '<BuriedUnsatHbonds '
            'report_all_heavy_atom_unsats="true" '
            'scorefxn="scorefxn" '
            'ignore_surface_res="false" '
            'use_ddG_style="true" '
            'dalphaball_sasa="0" '
            'probe_radius="1.4" '
            'burial_cutoff_apo="0.2" '
            'confidence="0" />'
        )
        n_unsat_hbonds = buns_filter.report_sm(pose)
    except Exception:
        n_unsat_hbonds = None

    # ── Interface residues + AA composition ──
    interface_residues = get_interface_residues(
        relaxed_pdb, binder_chain=binder_chain, target_chain=target_chain
    )
    n_interface = len(interface_residues)

    aa_counts = {aa: 0 for aa in 'ACDEFGHIKLMNPQRSTVWY'}
    for aa in interface_residues.values():
        if aa in aa_counts:
            aa_counts[aa] += 1

    hydrophobic_aa       = set('ACFILMPVWY')
    hydrophobic_count    = sum(aa_counts[aa] for aa in hydrophobic_aa)
    iface_hydrophobicity = round((hydrophobic_count / n_interface) * 100, 2) if n_interface > 0 else 0.0

    hbond_pct = round((n_hbonds / n_interface) * 100, 2) if n_interface > 0 else None
    unsat_hbond_pct = (
        round((n_unsat_hbonds / n_interface) * 100, 2)
        if (n_interface > 0 and n_unsat_hbonds is not None) else None
    )

    # ── Binder total energy ──
    chain_sel = ChainSelector(binder_chain)
    tem = pr.rosetta.core.simple_metrics.metrics.TotalEnergyMetric()
    tem.set_scorefunction(scorefxn)
    tem.set_residue_selector(chain_sel)
    binder_score = tem.calculate(pose)

    # ── Binder SASA fraction ──
    bsasa_metric = pr.rosetta.core.simple_metrics.metrics.SasaMetric()
    bsasa_metric.set_residue_selector(chain_sel)
    binder_sasa       = bsasa_metric.calculate(pose)
    interface_sasa_pct = round((dSASA / binder_sasa) * 100, 2) if binder_sasa > 0 else 0.0

    # ── Surface hydrophobicity (binder chain, LayerSelector) ──
    chain_poses = {
        pose.pdb_info().chain(pose.conformation().chain_begin(i)): p
        for i, p in zip(range(1, pose.num_chains() + 1), pose.split_by_chain())
    }

    surface_hydrophobicity = 0.0
    if binder_chain in chain_poses:
        binder_pose = chain_poses[binder_chain]
        layer_sel = pr.rosetta.core.select.residue_selector.LayerSelector()
        layer_sel.set_layers(pick_core=False, pick_boundary=False, pick_surface=True)
        surface_res = layer_sel.apply(binder_pose)

        exp_apol      = 0
        total_surface = 0
        for i in range(1, len(surface_res) + 1):
            if surface_res[i]:
                res = binder_pose.residue(i)
                if res.is_apolar() or res.name() in ('PHE', 'TRP', 'TYR'):
                    exp_apol += 1
                total_surface += 1

        surface_hydrophobicity = round(exp_apol / total_surface, 3) if total_surface > 0 else 0.0

    return {
        'dG':                       round(dG, 3),
        'dSASA':                    round(dSASA, 2),
        'dG_dSASA_ratio':           round(dG_dSASA, 3),
        'binder_score':             round(binder_score, 3),
        'shape_complementarity':    round(sc, 3),
        'packstat':                 round(packstat, 3),
        'n_hbonds':                 n_hbonds,
        'hbond_pct':                hbond_pct,
        'n_unsat_hbonds':           n_unsat_hbonds,
        'unsat_hbond_pct':          unsat_hbond_pct,
        'n_interface_residues':     n_interface,
        'interface_hydrophobicity': iface_hydrophobicity,
        'surface_hydrophobicity':   surface_hydrophobicity,
        'interface_sasa_pct':       interface_sasa_pct,
        'interface_aa_counts':      aa_counts,
    }


# ── Secondary structure + pLDDT ──────────────────────────────────────────

def calc_ss_percentage(input_pdb, binder_chain='A', distance_cutoff=4.0):
    """
    Secondary structure composition and pLDDT from B-factors.
    Runs DSSP on the UNRELAXED PDB to preserve prediction SS.
    """
    pose = pr.pose_from_pdb(input_pdb)

    dssp = Dssp(pose)
    dssp.insert_ss_into_pose(pose)

    # Find binder chain range
    binder_chain_num = None
    for chain_num in range(1, pose.num_chains() + 1):
        chain_id = pose.pdb_info().chain(pose.conformation().chain_begin(chain_num))
        if chain_id == binder_chain:
            binder_chain_num = chain_num
            break

    if binder_chain_num is None:
        raise ValueError(f"Binder chain '{binder_chain}' not found in {input_pdb}")

    binder_start = pose.conformation().chain_begin(binder_chain_num)
    binder_end   = pose.conformation().chain_end(binder_chain_num)

    # Interface residues for interface-specific SS
    interface_residues = get_interface_residues(
        input_pdb, binder_chain=binder_chain, distance_cutoff=distance_cutoff
    )
    interface_res_nums = set(interface_residues.keys())

    ss_counts           = defaultdict(int)
    ss_interface_counts = defaultdict(int)
    plddts_interface    = []
    plddts_ss           = []

    for pose_resid in range(binder_start, binder_end + 1):
        ss_char = pose.secstruct(pose_resid)

        if ss_char == 'H':
            ss_type = 'helix'
        elif ss_char == 'E':
            ss_type = 'sheet'
        else:
            ss_type = 'loop'

        ss_counts[ss_type] += 1

        bfactor = pose.pdb_info().bfactor(pose_resid, 1)
        if ss_type != 'loop':
            plddts_ss.append(bfactor)

        pdb_res_num = pose.pdb_info().number(pose_resid)
        if pdb_res_num in interface_res_nums:
            ss_interface_counts[ss_type] += 1
            plddts_interface.append(bfactor)

    total           = sum(ss_counts.values())
    total_interface = sum(ss_interface_counts.values())

    def pct(count, denom):
        return round((count / denom) * 100, 2) if denom > 0 else 0.0

    i_plddt  = round(sum(plddts_interface) / len(plddts_interface) / 100, 3) if plddts_interface else 0.0
    ss_plddt = round(sum(plddts_ss) / len(plddts_ss) / 100, 3) if plddts_ss else 0.0

    return {
        'binder_helix_pct':    pct(ss_counts['helix'], total),
        'binder_sheet_pct':    pct(ss_counts['sheet'], total),
        'binder_loop_pct':     pct(ss_counts['loop'],  total),
        'interface_helix_pct': pct(ss_interface_counts['helix'], total_interface),
        'interface_sheet_pct': pct(ss_interface_counts['sheet'], total_interface),
        'interface_loop_pct':  pct(ss_interface_counts['loop'],  total_interface),
        'i_plddt':             i_plddt,
        'ss_plddt':            ss_plddt,
    }


# ── Internal helpers ─────────────────────────────────────────────────────

def _clean_pdb(pdb_file):
    """Strip non-structural lines from a PDB file."""
    with open(pdb_file, 'r') as f:
        lines = [l for l in f if l.startswith(('ATOM', 'HETATM', 'MODEL', 'TER', 'END', 'LINK'))]
    with open(pdb_file, 'w') as f:
        f.writelines(lines)
