import warnings
import csv
import sys
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Superimposer, Select, NeighborSearch
from Bio.PDB.Polypeptide import is_aa, is_nucleic
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path

from evaluate_homodimer_functions import (
    init_pyrosetta,
    pr_relax,
    score_interface,
    calc_ss_percentage,
    calculate_clash_score,
    get_binder_sequence,
)


def is_pdb(path):
    return str(path).lower().endswith(".pdb")


def is_cif(path):
    # handles .cif and .mmcif
    p = str(path).lower()
    return p.endswith(".cif") or p.endswith(".mmcif")


class ResidueOnlySelect(Select):
    def accept_residue(self, residue):
        return residue.id[0] == " "


def strip_non_residues(pdb_in, pdb_out):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("x", pdb_in)
    io = PDBIO()
    io.set_structure(struct)
    io.save(pdb_out, ResidueOnlySelect())


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


def extract_iptm_from_npz(npz_path):
    """Extract the iptm value from the .npz file."""
    with np.load(npz_path) as d:
        if "iptm" not in d:
            return None
        return float(d["iptm"].squeeze())


def extract_ptm_from_npz(npz_path):
    """Extract the ptm value from the .npz file."""
    with np.load(npz_path) as d:
        if "ptm" not in d:
            return None
        return float(d["ptm"].squeeze())


def extract_aggregate_from_npz(npz_path):
    """Extract the aggregate_score from the .npz file."""
    with np.load(npz_path) as d:
        if "aggregate_score" not in d:
            return None
        return float(d["aggregate_score"].squeeze())


def get_chain_seq_from_fasta(in_fasta, chain_id="A"):
    """
    Return the sequence for a given chain from a 2-entry FASTA.
    """
    for rec in SeqIO.parse(in_fasta, "fasta"):
        # Return the seq that matches up with chain_id in it's header
        if chain_id in rec.id:
            return str(rec.seq)
    # Fallback: if nothing matches, return the first record
    for rec in SeqIO.parse(in_fasta, "fasta"):
        return str(rec.seq)
    raise ValueError(f"No sequences found in {in_fasta}")


def split_fasta_entries(
    in_fasta,
    out_dir,
    molecule="protein",
    mode="binder_target",
    target_seq=None,
    target_chain="B",
    target_seq_2=None,
    target_chain_2="C",
):
    """
    Split design sequences into per-complex FASTA files.

    Modes:
      homodimer:
        Split input sequence equally into chains A and B.

      binder_target:
        Input sequence contains binder + target. target_seq identifies the
        target portion, which is written using target_chain.

      binder_dna:
        Input sequence is the binder sequence. target_seq and target_seq_2
        are written as DNA target chains (typically B and C).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if mode not in ("homodimer", "binder_target", "binder_dna"):
        raise ValueError(f"Unknown mode: {mode}")

    if mode == "binder_target" and not target_seq:
        raise ValueError("target_seq is required for mode='binder_target'")

    if mode == "binder_dna":
        if not target_seq:
            raise ValueError("target_seq is required for mode='binder_dna'")

        if not target_seq_2:
            raise ValueError(
                "target_seq_2 is required for mode='binder_dna'"
            )

        if target_chain == target_chain_2:
            raise ValueError(
                "target_chain and target_chain_2 must be different "
                f"(both are currently {target_chain!r})"
            )

        if target_chain == "A" or target_chain_2 == "A":
            raise ValueError(
                "In binder_dna mode, chain A is reserved for the binder. "
                "Use B and C for target_chain and target_chain_2."
            )

    count = 0

    for rec in SeqIO.parse(in_fasta, "fasta"):
        base = rec.id.split(",")[0].strip()
        full_seq = str(rec.seq).upper()

        if not full_seq:
            raise ValueError(f"Empty sequence for FASTA entry {base!r}")

        if mode == "homodimer":
            if len(full_seq) % 2 != 0:
                raise ValueError(
                    f"{base} sequence length is odd; cannot split evenly: "
                    f"{len(full_seq)}"
                )

            midpoint = len(full_seq) // 2

            seq_a = full_seq[:midpoint]
            seq_b = full_seq[midpoint:]

            records = [
                SeqRecord(
                    Seq(seq_a),
                    id=f"{molecule}|name={base}_A",
                    description="",
                ),
                SeqRecord(
                    Seq(seq_b),
                    id=f"{molecule}|name={base}_B",
                    description="",
                ),
            ]

        elif mode == "binder_target":
            target_seq = str(target_seq).upper()

            if len(full_seq) <= len(target_seq):
                raise ValueError(
                    f"Input sequence for {base} is not longer than target_seq "
                    f"(len(full)={len(full_seq)}, "
                    f"len(target)={len(target_seq)})"
                )

            if target_chain == "B":
                seq_a = full_seq[:-len(target_seq)]
                seq_b = target_seq

            elif target_chain == "A":
                seq_a = target_seq
                seq_b = full_seq[len(target_seq):]

            else:
                raise ValueError(
                    f"target_chain must be 'A' or 'B'; got {target_chain!r}"
                )

            if not seq_a or not seq_b:
                raise ValueError(
                    f"One chain is empty for {base}; verify target_seq "
                    "and target_chain."
                )

            records = [
                SeqRecord(
                    Seq(seq_a),
                    id=f"{molecule}|name={base}_A",
                    description="",
                ),
                SeqRecord(
                    Seq(seq_b),
                    id=f"{molecule}|name={base}_B",
                    description="",
                ),
            ]

        elif mode == "binder_dna":
            binder_seq = full_seq
            dna_seq_1 = str(target_seq).upper()
            dna_seq_2 = str(target_seq_2).upper()

            if not binder_seq:
                raise ValueError(f"Binder sequence is empty for {base!r}")

            if not set(dna_seq_1) <= {"A", "T", "G", "C"}:
                raise ValueError(
                    f"target_seq for {base} contains non-ATGC characters: "
                    f"{dna_seq_1}"
                )

            if not set(dna_seq_2) <= {"A", "T", "G", "C"}:
                raise ValueError(
                    f"target_seq_2 for {base} contains non-ATGC characters: "
                    f"{dna_seq_2}"
                )

            records = [
                SeqRecord(
                    Seq(binder_seq),
                    id=f"{molecule}|name={base}_A",
                    description="",
                ),
                SeqRecord(
                    Seq(dna_seq_1),
                    id=f"dna|name={base}_{target_chain}",
                    description="",
                ),
                SeqRecord(
                    Seq(dna_seq_2),
                    id=f"dna|name={base}_{target_chain_2}",
                    description="",
                ),
            ]

        else:
            raise RuntimeError("Unreachable mode branch")

        out_path = out_dir / f"{base}.fasta"

        with open(out_path, "w") as out_handle:
            SeqIO.write(records, out_handle, "fasta")

        count += 1

    return count


def pick_best_chai_output(output_dir, metrics_csv, interface=True):
    """Choose the best cif file from the chai output to be passed back to mpnn."""
    confidence_metrics = {}
    rows = []
    for cif_file in Path(output_dir).glob("*.cif"):
        npz_file = cif_file.with_suffix(".npz")
        name = cif_file.stem
        plddt = extract_plddt_from_cif(cif_file)
        iptm = extract_iptm_from_npz(npz_file)
        rows.append({
            "file": name,
            "iptm": iptm,
            "plddt": plddt,
        })
        if plddt is None or iptm is None:
            continue
        confidence_metrics[cif_file] = {"plddt": plddt, "iptm": iptm}
    save_metrics_csv(rows, metrics_csv)
        
    if not confidence_metrics:
        return None

    if interface:
        best_chai = max(confidence_metrics, key=lambda k: confidence_metrics[k]["iptm"])
    else:
        best_chai = max(confidence_metrics, key=lambda k: confidence_metrics[k]["plddt"])
    return best_chai.name


def rank_filtered_from_csv(metrics_file, output_file, interface=True, min_shape_complementarity=None, max_unsat_hbond=None, max_surface_hydrophobicity=None):
    """
    Read an unfiltered candidate metrics CSV, apply quality filters, rank
    survivors by iptm (interface=True) or plddt (interface=False), and write
    the results to output_file.
    """
    df = pd.read_csv(metrics_file)
    df.columns = df.columns.str.strip()

    required_columns = ["shape_complementarity","n_unsat_hbonds","surface_hydrophobicity","iptm","plddt"]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(
            f"Missing required columns in {metrics_file}: {missing_columns}. "
            f"Found columns: {df.columns.tolist()}"
        )

    numeric_columns = ["shape_complementarity","n_unsat_hbonds","surface_hydrophobicity","iptm","plddt"]
    for col in numeric_columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    print(
        f"Filtering {len(df)} candidates from {metrics_file}",
        file=sys.stderr,
    )

    if min_shape_complementarity is not None:
        before = len(df)
        df = df[df["shape_complementarity"] >= min_shape_complementarity]
        print(
            f"SC >= {min_shape_complementarity}: {len(df)}/{before} remain",
            file=sys.stderr,
        )
    '''
    if max_unsat_hbond is not None:
        before = len(df)
        df = df[df["n_unsat_hbonds"] <= max_unsat_hbond]
        print(
            f"n_unsat_hbonds <= {max_unsat_hbond}: {len(df)}/{before} remain",
            file=sys.stderr,
        )
    '''
    if max_surface_hydrophobicity is not None:
        before = len(df)
        df = df[
            df["surface_hydrophobicity"] <= max_surface_hydrophobicity
        ]
        print(
            f"surface_hydrophobicity <= {max_surface_hydrophobicity}: "
            f"{len(df)}/{before} remain",
            file=sys.stderr,
        )
    sort_col = "iptm" if interface else "plddt"

    df = df.sort_values(by=sort_col,ascending=False,na_position="last").reset_index(drop=True)

    df["rank"] = df.index + 1

    # This writes filtered/ranked candidates to a NEW file.
    df.to_csv(output_file, index=False)

    print(
        f"Wrote {len(df)} passing candidates to {output_file}",
        file=sys.stderr,
    )

    return df


def update_convergence(best_iptm, best_plddt, cif_path, npz_path, interface=True,
                       iptm_margin=0.01, plddt_margin=0.005):
    """Determine whether prediction is improving or not."""
    new_plddt = extract_plddt_from_cif(cif_path)
    new_iptm = extract_iptm_from_npz(npz_path)

    if new_plddt is None or new_iptm is None:
        return False, best_iptm, best_plddt, new_iptm, new_plddt

    if interface:
        improved = (new_iptm - best_iptm) >= iptm_margin
    else:
        improved = (new_plddt - best_plddt) >= plddt_margin

    return improved, max(best_iptm, new_iptm), max(best_plddt, new_plddt), new_iptm, new_plddt


def save_metrics_csv(rows, csv_path):
    """Save plddt and iptm for each predicted structure."""
    fieldnames = ["file", "iptm", "plddt"]
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError(f"Length mismatch: {len(seq1)} vs {len(seq2)}")
    return sum(a != b for a, b in zip(seq1, seq2))


def extract_seq_from_pdb(pdb_path, chain_mode="pdb-atom"):
    """
    Return concatenated sequence from a PDB file.
    """
    records = list(SeqIO.parse(pdb_path, chain_mode))
    if not records:
        raise ValueError(f"No sequences found in {pdb_path} using mode {chain_mode}")
    # If only care about one chain, take the first
    return str(records[0].seq)


def extract_chain_seq_from_cif(cif_path, chain_id, mode="cif-atom"):
    """
    Return the sequence of a specific chain from a CIF file.
    """
    records = list(SeqIO.parse(cif_path, mode))
    if not records:
        raise ValueError(f"No sequences found in {cif_path} using mode {mode}")

    chain_id = chain_id.strip().upper()

    for rec in records:
        # rec.id might be like '????:A' or similar
        chain_part = rec.id.split(":")[-1].strip().upper()
        if chain_part == chain_id:
            return str(rec.seq)

    raise ValueError(f"Chain {chain_id} not found in {cif_path}")


def extract_chain_seq_from_pdb(pdb_path, chain_id, chain_mode="pdb-atom"):
    """
    Return the sequence of a specific chain from a PDB file.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)
    model = structure[0]

    chain_id = chain_id.strip().upper()

    if chain_id not in model:
        available = ", ".join(
            repr(chain.id)
            for chain in model
        )
        raise ValueError(
            f"Chain {chain_id!r} not found in {pdb_path}. "
            f"Available chains: {available}"
        )

    chain = model[chain_id]

    amino_acid_map = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
        "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
        "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
        "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
        "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }

    nucleotide_map = {
        "DA": "A",
        "DC": "C",
        "DG": "G",
        "DT": "T",
        "A": "A",
        "C": "C",
        "G": "G",
        "U": "U",
    }

    sequence = []

    for residue in chain:
        residue_name = residue.get_resname().strip().upper()

        if is_aa(residue, standard=True):
            sequence.append(amino_acid_map[residue_name])

        elif is_nucleic(residue, standard=False):
            sequence.append(nucleotide_map.get(residue_name, "N"))

    if not sequence:
        observed_residues = sorted({
            residue.get_resname().strip()
            for residue in chain
        })

        raise ValueError(
            f"Chain {chain_id!r} exists in {pdb_path}, but no recognizable "
            f"protein or nucleic-acid residues were found. "
            f"Observed residue names: {observed_residues}"
        )

    return "".join(sequence)


def convert_cif_to_pdb(cif_path, pdb_path):
    """Convert mmCIF to PDB using Biopython."""
    parser = MMCIFParser(QUIET=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        structure = parser.get_structure(Path(cif_path).stem, cif_path)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path)


def pdb_to_fasta(pdb_path, out_fasta, chain_mode="pdb-atom", molecule="protein", name=None):
    """
    Create a FASTA file from a PDB.

    - Writes a single FASTA entry:
        >{molecule}|name={name}
        {sequence}
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)
    model = structure[0]

    if name is None:
        name = Path(pdb_path).stem

    records = []

    for chain in model:
        raw_chain_id = chain.id
        chain_id = raw_chain_id.strip() or "blank"

        try:
            sequence = extract_chain_seq_from_pdb(
                pdb_path=pdb_path,
                chain_id=raw_chain_id,
                chain_mode=chain_mode,
            ).upper()
        except ValueError as error:
            print(
                f"Skipping chain {raw_chain_id!r}: {error}",
                file=sys.stderr,
            )
            continue

        if not sequence:
            continue

        molecule = (
            "dna"
            if set(sequence) <= {"A", "T", "G", "C"}
            else "protein"
        )

        header = f"{molecule}|name={name}_{chain_id}"

        records.append(
            SeqRecord(
                Seq(sequence),
                id=header,
                description="",
            )
        )

    if not records:
        raise ValueError(
            f"No protein or nucleic-acid sequences could be extracted "
            f"from {pdb_path}"
        )

    output_path = Path(out_fasta)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as handle:
        SeqIO.write(records, handle, "fasta")


def interface_eval(input_path, relaxed_pdb, binder_chain, target_chain, max_iter=200):
    """
    Run BindCraft-style FastRelax + interface scoring and return core metrics.
    """
    # Initialize PyRosetta once
    init_pyrosetta(silent=True)

    # If input is CIF, convert to a temporary PDB first
    if str(input_path).lower().endswith(".cif"):
        tmp_pdb = str(Path(relaxed_pdb).with_suffix(".input.pdb"))
        convert_cif_to_pdb(input_path, tmp_pdb)
        pdb_in = tmp_pdb
    else:
        pdb_in = input_path

    # Run FastRelax
    pr_relax(pdb_in, relaxed_pdb, max_iter=max_iter)

    # Interface scoring
    iface_scores = score_interface(
        relaxed_pdb, binder_chain=binder_chain, target_chain=target_chain)
    iface_AA = None
    iface_residues = None

    # Additional metrics
    binder_seq = get_binder_sequence(relaxed_pdb, binder_chain=binder_chain)
    binder_length = len(binder_seq) if binder_seq else 0
    shape_complementarity = iface_scores.get("shape_complementarity", "nan")
    unsat_hbond_pct = iface_scores.get("unsat_hbond_pct", "nan")
    n_interface_residues = iface_scores.get("n_interface_residues", "nan")
    surface_hydrophobicity = iface_scores.get("surface_hydrophobicity", "nan")
    interface_sasa_pct = iface_scores.get("interface_sasa_pct", "nan")

    # Secondary structure fractions
    ss = calc_ss_percentage(relaxed_pdb, binder_chain=binder_chain)
    binder_helix_pct    = ss.get('binder_helix_pct', 'nan')
    binder_sheet_pct    = ss.get('binder_sheet_pct', 'nan')
    binder_loop_pct     = ss.get('binder_loop_pct', 'nan')
    interface_helix_pct = ss.get('interface_helix_pct', 'nan')
    interface_sheet_pct = ss.get('interface_sheet_pct', 'nan')
    interface_loop_pct  = ss.get('interface_loop_pct', 'nan')

    return {
        "dG": iface_scores.get("dG", "nan"),
        "dSASA": iface_scores.get("dSASA", "nan"),
        "dG_dSASA_ratio": iface_scores.get("dG_dSASA_ratio", "nan"),
        "interface_hydrophobicity": iface_scores.get("interface_hydrophobicity", "nan"),
        "n_unsat_hbonds": iface_scores.get("n_unsat_hbonds", "nan"),
        "binder_seq": binder_seq,
        "binder_length": binder_length,
        "shape_complementarity": shape_complementarity,
        "unsat_hbond_pct": unsat_hbond_pct,
        "n_interface_residues": n_interface_residues,
        "surface_hydrophobicity": surface_hydrophobicity,
        "interface_sasa_pct": interface_sasa_pct,
        "binder_helix_pct": binder_helix_pct,
        "binder_sheet_pct": binder_sheet_pct,
        "binder_loop_pct": binder_loop_pct,
        "interface_helix_pct": interface_helix_pct,
        "interface_sheet_pct": interface_sheet_pct,
        "interface_loop_pct": interface_loop_pct
    }


def load_structure(path):
    p = str(path).lower()
    if p.endswith(".cif") or p.endswith(".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
    return parser.get_structure("s", path)


def complex_rmsd(path1, path2):
    """
    Compute complex RMSD between two structures (PDB or CIF).
    """
    ref_struct = load_structure(path1)
    mobile_struct = load_structure(path2)

    ref_atoms = [atom for atom in ref_struct.get_atoms() if atom.get_name() == "CA"]
    mobile_atoms = [atom for atom in mobile_struct.get_atoms() if atom.get_name() == "CA"]

    if len(ref_atoms) != len(mobile_atoms):
        raise ValueError(
            f"CA atom count mismatch: {len(ref_atoms)} vs {len(mobile_atoms)} "
            f"between {path1} and {path2}"
        )

    super_imposer = Superimposer()
    super_imposer.set_atoms(ref_atoms, mobile_atoms)
    return super_imposer.rms


def get_interface_residues(pdb_path, binder_chain, target_chain, cutoff):
    """
    Return PDB-numbered binder residues with a certain distance of target residues.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_path)
    model = structure[0]
    if binder_chain not in model:
        raise ValueError(
            f"Binder chain {binder_chain!r} not found in {pdb_path}"
        )
    if target_chain not in model:
        raise ValueError(
            f"Target chain {target_chain!r} not found in {pdb_path}"
        )
    binder_residues = [
        residue
        for residue in model[binder_chain]
        if is_aa(residue, standard=True)
    ]
    target_atoms = [
        atom
        for residue in model[target_chain]
        if is_aa(residue, standard=True)
        for atom in residue
        if atom.element != "H"
    ]
    if not target_atoms:
        raise ValueError(
            f"No protein heavy atoms found in target chain {target_chain!r}"
        )
    neighbor_search = NeighborSearch(target_atoms)
    interface_positions = set()

    for sequence_position, residue in enumerate(binder_residues, start=1):
        for atom in residue:
            if atom.element == "H":
                continue

            contacts = neighbor_search.search(
                atom.coord,
                cutoff,
                level="R",
            )

            if contacts:
                # Return 1-based position within the binder chain
                interface_positions.add(sequence_position)
                break
    return ",".join(str(position) for position in sorted(interface_positions))


def fixed_residues(fixed_residue_file, binder_chain, target_sequence, target_chain):
    """
    Build a list of residue IDs to fix in MPNN.
    """

    if not binder_chain or not binder_chain.isalpha():
        raise ValueError(
            f"Binder chain must contain letters, got: {binder_chain!r}"
        )

    if not target_chain or not target_chain.isalpha():
        raise ValueError(
            f"Target chain must contain letters, got: {target_chain!r}"
        )

    binder_positions = set()

    with open(fixed_residue_file) as handle:
        for line in handle:
            line = line.strip()

            # Ignore blank lines and  section headers.
            if not line or line.startswith("#"):
                continue

            for item in line.split(","):
                item = item.strip()

                # Prevents empty MPNN IDs from malformed commas.
                if not item:
                    continue

                try:
                    position = int(item)
                except ValueError as error:
                    raise ValueError(
                        f"Invalid residue position {item!r} in "
                        f"{fixed_residue_file}"
                    ) from error

                if position < 1:
                    raise ValueError(
                        f"Residue positions must be >= 1; found {position} "
                        f"in {fixed_residue_file}"
                    )

                binder_positions.add(position)

    target_sequence = "".join(target_sequence.split()).upper()

    if not target_sequence:
        raise ValueError("Target sequence is empty.")

    fixed_ids = [
        f"{binder_chain}{position}"
        for position in sorted(binder_positions)
    ]

    fixed_ids.extend(
        f"{target_chain}{position}"
        for position in range(1, len(target_sequence) + 1)
    )

    return "[" + ",".join(f'"{residue_id}"' for residue_id in fixed_ids) + "]"


def main():
    cmd = sys.argv[1]
    if cmd == "split-fasta":
        in_fasta = sys.argv[2]
        out_dir = sys.argv[3]
        mode = sys.argv[4]
        target_seq = sys.argv[5] if len(sys.argv) > 5 else None
        target_chain = sys.argv[6].strip().upper()
        target_seq_2 = sys.argv[7] if len(sys.argv) > 7 else None
        target_chain_2 = sys.argv[8].strip().upper() if len(sys.argv) > 8 else None
        molecule = sys.argv[9] if len(sys.argv) > 9 else "protein"
        n = split_fasta_entries(
                in_fasta, 
                out_dir, 
                molecule=molecule, 
                mode=mode, 
                target_seq=target_seq,
                target_chain=target_chain, 
                target_seq_2=target_seq_2, 
                target_chain_2=target_chain_2)
        print(n)
    elif cmd == "pick-best":
        interface = sys.argv[4].lower() == "true"
        best = pick_best_chai_output(sys.argv[2], sys.argv[3], interface=interface)
        print(best)
    elif cmd == "rank-filtered-csv":
        if len(sys.argv) < 4:
            raise ValueError(
                "Usage: chai_mpnn_utils.py rank-filtered-csv "
                "<input_metrics_csv> <output_filtered_csv> "
                "[interface=true|false] [min_sc] [max_unsat] [max_surf_h]"
            )   

        metrics_csv = sys.argv[2]
        filtered_csv = sys.argv[3]
        interface = sys.argv[4].lower() == "true" if len(sys.argv) > 4 else True

        def maybe_float(value):
            if value is None or value.lower() == "none":
                return None
            return float(value)

        min_sc = maybe_float(sys.argv[5]) if len(sys.argv) > 5 else None
        max_unsat = maybe_float(sys.argv[6]) if len(sys.argv) > 6 else None
        max_surf_h = maybe_float(sys.argv[7]) if len(sys.argv) > 7 else None

        rank_filtered_from_csv(metrics_file=metrics_csv,output_file=filtered_csv,interface=interface,min_shape_complementarity=min_sc,max_unsat_hbond=max_unsat,max_surface_hydrophobicity=max_surf_h)
    elif cmd == "update-convergence":
        best_iptm = float(sys.argv[2])
        best_plddt = float(sys.argv[3])
        cif_path = sys.argv[4]
        npz_path = sys.argv[5]
        interface = sys.argv[6].lower() == "true"
        iptm_margin = float(sys.argv[7])
        plddt_margin = float(sys.argv[8])

        improved, new_best_iptm, new_best_plddt, cur_iptm, cur_plddt = update_convergence(
            best_iptm, best_plddt, cif_path, npz_path,
            interface=interface,
            iptm_margin=iptm_margin,
            plddt_margin=plddt_margin,
        )
        print(int(improved), new_best_iptm, new_best_plddt)
    elif cmd == "metrics":
        cif_path = sys.argv[2]
        npz_path = sys.argv[3]
        iptm = extract_iptm_from_npz(npz_path)
        plddt = extract_plddt_from_cif(cif_path)
        ptm = extract_ptm_from_npz(npz_path)
        agg = extract_aggregate_from_npz(npz_path)
        print(iptm, plddt, ptm, agg)
    elif cmd == "hamming-dist":
        path1 = sys.argv[2]
        path2 = sys.argv[3]
        # Default to binder chain for binder_target mode
        chain_id = sys.argv[4] if len(sys.argv) > 4 else "A"

        if is_pdb(path1) and is_cif(path2):
            seq1 = extract_chain_seq_from_pdb(path1, chain_id, chain_mode="pdb-atom")
            seq2 = extract_chain_seq_from_cif(path2, chain_id, mode="cif-atom")
        elif is_cif(path1) and is_cif(path2):
            seq1 = extract_chain_seq_from_cif(path1, chain_id, mode="cif-atom")
            seq2 = extract_chain_seq_from_cif(path2, chain_id, mode="cif-atom")
        else:
            raise ValueError(
                f"hamming-dist expects either PDB+CIF or CIF+CIF; got {path1} and {path2}"
            )

        dist = hamming_distance(seq1, seq2)
        print(dist)
    elif cmd == "interface-eval":
        input_pdb = sys.argv[2]
        relaxed_pdb = sys.argv[3]
        binder_chain = sys.argv[4]
        target_chain = sys.argv[5]
        metrics = interface_eval(input_pdb, relaxed_pdb, binder_chain, target_chain, max_iter=200)
        # Print as space-separated values in a fixed order
        order = [
            "dG", "dSASA", "dG_dSASA_ratio", "interface_hydrophobicity", "n_unsat_hbonds",
            "binder_seq", "binder_length", "shape_complementarity", "unsat_hbond_pct",
            "n_interface_residues", "surface_hydrophobicity", "interface_sasa_pct",
            "binder_helix_pct", "binder_sheet_pct", "binder_loop_pct",
            "interface_helix_pct", "interface_sheet_pct", "interface_loop_pct"
        ]
        print(*(metrics[k] for k in order))
    elif cmd == "complex-rmsd":
        path1 = sys.argv[2]
        path2 = sys.argv[3]
        rmsd_val = complex_rmsd(path1, path2)
        print(rmsd_val)
    elif cmd == "get-chain-seq":
        pdb_path = sys.argv[2]
        chain_id = sys.argv[3] if len(sys.argv) > 3 else "B"
        mode = sys.argv[4] if len(sys.argv) > 4 else "pdb-atom"
        seq = extract_chain_seq_from_pdb(pdb_path, chain_id, chain_mode=mode)
        print(seq)
    elif cmd == "pdb-to-fasta":
        pdb_path = sys.argv[2]
        out_fasta = sys.argv[3] if len(sys.argv) > 3 else pdb_path.replace(".pdb", ".fasta")
        chain_mode = sys.argv[4] if len(sys.argv) > 4 else "pdb-atom"
        molecule = sys.argv[5] if len(sys.argv) > 5 else "protein"
        name = sys.argv[6] if len(sys.argv) > 6 else None
        pdb_to_fasta(pdb_path, out_fasta, chain_mode=chain_mode, molecule=molecule, name=name)
    elif cmd == "strip-non-residues":
        pdb_in = sys.argv[2]
        pdb_out = sys.argv[3]
        strip_non_residues(pdb_in, pdb_out)
    elif cmd == "interface-residues":
        pdb_path = sys.argv[2]
        binder_chain = sys.argv[3]
        target_chain = sys.argv[4]
        cutoff = float(sys.argv[5]) if len(sys.argv) > 5 else 7.5
        residues = get_interface_residues(
            pdb_path=pdb_path,
            binder_chain=binder_chain,
            target_chain=target_chain,
            cutoff=cutoff,
        )
        print(residues)
    elif cmd == "fixed-residues":
        fixed_residue_file = sys.argv[2]
        binder_chain = sys.argv[3]
        target_sequence = sys.argv[4]
        target_chain = sys.argv[5]
        print(fixed_residues(
                fixed_residue_file=fixed_residue_file,
                binder_chain=binder_chain,
                target_sequence=target_sequence,
                target_chain=target_chain,
            )
        )
    elif cmd == "cif-to-pdb":
        if len(sys.argv) != 4:
            raise ValueError(
                "Usage: chai_mpnn_utils.py cif-to-pdb <input.cif> <output.pdb>"
            )
        cif_path = sys.argv[2]
        pdb_path = sys.argv[3]
        convert_cif_to_pdb(cif_path, pdb_path)
        print(pdb_path)


if __name__ == "__main__":
    main()
