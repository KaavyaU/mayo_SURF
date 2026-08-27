#!/usr/bin/env python3

import argparse
import random
import tarfile
import time
from collections import Counter
from pathlib import Path

import numpy as np
import requests
from Bio import SeqIO


def run_mmseqs2(
    out_path: str,
    sequence: str,
    use_filter: bool = True,
    host_url: str = "https://api.colabfold.com",
) -> str:
    """
    Submit one protein sequence to the ColabFold MMseqs2 API.

    Returns:
        Path to the extracted UniRef A3M alignment.
    """
    out_dir = Path(out_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    tar_file = out_dir / "out.tar.gz"
    uniref_a3m = out_dir / "uniref.a3m"

    # Reuse an existing valid alignment.
    if uniref_a3m.is_file() and uniref_a3m.stat().st_size > 0:
        return str(uniref_a3m)

    def submit(sequence_to_search: str, mode: str) -> dict:
        response = requests.post(
            f"{host_url}/ticket/msa",
            data={
                "q": f">query\n{sequence_to_search}\n",
                "mode": mode,
            },
            timeout=60,
        )
        response.raise_for_status()
        return response.json()

    def get_status(job_id: str) -> dict:
        response = requests.get(
            f"{host_url}/ticket/{job_id}",
            timeout=60,
        )
        response.raise_for_status()
        return response.json()

    def download(job_id: str) -> None:
        response = requests.get(
            f"{host_url}/result/download/{job_id}",
            timeout=300,
        )
        response.raise_for_status()

        with open(tar_file, "wb") as handle:
            handle.write(response.content)

    # Download the archive only if one is not already cached.
    if not tar_file.is_file() or tar_file.stat().st_size == 0:
        mode = "env" if use_filter else "env-nofilter"

        while True:
            ticket = submit(sequence, mode)
            ticket_status = ticket.get("status")

            while ticket_status in {"UNKNOWN", "RATELIMIT"}:
                wait_seconds = 10 + random.randint(0, 5)
                print(
                    f"MMseqs2 server status: {ticket_status}. "
                    f"Waiting {wait_seconds} s..."
                )
                time.sleep(wait_seconds)

                ticket = submit(sequence, mode)
                ticket_status = ticket.get("status")

            if ticket_status == "MAINTENANCE":
                raise RuntimeError(
                    "The ColabFold MMseqs2 API is under maintenance. "
                    "Try again later."
                )

            if ticket_status == "ERROR":
                raise RuntimeError(
                    "The ColabFold MMseqs2 API rejected the sequence search."
                )

            job_id = ticket.get("id")

            if not job_id:
                raise RuntimeError(
                    f"MMseqs2 API returned no job ID: {ticket}"
                )

            print(f"MMseqs2 job submitted: {job_id}")

            while ticket_status in {"UNKNOWN", "PENDING", "RUNNING"}:
                wait_seconds = 5 + random.randint(0, 5)
                print(
                    f"MMseqs2 job status: {ticket_status}. "
                    f"Waiting {wait_seconds} s..."
                )
                time.sleep(wait_seconds)

                ticket = get_status(job_id)
                ticket_status = ticket.get("status")

            if ticket_status == "COMPLETE":
                print("MMseqs2 search complete; downloading results...")
                download(job_id)
                break

            if ticket_status == "ERROR":
                raise RuntimeError("MMseqs2 search failed.")

            raise RuntimeError(
                f"Unexpected final MMseqs2 job status: {ticket}"
            )

    # Extract only uniref.a3m from the API archive.
    if not uniref_a3m.is_file() or uniref_a3m.stat().st_size == 0:
        with tarfile.open(tar_file) as archive:
            matches = [
                member
                for member in archive.getmembers()
                if member.name.endswith("uniref.a3m")
            ]

            if not matches:
                archive_contents = "\n".join(
                    member.name
                    for member in archive.getmembers()
                )

                raise FileNotFoundError(
                    "Could not find uniref.a3m in the MMseqs2 archive.\n"
                    f"Archive contents:\n{archive_contents}"
                )

            member = matches[0]

            # Avoid retaining nested directories from the archive.
            member.name = Path(member.name).name

            archive.extract(
                member,
                path=out_dir,
                filter="data",
            )

    if not uniref_a3m.is_file() or uniref_a3m.stat().st_size == 0:
        raise FileNotFoundError(
            f"Expected extracted alignment was not found: {uniref_a3m}"
        )

    return str(uniref_a3m)


def select_fasta_record(records, chain: str | None):
    """
    Select exactly one sequence record.

    `--chain A` matches records named:
      - A
      - target_A
      - target:A
    """
    if not records:
        raise ValueError("The input FASTA contains no sequence records.")

    if chain is None:
        if len(records) == 1:
            return records[0]

        available = "\n".join(f"  - {record.id}" for record in records)

        raise ValueError(
            "The input FASTA contains multiple records. "
            "Use --chain to select one.\n"
            f"Available FASTA IDs:\n{available}"
        )

    chain = chain.strip()

    matches = [
        record
        for record in records
        if (
            record.id == chain
            or record.name == chain
            or record.description.split()[0] == chain
            or record.id.endswith(f"_{chain}")
            or record.id.endswith(f":{chain}")
        )
    ]

    if not matches:
        available = "\n".join(f"  - {record.id}" for record in records)

        raise ValueError(
            f"Could not find chain {chain!r}.\n"
            f"Available FASTA IDs:\n{available}"
        )

    if len(matches) > 1:
        matched_ids = ", ".join(record.id for record in matches)

        raise ValueError(
            f"Chain {chain!r} matches multiple FASTA records: {matched_ids}.\n"
            "Use the exact FASTA record identifier with --chain."
        )

    return matches[0]


def parse_a3m(a3m_path: str) -> tuple[str, list[str]]:
    """
    Parse an A3M alignment and remove lowercase insertion columns.

    fasta-blast accepts MMseqs2/ColabFold comment lines beginning with #.
    """
    records = list(SeqIO.parse(a3m_path, "fasta-blast"))

    if not records:
        raise ValueError(f"No sequences found in A3M file: {a3m_path}")

    query_seq = str(records[0].seq)

    keep_cols = [
        index
        for index, residue in enumerate(query_seq)
        if not residue.islower()
    ]

    aligned = []

    for record in records:
        sequence = str(record.seq)

        stripped = "".join(
            sequence[index]
            for index in keep_cols
            if index < len(sequence)
        )

        aligned.append(stripped.upper())

    return aligned[0], aligned


def compute_conservation(query: str, aligned: list[str]) -> list[dict]:
    """
    Compute amino-acid conservation at each query residue.

    Score = frequency of the most frequent non-gap amino acid in that
    alignment column.
    """
    results = []

    for position in range(len(query)):
        column = [
            sequence[position]
            for sequence in aligned
            if position < len(sequence)
        ]

        residues = [
            amino_acid
            for amino_acid in column
            if amino_acid != "-"
        ]

        if not residues:
            results.append(
                {
                    "position": position + 1,
                    "query_aa": query[position],
                    "top_aa": "-",
                    "top_freq": 0.0,
                    "counts": {},
                    "n_seqs": 0,
                }
            )
            continue

        counts = Counter(residues)
        top_aa, top_count = counts.most_common(1)[0]

        results.append(
            {
                "position": position + 1,
                "query_aa": query[position],
                "top_aa": top_aa,
                "top_freq": top_count / len(residues),
                "counts": dict(counts),
                "n_seqs": len(residues),
            }
        )

    return results


def select_fixed_positions(
    results: list[dict],
    percentile: int,
) -> list[int]:
    """
    Select the top `percentile` percent positions by conservation score.
    """
    if not 0 <= percentile <= 100:
        raise ValueError("--percentile must be an integer from 0 through 100.")

    ranked = sorted(
        results,
        key=lambda result: result["top_freq"],
        reverse=True,
    )

    number_to_keep = max(
        1,
        int(np.ceil(len(ranked) * percentile / 100)),
    )

    return [
        result["position"]
        for result in ranked[:number_to_keep]
    ]


def write_conserved_positions(
    conserved_positions: list[int],
    percentile: int,
    additional_positions: list[int],
    interface_positions: list[int],
    out_path: str,
) -> None:
    """
    Write conserved, manually-added, and interface residues in separate blocks.
    """
    output_file = Path(out_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w") as handle:
        handle.write(
            f"# Top {percentile}% most conserved positions "
            f"({len(conserved_positions)} residues)\n"
        )
        handle.write(
            ", ".join(map(str, sorted(set(conserved_positions)))) + "\n\n"
        )

        if additional_positions:
            handle.write(
                f"# Additional manually specified residues "
                f"({len(set(additional_positions))} residues)\n"
            )
            handle.write(
                ", ".join(map(str, sorted(set(additional_positions)))) + "\n\n"
            )

        if interface_positions:
            handle.write(
                f"# Interface residues "
                f"({len(set(interface_positions))} residues)\n"
            )
            handle.write(
                ", ".join(map(str, sorted(set(interface_positions)))) + "\n"
            )

def parse_residue_list(value: str) -> list[int]:
    """
    Parse comma-separated 1-indexed positions, e.g. '1,139'.
    Empty input returns an empty list.
    """
    if not value or not value.strip():
        return []

    try:
        residues = [
            int(item.strip())
            for item in value.split(",")
            if item.strip()
        ]
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "Expected comma-separated residue positions, e.g. '1,139'"
        ) from error

    if any(position < 1 for position in residues):
        raise argparse.ArgumentTypeError(
            "Residue positions must be positive 1-indexed integers."
        )

    return sorted(set(residues))


def main(args: argparse.Namespace) -> None:
    records = list(SeqIO.parse(args.fasta_file, "fasta"))
    record = select_fasta_record(records, args.chain)

    sequence = str(record.seq).upper()

    if not sequence.isalpha():
        raise ValueError(
            f"Selected record {record.id!r} contains non-alphabetical "
            "characters. Provide a protein-only FASTA sequence."
        )

    safe_tag = record.id.replace("/", "_").replace("\\", "_")
    alignment_dir = Path(args.alignment_dir) / safe_tag

    print(f"Selected record: {record.id}")
    print(f"Sequence length: {len(sequence)} residues")

    print("Generating or reusing MMseqs2 alignment...")
    a3m_path = run_mmseqs2(
        out_path=str(alignment_dir),
        sequence=sequence,
        use_filter=not args.no_filter,
        host_url=args.host_url,
    )

    print(f"Using A3M: {a3m_path}")
    print("Parsing alignment...")

    query, aligned = parse_a3m(a3m_path)

    print(f"Query length: {len(query)} residues")
    print(f"Sequences in MSA: {len(aligned)}")

    print("Computing conservation...")
    results = compute_conservation(query, aligned)

    if args.percentile == 0:
        conserved_positions = []
    else:
        conserved_positions = select_fixed_positions(
            results,
            percentile=args.percentile,
        )

    additional_positions = args.additional_residues
    interface_positions = args.interface_residues

    if args.output is None:
        output_file = alignment_dir / (
            f"conserved_positions_top_{args.percentile}.txt"
        )
    else:
        output_file = Path(args.output)

    write_conserved_positions(
        conserved_positions=conserved_positions,
        percentile=args.percentile,
        interface_positions=interface_positions,
        additional_positions=additional_positions,
        out_path=str(output_file),
    )
    print(
        f"Selected {len(conserved_positions)} conservation-derived positions "
        f"(top {args.percentile}% by conservation)."
    )   

    if additional_positions:
        print(f"Added {len(additional_positions)} manually specified positions.")

    if interface_positions:
        print(f"Added {len(interface_positions)} interface positions.")

    print(f"Conserved-position file: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Search one protein chain with the ColabFold/MMseqs2 API and "
            "write its most conserved residue positions."
        )
    )

    parser.add_argument(
        "fasta_file",
        help="Protein FASTA file containing one or more records.",
    )

    parser.add_argument(
        "--chain",
        default=None,
        help=(
            "Chain or FASTA record ID to analyze. Required if the input "
            "FASTA contains multiple records."
        ),
    )

    parser.add_argument(
        "--percentile",
        type=int,
        default=30,
        help=(
            "Retain the top N percent of positions ranked by conservation. "
            "Default: 30."
        ),
    )

    parser.add_argument(
        "--alignment-dir",
        default="alignments",
        help="Directory to cache MMseqs2 archive and A3M output.",
    )

    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Output text-file path. Defaults to a conserved-position file "
            "inside the chain-specific alignment directory."
        ),
    )

    parser.add_argument(
        "--interface-residues",
        type=parse_residue_list,
        default=[],
        help=(
            "Structure-derived interface positions, comma-separated; "
            "e.g. '7,11,14'. Written separately from conserved positions."
        ),
    )

    parser.add_argument(
        "--additional-residues",
        type=parse_residue_list,
        default=[],
    )

    parser.add_argument(
        "--no-filter",
        action="store_true",
        help="Use the MMseqs2 env-nofilter search mode.",
    )

    parser.add_argument(
        "--host-url",
        default="https://api.colabfold.com",
        help="ColabFold MMseqs2 API endpoint.",
    )

    main(parser.parse_args())
