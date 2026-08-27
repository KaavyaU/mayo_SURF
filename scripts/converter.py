import os
import sys
import argparse
import warnings
from Bio.PDB import MMCIFParser, PDBIO


def convert_cif_to_pdb(cif_path: str, pdb_path: str):
    parser = MMCIFParser(QUIET=True)

    structure_id = os.path.splitext(os.path.basename(cif_path))[0]

    print(f"Parsing: {cif_path}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = parser.get_structure(structure_id, cif_path)

    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path)

    print(f"Written: {pdb_path}")


def main(args: argparse.Namespace):
    if os.path.isdir(args.input):
        # Directory mode: convert all .cif files found in input dir
        cif_files = [
            f for f in os.listdir(args.input)
            if f.lower().endswith(".cif")
        ]
        if not cif_files:
            print(f"[ERROR] No .cif files found in directory: {args.input}")
            sys.exit(1)

        # Determine output directory
        if args.output:
            out_dir = args.output
        else:
            out_dir = os.path.join(args.input, "pdb")

        os.makedirs(out_dir, exist_ok=True)
        print(f"Output directory: {out_dir}")

        for fname in sorted(cif_files):
            cif_path = os.path.join(args.input, fname)
            stem = os.path.splitext(fname)[0]
            if stem.endswith("_0"):
                stem = stem[:-2]
            pdb_name = stem + ".pdb"
            pdb_path = os.path.join(out_dir, pdb_name)
            convert_cif_to_pdb(cif_path, pdb_path)

    elif os.path.isfile(args.input):
        # Single file mode
        if args.output:
            pdb_path = args.output
        else:
            pdb_path = os.path.splitext(args.input)[0] + ".pdb"

        os.makedirs(os.path.dirname(os.path.abspath(pdb_path)), exist_ok=True)
        convert_cif_to_pdb(args.input, pdb_path)

    else:
        print(f"[ERROR] Input not found: {args.input}")
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert .cif (mmCIF) file(s) to .pdb format.",
        epilog="""
Examples:
  # Single file — output alongside input
  python converter.py structure.cif

  # Single file — explicit output path
  python converter.py structure.cif --output ./pdbs/my_structure.pdb

  # Directory — converted files written to ./cifs/pdb/ (auto-created)
  python converter.py ./cifs/

  # Directory — converted files written to a specific output directory
  python converter.py ./cifs/ --output ./pdbs/
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "input",
        help="Path to a .cif file, or a directory containing .cif files",
    )
    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Output path. For a file input: path to the .pdb file "
            "(default: same location as input with .pdb extension). "
            "For a directory input: path to the output directory "
            "(default: <input_dir>/pdb/)"
        ),
    )
    args = parser.parse_args()
    main(args)
