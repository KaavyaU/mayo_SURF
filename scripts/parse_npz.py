from pathlib import Path
import csv
import numpy as np

def extract_iptm_from_npz(npz_path):
    with np.load(npz_path) as d:
        if "iptm" not in d:
            return None
        return float(d["iptm"].squeeze())

def short_name(path):
    # 1_1_0.npz -> 1_1
    parts = Path(path).stem.split("_")
    return "_".join(parts[:2]) if len(parts) >= 2 else Path(path).stem

def main(in_dir, out_csv):
    in_dir = Path(in_dir)
    out_csv = Path(out_csv)

    rows = []
    for npz_file in sorted(in_dir.glob("*.npz")):
        rows.append({
            "file": short_name(npz_file),
            "iptm": extract_iptm_from_npz(npz_file),
            "npz_file": npz_file.name,
        })

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["file", "iptm", "npz_file"])
        writer.writeheader()
        writer.writerows(rows)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True)
    parser.add_argument("-o", "--output_csv", required=True)
    args = parser.parse_args()
    main(args.input_dir, args.output_csv)
