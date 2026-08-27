from pathlib import Path
from itertools import combinations_with_replacement
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse

def make_pairwise_fastas(in_fasta, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    records = list(SeqIO.parse(in_fasta, "fasta"))
    count = 0

    for rec1, rec2 in combinations_with_replacement(records, 2):
        base1 = rec1.id.split(",")[0].strip()
        base2 = rec2.id.split(",")[0].strip()

        out_name = f"{base1}_{base2}.fasta"
        out_path = out_dir / out_name

        pair = [
            SeqRecord(rec1.seq, id=f"protein|name={rec1.id}", description=""),
            SeqRecord(rec2.seq, id=f"protein|name={rec2.id}", description="")
        ]

        with open(out_path, "w") as handle:
            SeqIO.write(pair, handle, "fasta")

        count += 1

    return count


def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", action="store", help="Input fasta path", required=True)
    parser.add_argument("-o", action="store", help="Output directory path", required=True)
    
    args = parser.parse_args()
    
    make_pairwise_fastas(Path(args.i), Path(args.o))


if __name__ == "__main__":
    main()

