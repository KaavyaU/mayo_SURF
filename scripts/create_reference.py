import argparse  
from pathlib import Path

'''
USAGE: python3 create_reference.py -i [design fasta path] -o [references fasta path] -s1 [string1] -s2 [string2] -s3 [string3]

string1, string2, string3 are optional parameters (defaults specified below)
string1: TAAACTGCTGACTGGTGAAGGTGGAGGAGGATCAGGAGGAGGTGGTTCAGGTGGTGGAGGCAGC
string2: CCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCTATTATTACCAAGCGGGCAGTGAGCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGCGGATAACAATTTCATGAAAAAGAGGAGAAAGGATCTATGAAGAAGACGGCCATTGCCATTGCGGTTGCGCTGGCCGGATTTGCAACTGTGGCACAGGCA
string3: GGAGGAGGTGGTAGTGGAGGTGGAGGTTCAGGTGGTGGTGGATCTCTGCTGACCTTAGCCAGCCGTCAGCAGCTGATTGATT
'''

def fasta_to_dict(fasta):
    '''Convert a fasta file to a dictionary.'''
    # Create a dictionary to store fasta contents
    records = {}
    header = None
    seq_parts = []
    # Open the fasta file to read the contents
    with open(fasta) as f:
        for line in f:
            # Get rid of any whitespace around the header
            line = line.strip()
            if not line:
                continue
            # If the line is a header, store it and the corresponding sequence in the dictionary
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records[header] = "".join(seq_parts)
    # Return the python dictionary
    return records


def dict_to_fasta(records, fasta):
    '''Write the contents of a dictionary to a fasta file.'''
    # Open a fasta file for writing
    with open(fasta, "w") as f:
        for header, seq in records.items():
            # Write the header to the fasta file first
            f.write(f">{header}\n")
            # Write the corresponding sequence to the fasta file next
            f.write(f"{seq}\n")


def make_reference_fasta(designs, references, string1, string2, string3):
    '''Calls helper functions to create a reference fasta file for the designs.'''
    # Parse contents of designs fasta file
    design_seqs = fasta_to_dict(designs)
    keys = list(design_seqs.keys())
    # Iterate through every combination of designs and store the resulting references 
    reference_dict = {}
    for design1 in keys:
        for design2 in keys:
            reference_seq = string1 + design_seqs[design1] + string2 + design_seqs[design2] + string3
            reference_dict[f"{design1}_{design2}"] = reference_seq
    # Create a new fasta file for the final design references
    dict_to_fasta(reference_dict, references)



def main():
    # Initialize the command line argument parser
    parser = argparse.ArgumentParser()

    # Add final designs fasta path (input), references fasta path (output), and string arguments
    parser.add_argument("-i", "--input", action="store", help="Final designs input fasta file path", required=True)
    parser.add_argument("-o", "--output", action="store", help="Final references output fasta file path", required=True)
    parser.add_argument("-s1", "--string1", action="store", help="pre-Dropin1 string", default="TAAACTGCTGACTGGTGAAGGTGGAGGAGGATCAGGAGGAGGTGGTTCAGGTGGTGGAGGCAGC")
    parser.add_argument("-s2", "--string2", action="store", help="between drop-ins string", default="CCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCTATTATTACCAAGCGGGCAGTGAGCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGCGGATAACAATTTCATGAAAAAGAGGAGAAAGGATCTATGAAGAAGACGGCCATTGCCATTGCGGTTGCGCTGGCCGGATTTGCAACTGTGGCACAGGCA")
    parser.add_argument("-s3", "--string3", action="store", help="post-Dropin2 string", default="GGAGGAGGTGGTAGTGGAGGTGGAGGTTCAGGTGGTGGTGGATCTCTGCTGACCTTAGCCAGCCGTCAGCAGCTGATTGATT")

    # Parse command line arguments
    args = parser.parse_args()

    # Store arguments as python variables
    designs = Path(args.input)
    references = Path(args.output)
    string1 = args.string1
    string2 = args.string2
    string3 = args.string3

    # Call reference fasta function
    make_reference_fasta(designs, references, string1, string2, string3)


if __name__ == "__main__":
    main()
