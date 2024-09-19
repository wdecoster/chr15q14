# this script will parse a VCF, look at the ALT alleles, and determine the length of the CT stretch
# this first approach will be to subtract all other known motifs from the sequence

from argparse import ArgumentParser
import pandas as pd


def main():
    args = get_args()
    df = pd.read_csv(args.summary, sep="\t")
    df["CT_dimer_count"] = df["sequence"].apply(count_ct_by_subtracting_motifs)
    print(df.to_csv(sep="\t", index=False))


def count_ct_by_subtracting_motifs(seq):
    motifs = ["CCCTCT", "CCCCT", "CCTT", "CCCT", "CTTT"]
    for motif in motifs:
        seq = seq.replace(motif, "")
    return seq.count("CT")


def get_args():
    parser = ArgumentParser(
        description="Determine the length of the CT stretch in a VCF"
    )
    parser.add_argument("summary", help="Summary file to parse", type=str)
    parser.add_argument(
        "-o", "--output", help="output file", type=str, default="ct-stretch.tsv"
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
