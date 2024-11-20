"""
This script takes in an overview file created in the workflow, and multiple pairs of lists and output files to split the overview file on the basis of the lists.

the cli looks like:

python {params.script} -i {input} \
--all {params.all} --all_out {output.all} \
--hapA {params.hapA} --hapA_out {output.hapA} \
--hapB {params.hapB} --hapB_out {output.hapB} \
--duplicates {params.duplicates} --duplicates_out {output.duplicates} 2> {log}

where the input is the overview file, and the output files are the split overview files. The lists are the haplotype carriers, relatives, and duplicates, respectively, which are lists of identifiers, not a file.
aSTRonaut also only needs the columns name, sequence and group, so we can drop the rest of the columns.
"""

from argparse import ArgumentParser
import pandas as pd

def main():
    args = get_args()
    df = pd.read_csv(args.input, sep="\t", usecols=["name", "sequence", "group"])
    df["group"] = df["group"].apply(lambda x: "case" if x == "aFTLD-U" else x)
    for l in [args.all, args.hapA, args.hapB]:
        check_valid_list_of_identifiers(df, l)
    df[df["name"].isin(args.all)].to_csv(args.all_out, sep="\t", index=False)
    df[df["name"].isin(args.hapA)].to_csv(args.hapA_out, sep="\t", index=False)
    df[df["name"].isin(args.hapB)].to_csv(args.hapB_out, sep="\t", index=False)

def check_valid_list_of_identifiers(df, list_of_identifiers):
    for identifier in list_of_identifiers:
        if identifier not in df["name"].values:
            raise ValueError(f"{identifier} is not in the overview file")
        


def get_args():
    parser = ArgumentParser("Split the overview file to create tables for aSTRonaut")
    parser.add_argument("-i", "--input", help="Input overview file")
    parser.add_argument("--all", help="List of all samples", nargs="+")
    parser.add_argument("--all_out", help="Output file for all samples")
    parser.add_argument("--hapA", help="List of haplotype A carriers", nargs="+")
    parser.add_argument("--hapA_out", help="Output file for haplotype A carriers")
    parser.add_argument("--hapB", help="List of haplotype B carriers", nargs="+")
    parser.add_argument("--hapB_out", help="Output file for haplotype B carriers")
    return parser.parse_args()


if __name__ == "__main__":
    main()
