#!/usr/bin/env python
import pandas as pd
from argparse import ArgumentParser
from os.path import basename

def main():
    args = get_args()
    cramino = pd.concat(
        [pd.read_csv(f, sep="\t").set_index("File name") for f in args.input],
        axis="columns",
    )
    cramino.columns = [
        c.replace(".cram", "").replace("_Prom", "") for c in cramino.columns
    ]
    cramino = cramino.transpose().reset_index(names="identifier")
    cramino["name"] = args.names
    cramino.to_csv(args.output, sep="\t", index=False)


def get_args():
    parser = ArgumentParser("")
    parser.add_argument("-i", "--input", nargs="+", help="input files")
    parser.add_argument("-n", "--names", nargs="+", help="names")
    parser.add_argument("-o", "--output", help="output files")
    args = parser.parse_args()
    if args.names:
        if not len(args.input) == len(args.names):
            raise Exception("Expected same number of input files and names")
    else:
        args.names = [basename(f).replace('.cramino', '') for f in args.input]
    return args


if __name__ == "__main__":
    main()
