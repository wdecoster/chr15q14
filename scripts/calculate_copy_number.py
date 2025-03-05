import numpy as np
import gzip
from argparse import ArgumentParser
from os.path import basename

# Calculate copy number in a region
# This script takes a mosdepth .regions.bed.gz file and a region and returns the median copy number in that region
# It also takes a control mosdepth .regions.bed.gz file to normalise the copy number
# The region is specified as a string in the format "chr:start-end"
# Note that the mosdepth files are supposed to have only 1 chromosome
# The output is a single number, the normalised copy number in the region


class Region(object):
    def __init__(self, region, fasta=None):
        self.chromosome, interval = region.replace(",", "").split(":")
        self.begin, self.end = [int(i) for i in interval.split("-")]


def main():
    args = get_args()

    norm_window = Region(args.norm_window)
    call_window = Region(args.call_window)

    norm_factor = get_normalization_factor(args.control, norm_window)
    copy_number = get_copy_number(args.mosdepth, norm_factor, window=call_window)
    print(f"{basename(args.mosdepth).replace('.regions.bed.gz', '')}\t{copy_number}")


def get_normalization_factor(bed, norm_window):
    return np.median(
        np.array(
            [
                float(l.split("\t")[3])
                for l in gzip.open(bed, "rt")
                if norm_window.begin < int(l.split("\t")[2]) < norm_window.end
            ]
        )
    )


def get_copy_number(bed, norm_factor, window):
    return np.median(
        np.array(
            [
                float(l.split("\t")[3]) / norm_factor
                for l in gzip.open(bed, "rt")
                if window.begin < int(l.split("\t")[2]) < window.end
            ]
        )
    )


def get_args():
    parser = ArgumentParser(description="calculate copy number in an interval")
    parser.add_argument(
        "-i", "--mosdepth", help="mosdepth .regions.bed.gz file", required=True
    )
    parser.add_argument(
        "--control",
        help="mosdepth .regions.bed.gz file to use to normalise",
        required=True,
    )
    parser.add_argument(
        "-c", "--call_window", help="region to call copy number in", required=True
    )
    parser.add_argument(
        "-n", "--norm_window", help="region on which to normalise", required=True
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
