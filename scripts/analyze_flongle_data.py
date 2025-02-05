import pysam
import plotly.graph_objects as go
import pandas as pd
from argparse import ArgumentParser
import sys
from itertools import cycle
import pypoars
import tqdm
import concurrent.futures
import os
from fuzzysearch import find_near_matches


class Genotype(object):
    def __init__(self, individual, num_reads, reads_initially, seq):
        self.individual = individual
        self.num_reads = num_reads
        self.reads_initially = reads_initially
        self.seq = seq if seq else "-"
        self.ct_dimer_count = self.count_ct_dimer() if seq else 0
        self.hexamer_count = seq.count("CCCTCT") if seq else 0
        self.length = len(seq) if seq else 0

    def __repr__(self):
        return f"{self.individual}\t{self.num_reads}\t{self.reads_initially}\t{self.length}\t{self.ct_dimer_count}\t{self.hexamer_count}\t{self.seq}"

    def count_ct_dimer(self):
        _seq = self.seq
        for motif in ["CCCTCT", "CCCCT", "CCTTT", "CCTT", "CCCT", "CTTT"]:
            _seq = _seq.replace(motif, "-")
        return _seq.count("CT")


def main():
    args = get_args()
    df = get_data(args)
    genotypes = genotype_samples(df, args)
    for genotype in sorted(genotypes, key=lambda x: x.ct_dimer_count):
        print(genotype)
    if not args.noplot:
        df = df[df["length"].between(args.minlength, args.maxlength)]
        with open(args.output, "w") as out:
            out.write(plot_genotypes(genotypes).to_html())
            # out.write(ridges_plot(df).to_html())
            out.write(scatter_plot(df, motif="CCCTCT count").to_html())
            out.write(scatter_plot(df, motif="CT count").to_html())
            out.write(scatter_motifs(df).to_html())


def genotype_samples(df, args):
    num_samples = len(df["sample"].unique())
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        genotypes = list(
            tqdm.tqdm(
                executor.map(genotype_sample, next_sample(df), [args] * num_samples),
                total=num_samples,
            )
        )
    print(
        "Individual\tNum reads expanded\tNum reads initially\tLength\tCT dimer count\tCCCTCT hexamer count\tConsensus sequence"
    )
    return genotypes


def next_sample(df):
    for sample in df["sample"].unique():
        yield df[df["sample"] == sample]


def genotype_sample(df_sample, args, min_reads=100):
    reads_before = len(df_sample)
    df_sample = df_sample[df_sample["length"].between(args.minlength, args.maxlength)]

    name = df_sample["sample"].iloc[0]
    if (num_reads := len(df_sample)) < min_reads:
        return Genotype(name, num_reads, reads_before, None)
    else:
        # take min_reads number of reads, sorted by length
        seqs = df_sample.sort_values(by="length")["seq"].to_list()[:num_reads]
        # use those min_reads number of longest reads to create a consensus sequence
        consensus = pypoars.poa_consensus(seqs)
        return Genotype(name, num_reads, reads_before, consensus)


def get_data(args):
    records = []
    for cram in args.crams:
        name = (
            os.path.basename(cram)
            .replace("masked_rm_map-sminimap2-", "")
            .replace("_hg38s.cram", "")
            .split(".")[0]
            .replace("_v7", "")
        )
        with pysam.AlignmentFile(cram) as f:
            for r in f.fetch("chr15", 34419288, 34419527):
                records.append(parse_read(r, args.full, name))

    df = pd.DataFrame(
        records,
        columns=["sample", "strand", "length", "seq", "CCCTCT count", "CT count"],
    )
    return df


def parse_read(read, full, name):
    seq = read.query_sequence if full else non_ref_bases(read)
    if not seq:
        return name, get_strand(read), 0, None, 0, 0
    fragment_length = read.query_length if full else len(seq)
    return (
        name,
        get_strand(read),
        fragment_length,
        seq,
        get_ccctct(seq),
        count_ct_by_subtracting_motifs(seq),
    )


def non_ref_bases(read, minlength=50):
    """
    This function slices out the bases from a read that do not match the reference genome, by parsing the CIGAR string.
    Only cigar operations longer than minlength are considered, and only insertions and softclips are returned.
    The script iterates over the cigar string, while moving the cursor in the read sequence
    """
    if not read.query_sequence:
        return ""
    non_ref = ""
    read_position = 0
    for operation, length in read.cigartuples:
        if operation in [0, 7, 8]:
            # operation 0 is match (M), 7 is match with sequence (=), 8 is alignment match (X)
            read_position += length
        elif operation == 3:
            # operation 3 is refskip (N), this shouldn't happen
            sys.stderr.write("Warning: unexpected refskip cigar operation in read\n")
            continue
        elif operation == 2:
            # operation 2 is deletion (D). This script does not care about repeat contractions.
            continue
        elif operation in [1, 4]:
            # operation 4 is softclip (S), operation 1 is insertion (I)
            if length >= minlength:
                non_ref += read.query_sequence[read_position : read_position + length]
            read_position += length
    # attempt to trim off the reference sequences that may be caugth in the softclipped region
    # however, it is not likely that perfect matches are found for every read
    # therefore, using a fuzzy search allowing for a few mismatches (~5%) is used
    right_seq = "GAGACGGAGTTTCTCTCTTGTTGCCCAGGCTGGAGTGCATGTTGCTGTGCACTTTGAGGGCAGGAACTG"
    matches = find_near_matches(right_seq, non_ref, max_l_dist=4)
    if len(matches) > 1:
        return None
    if len(matches) == 1:
        non_ref = non_ref[: matches[0].start]

    left_seq = "TTATCAGGGCCTCTCTTCGCAGGCAGTGGGGCCTCATCCACAACCCTGGAAAAGAACTGGAAAGCGTTGCTCAGCCAGGTACGGAGGGCAGGGCCATGTGGGACTCCCGTCTCCAGGCCCCCTCTCCCCAGCTCCCG"
    matches = find_near_matches(left_seq, non_ref, max_l_dist=7)
    if len(matches) > 1:
        return None
    if len(matches) == 1:
        non_ref = non_ref[matches[0].end :]
    return non_ref


def get_strand(read):
    return "-" if read.is_reverse else "+"


def get_ccctct(seq):
    return seq.count("CCCTCT") if seq else 0


def count_ct_by_subtracting_motifs(seq):
    if not seq:
        return 0
    for motif in ["CCCTCT", "CCCCT", "CCTTT", "CCTT", "CCCT", "CTTT"]:
        seq = seq.replace(motif, "-")
    return seq.count("CT")


def ridges_plot(df):
    fig_ridges = go.Figure()
    for sample in df["sample"].unique():
        for strand in ["+", "-"]:
            df_strand = df[(df["sample"] == sample) & (df["strand"] == strand)]
            fig_ridges.add_trace(
                go.Violin(
                    x=df_strand["length"],
                    line_color="blue" if strand == "+" else "red",
                    name=f"{sample}_{strand}",
                )
            )
    fig_ridges.update_traces(
        orientation="h", side="positive", width=3, points=False, spanmode="hard"
    )
    fig_ridges.update_layout(
        plot_bgcolor="white",
        xaxis_zerolinecolor="black",
        yaxis_zerolinecolor="black",
        showlegend=False,
    )
    fig_ridges.update_xaxes(
        showline=True, linewidth=1, linecolor="black", mirror=True, title="Read length"
    )
    fig_ridges.update_yaxes(showline=True, linewidth=1, linecolor="black", mirror=True)
    return fig_ridges


def scatter_plot(df, motif):
    fig_scatter = go.Figure()
    colors = ["blue", "red", "green", "purple", "orange"]
    if len(df["sample"].unique()) > len(colors):
        sys.stderr.write(
            "Warning: more samples than colors, some samples will have the same color.\n"
        )
        colors = cycle(colors)
    for sample, color in zip(df["sample"].unique(), colors):
        df_sample = df[df["sample"] == sample]
        fig_scatter.add_trace(
            go.Scatter(
                x=df_sample[motif],
                y=df_sample["length"],
                mode="markers",
                name=sample,
                marker=dict(size=3, color=color),
                hovertext=df_sample["sample"],
            )
        )
    fig_scatter.update_layout(
        plot_bgcolor="white", xaxis_zerolinecolor="black", yaxis_zerolinecolor="black"
    )
    fig_scatter.update_xaxes(
        showline=True,
        linewidth=1,
        linecolor="black",
        mirror=True,
        title=motif,
        rangemode="nonnegative",
    )
    fig_scatter.update_yaxes(
        showline=True,
        linewidth=1,
        linecolor="black",
        mirror=True,
        title="Read length",
        rangemode="nonnegative",
    )
    return fig_scatter


def scatter_motifs(df):
    """
    Make a scatter plot of CTs vs CCCTCTs, colored by read length, with shape indicating the sample
    """
    df = df.sort_values("length")
    fig_scatter = go.Figure()
    min_length = df["length"].min()
    # for max_length, use the 95th percentile to avoid outliers
    max_length = df["length"].quantile(0.95)
    for sample in df["sample"].unique():
        df_sample = df[df["sample"] == sample]
        fig_scatter.add_trace(
            go.Scatter(
                x=df_sample["CCCTCT count"],
                y=df_sample["CT count"],
                mode="markers",
                name=sample,
                marker=dict(
                    size=3,
                    color=df_sample["length"],
                    colorscale="Reds",
                    cmin=min_length,
                    cmax=max_length,
                ),
                hovertext=df_sample["sample"]
                + "<br>Length: "
                + df_sample["length"].astype(str),
            )
        )
    fig_scatter.update_layout(
        coloraxis=dict(
            colorscale="Reds",
            cmin=min_length,
            cmax=max_length,
            colorbar=dict(title="Read length"),
        )
    )
    fig_scatter.update_layout(
        plot_bgcolor="white", xaxis_zerolinecolor="black", yaxis_zerolinecolor="black"
    )
    fig_scatter.update_xaxes(
        showline=True,
        linewidth=1,
        linecolor="black",
        mirror=True,
        title="CCCTCT count",
        rangemode="nonnegative",
    )
    fig_scatter.update_yaxes(
        showline=True,
        linewidth=1,
        linecolor="black",
        mirror=True,
        title="CT count",
        rangemode="nonnegative",
    )
    return fig_scatter


def plot_genotypes(genotypes):
    """
    Make a scatter plot showing for each sample the CT dimer count vs the CCCTCT hexamer count
    """
    fig_scatter = go.Figure(
        go.Scatter(
            x=[genotype.ct_dimer_count for genotype in genotypes],
            y=[genotype.hexamer_count for genotype in genotypes],
            mode="markers",
            text=[genotype.individual for genotype in genotypes],
            marker=dict(size=4, color="black"),
            hovertext=[genotype.individual for genotype in genotypes],
        )
    )
    fig_scatter.update_layout(
        plot_bgcolor="white",
        margin=dict(l=0, r=0, t=50, b=0),
        title="Genotypes from spanning PCR",
        font=dict(size=16),
        height=400,
        width=800,
    )
    fig_scatter.update_xaxes(
        showline=True,
        linewidth=1,
        linecolor="black",
        mirror=True,
        title="CT dimers",
    )
    fig_scatter.update_yaxes(
        showline=True,
        linewidth=1,
        linecolor="black",
        mirror=True,
        title="CCCTCT hexamers",
    )
    fig_scatter.add_vline(x=190, line_width=2, line_dash="dash", line_color="black")
    return fig_scatter


def get_args():
    parser = ArgumentParser()
    parser.add_argument("--crams", nargs="+", help="Input CRAM file(s)")
    parser.add_argument(
        "--minlength", type=int, default=50, help="Minimum fragment length"
    )
    parser.add_argument(
        "--maxlength", type=int, default=1200, help="Maximum fragment length"
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="Use entire read sequence, not just non-reference bases",
    )
    parser.add_argument("--output", default="read_lengths.html", help="Output file")
    parser.add_argument(
        "--threads", type=int, default=4, help="Number of threads to use"
    )
    parser.add_argument("--noplot", action="store_true", help="Don't make plots")
    return parser.parse_args()


if __name__ == "__main__":
    main()
