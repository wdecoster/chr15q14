import pysam
import plotly.graph_objects as go
import pandas as pd
from argparse import ArgumentParser
import sys
from itertools import cycle


class Genotype(object):
    def __init__(self, individual, num_reads, seq):
        self.individual = individual
        self.num_reads = num_reads
        self.seq = seq


def main():
    args = get_args()
    df = get_data(args)
    with open(args.output, "w") as out:
        out.write(ridges_plot(df).to_html())
        out.write(scatter_plot(df, motif="CCCTCT count").to_html())
        out.write(scatter_plot(df, motif="CT count").to_html())
        out.write(scatter_motifs(df).to_html())
    genotype(df)

def genotype(df, min_reads=200):
    genotypes = []
    for sample in df["sample"].unique():
        df_sample = df[df["sample"] == sample]
        if (num_reads := len(df_sample)) < min_reads:
            genotypes.append(Genotype(sample, num_reads, None))
        else:
            # take min_reads number of reads, sorted by length
            seqs = df_sample.sort_values(by="length")["seq"].to_list()[:num_reads]
            # use those min_reads number of longest reads to create a consensus sequence, excluding outliers

def get_data(args):
    records = []
    for cram in args.crams:
        name = (
            cram.replace("masked_rm_map-sminimap2-", "")
            .replace("_hg38s.cram", "")
            .split(".")[0]
            .replace("_v7", "")
        )
        with pysam.AlignmentFile(cram) as f:
            for r in f.fetch("chr15", 34419288, 34419527):
                records.append(parse_read(r, args.full, name))

    df = pd.DataFrame(
        records, columns=["sample", "strand", "length", "seq", "CCCTCT count", "CT count"]
    )
    df = df[df["length"].between(args.minlength, args.maxlength)]
    return df


def parse_read(read, full, name):
    seq = read.query_sequence if full else non_ref_bases(read)
    fragment_length = read.query_length if full else len(seq)
    return name, get_strand(read), fragment_length, seq, get_ccctct(seq), count_ct_by_subtracting_motifs(seq)

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
        sys.stderr.write("Warning: more samples than colors, some samples will have the same color.\n")
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


def get_args():
    parser = ArgumentParser()
    parser.add_argument("--crams", nargs="+", help="Input CRAM files")
    parser.add_argument("--minlength", type=int, default=50, help="Minimum read length")
    parser.add_argument("--maxlength", type=int, default=1200, help="Maximum read length")
    parser.add_argument("--full", action="store_true", help="Use entire read sequence, not just non-reference bases")
    parser.add_argument("--output", default="read_lengths.html", help="Output file")
    return parser.parse_args()


if __name__ == "__main__":
    main()
