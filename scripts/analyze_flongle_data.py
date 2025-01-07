import pysam
import glob
import plotly.graph_objects as go
import pandas as pd
from argparse import ArgumentParser


def get_strand(read):
    return "-" if read.is_reverse else "+"


def get_ccctct(read):
    return read.query_sequence.count("CCCTCT") if read.query_sequence else 0


def count_ct_by_subtracting_motifs(read):
    seq = read.query_sequence
    if not seq:
        return 0
    motifs = ["CCCTCT", "CCCCT", "CCTT", "CCCT", "CTTT"]
    for motif in motifs:
        seq = seq.replace(motif, "")
    return seq.count("CT")


def get_data(crams, minlength=200, maxlength=1000):
    records = []
    for cram in crams:
        name = (
            cram.replace("masked_rm_map-sminimap2-", "")
            .replace("_hg38s.cram", "")
            .split(".")[0]
        )
        with pysam.AlignmentFile(cram) as f:
            for r in f.fetch("chr15", 34419288, 34419527):
                records.append(
                    (
                        name,
                        get_strand(r),
                        r.query_length,
                        get_ccctct(r),
                        count_ct_by_subtracting_motifs(r),
                    )
                )
    df = pd.DataFrame(
        records, columns=["sample", "strand", "length", "CCCTCT count", "CT count"]
    )
    df = df[(df["length"] > minlength) & (df["length"] < maxlength)]
    df["CCCTCT fraction"] = df["CCCTCT count"] / df["length"]
    return df


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
    fig_ridges.update_layout(xaxis_showgrid=False, xaxis_zeroline=False)
    fig_ridges.update_layout(
        plot_bgcolor="white", xaxis_zerolinecolor="black", yaxis_zerolinecolor="black"
    )
    fig_ridges.update_xaxes(
        showline=True, linewidth=1, linecolor="black", mirror=True, title="Read length"
    )
    fig_ridges.update_yaxes(showline=True, linewidth=1, linecolor="black", mirror=True)
    fig_ridges.update_layout(showlegend=False)
    return fig_ridges


def scatter_plot(df, motif="ccctct fraction"):
    fig_scatter = go.Figure()
    colors = ["blue", "red", "green", "purple", "orange"]
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

def main():
    args = get_args()
    df = get_data(args.crams, minlength=args.minlength, maxlength=args.maxlength)
    with open(args.output, "w") as f:
        f.write(ridges_plot(df).to_html())
        f.write(scatter_plot(df, motif="CCCTCT count").to_html())
        f.write(scatter_plot(df, motif="CT count").to_html())
        f.write(scatter_motifs(df).to_html())

def get_args():
    parser = ArgumentParser()
    parser.add_argument("--crams", nargs="+", help="Input CRAM files")
    parser.add_argument("--minlength", type=int, default=250, help="Minimum read length")
    parser.add_argument("--maxlength", type=int, default=1000, help="Maximum read length")
    parser.add_argument("--output", default="read_lengths.html", help="Output file")
    return parser.parse_args()


if __name__ == "__main__":
    main()