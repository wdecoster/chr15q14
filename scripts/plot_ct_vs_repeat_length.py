from random import random
import pandas as pd
from argparse import ArgumentParser
import plotly.express as px
from cyvcf2 import VCF
from os.path import basename
import sys


def main():
    args = get_args()
    res = []
    for vcf in args.vcfs:
        name = basename(vcf).replace(".vcf.gz", "").replace("_FCX", "")
        try:
            for variant in VCF(vcf):
                if args.showboth:
                    for alt, sup, stddev in zip(
                        variant.ALT,
                        variant.format("SUP")[0],
                        variant.INFO.get("STDEV"),
                    ):
                        if sup >= args.support:
                            res.append(
                                {
                                    "name": name,
                                    "variant": f"{variant.CHROM}:{variant.POS}",
                                    "length": len(alt),
                                    "%CT": kmer_fraction(alt, kmer="CT", jitter=True),
                                    "%CCTT": kmer_fraction(alt, kmer="CCTT"),
                                    "%CCCTCT": kmer_fraction(alt, kmer="CCCTCT"),
                                    "%AT": kmer_fraction(alt, kmer="AT"),
                                    "stddev": stddev,
                                    "support": sup,
                                    "sequence": alt,
                                }
                            )
                else:
                    # if not showing both alleles, choose one
                    alts = [
                        (alt, sup, stddev)
                        for alt, sup, stddev in zip(
                            variant.ALT,
                            variant.format("SUP")[0],
                            variant.INFO.get("STDEV"),
                        )
                        if sup >= args.support
                    ]
                    if alts:
                        # if there are multiple alts, choose the on with the highest CT count, and if there are multiple with the same count, choose the longest
                        sys.stderr.write(f"Warning: Multiple possible alts for {name} at {variant.CHROM}:{variant.POS}\n")
                        alt, sup, stddev = max(alts, key=lambda x: (count_ct_by_subtracting_motifs(x[0]), len(x[0])))
                        res.append(
                            {
                                "name": name,
                                "variant": f"{variant.CHROM}:{variant.POS}",
                                "length": len(alt),
                                "%CT": kmer_fraction(alt, kmer="CT", jitter=True),
                                "%CCTT": kmer_fraction(alt, kmer="CCTT"),
                                "%CCCTCT": kmer_fraction(alt, kmer="CCCTCT"),
                                "%AT": kmer_fraction(alt, kmer="AT"),
                                "CT_dimer_count": count_ct_by_subtracting_motifs(alt),
                                "stddev": stddev,
                                "support": sup,
                                "sequence": alt,
                            }
                        )

        except Exception:
            sys.stderr.write(f"Problem parsing {name}\n\n\n")
            raise
    df = pd.DataFrame(
        res,
        columns=[
            "name",
            "variant",
            "length",
            "%CT",
            "%CCTT",
            "%CCCTCT",
            "%AT",
            "CT_dimer_count",
            "stddev",
            "support",
            "sequence",
        ],
    )

    if args.sample:
        df["highlight"] = df["name"].isin(args.sample)
    else:
        df["highlight"] = False
    if args.sampleinfo:
        sampleinfo = pd.read_csv(
            args.sampleinfo,
            sep="\t",
            usecols=["individual", "cohort", "sex", "haplotype"],
        ).rename(columns={"cohort": "group", "individual": "name"})
        df = df.merge(sampleinfo, on="name", how="left").drop_duplicates()
    if args.copy_number:
        copy_number = pd.read_csv(args.copy_number, sep="\t", header=None, names=["name", "copy number"])
        copy_number["copy number"] = copy_number["copy number"].round(2)
        df = df.merge(copy_number, on="name", how="left").drop_duplicates()
    df.to_csv(args.overview, index=False, sep="\t")
    # filter on the minimum repeat length for plotting
    df = df.loc[df["length"] >= args.min_length]
    with open(args.output, "w") as output:
        for locus in df["variant"].unique():
            title = (
                f"Consensus repeat length vs %CT at {locus}"
                if df["variant"].nunique() > 1
                else "Consensus repeat length vs %CT"
            )
            fig, upper_limit = make_scatter_plot(
                df.loc[df["variant"] == locus], title, args
            )

            fig.write_html(
                output,
                include_plotlyjs="cdn",
                full_html=True if output.tell() == 0 else False,
            )
            if args.haplotypes:
                fig, _ = make_scatter_plot(
                    df.loc[
                        (df["variant"] == locus) & (df["haplotype"] == "major")
                    ].drop(columns=["haplotype"]),
                    title + " for haplotype A",
                    args,
                    upper_limit=upper_limit,
                )
                fig.write_html(
                    output,
                    include_plotlyjs="cdn",
                    full_html=False,
                )
                fig, _ = make_scatter_plot(
                    df.loc[
                        (df["variant"] == locus) & (df["haplotype"] == "minor")
                    ].drop(columns=["haplotype"]),
                    title + " for haplotype B",
                    args,
                    upper_limit=upper_limit,
                )
                fig.write_html(
                    output,
                    include_plotlyjs="cdn",
                    full_html=False,
                )

def count_ct_by_subtracting_motifs(seq):
    motifs = ["CCCTCT", "CCCCT", "CCTT", "CCCT", "CTTT"]
    for motif in motifs:
        seq = seq.replace(motif, "")
    return seq.count("CT")


def make_scatter_plot(df, title, args, upper_limit=None):
    if upper_limit is None:
        # round up to the nearest 500
        upper_limit = (df["length"].max() // 500 + 1) * 500
    # color the plot based on the groups column if it exists
    fig = px.scatter(
        df,
        x="%CT",
        y="length",
        color="group" if "group" in df.columns else "black",
        color_discrete_map={
            "1000G": "teal",
            "aFTLD-U": "red",
            "in-house control": "black",
        },
        #symbol="haplotype" if "haplotype" in df.columns else None,
        #symbol_sequence=["x", "cross", "circle"],
        hover_data=[
            "name",
            "%CCTT",
            "%CCCTCT",
            "%AT",
            "support",
            "stddev",
            "sex",
            "copy number",
            "haplotype" if "haplotype" in df.columns else None,
        ],
        title=title,
        error_y="stddev" if args.stddev else None,
    )
    if args.xline and args.yline:
        fig.add_shape(
            type="line",
            x0=args.xline,
            y0=args.yline,
            x1=1,
            y1=args.yline,
            line_dash="dot",
            line=dict(color="black", width=1),
        )
        fig.add_shape(
            type="line",
            x0=args.xline,
            y0=args.yline,
            x1=args.xline,
            y1=upper_limit,
            line_dash="dot",
            line=dict(color="black", width=1),
        )
    elif args.xline:
        fig.add_hline(y=args.xline, line_dash="dot", line_color="black", line_width=1)
    elif args.yline:
        fig.add_vline(x=args.yline, line_dash="dot", line_color="black", line_width=1)
    if args.arrow:
        for sample, color in zip(
            args.arrow.split(","), ["blue", "green", "orange", "purple"]
        ):
            if sample in df["name"].values:
                fig.add_annotation(
                    x=df.loc[df["name"] == sample, "%CT"].values[0] + 0.01,
                    y=df.loc[df["name"] == sample, "length"].values[0]
                    - 0.02 * upper_limit,
                    showarrow=True,
                    arrowhead=2,
                    arrowcolor=color,
                    arrowwidth=0.8,
                    ax=20,
                    ay=20,
                )
    fig.update_traces(marker=dict(size=9, opacity=0.7))
    fig.update_layout(
        plot_bgcolor="white",
        font=dict(size=20),
        legend=dict(
            title="",
            itemsizing="constant",
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.05,
        ),
        width=1000,
        height=600,
    )
    fig.update_xaxes(
        showline=True,
        linewidth=2,
        linecolor="black",
        mirror=True,
        range=[0, 1],
        title_text="%CT",
    )

    fig.update_yaxes(
        showline=True,
        linewidth=2,
        linecolor="black",
        mirror=True,
        range=[0, upper_limit],
        title_text="Repeat length",
    )

    # I don't like how the legend now looks, with both the group and haplotype in the same legend, resulting in 9 items
    # I want to have the group in the legend, and the haplotype as a separate legend
    # I want a separate legend for the group, and a separate legend for the haplotype
    if "haplotype" in df.columns:
        symbol_map = {"major": "x", "minor": "cross", "none": "circle"}

        # for each trace update marker symbol to list of symbols that correspond to haplotype

        for trace in fig.data:
            trace_haplotypes = df[df["group"] == trace.name]["haplotype"]
            trace.update(marker_symbol=[symbol_map[hap] for hap in trace_haplotypes])

        # add a legend for the haplotype, with the correct symbols
        fig.for_each_trace(lambda trace: trace.update(showlegend=False))
        import plotly.graph_objects as go
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(symbol="circle", size=10, color="teal"),
                name="1000G",
                legendgroup="group",
                legendgrouptitle=dict(text="Group")
            )
        )
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(symbol="circle", size=10, color="red"),
                name="aFTLD-U",
                legendgroup="group",
                legendgrouptitle=dict(text="Group"),
            )
        )

        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(symbol="circle", size=10, color="black"),
                name="in-house control",
                legendgroup="group",
                legendgrouptitle=dict(text="Group"),
            )
        )
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(symbol="x", size=10, color="grey"),
                name="Haplotype A",
                legendgroup="hap",
                legendgrouptitle=dict(text="Haplotype")
            )
        )
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(symbol="cross", size=10, color="grey"),
                name="Haplotype B",
                legendgroup="hap",
                legendgrouptitle=dict(text="Haplotype")
            )
        )
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(symbol="circle", size=10, color="grey"),
                name="No associated haplotype",
                legendgroup="hap",
                legendgrouptitle=dict(text="Haplotype")
            )
        )


    if df["length"].max() > upper_limit:
        sys.stderr.write(
            "Warning: Some samples have repeat lengths longer than the current y-axis limit\n"
        )
    return fig, upper_limit


def kmer_fraction(seq, kmer, jitter=False):
    """count the frequency of either CT or GA in a sequence
    to avoid overlapping points, add a small random number to the frequency
    round the result to 3 decimals
    """
    # get the reverse complement of the kmer
    kmer_rc = "".join(
        {"A": "T", "T": "A", "C": "G", "G": "C"}[base] for base in kmer.upper()
    )
    if jitter:
        # get a random number between -0.001 and 0.001
        rand = (random() - 0.5) / 500
    else:
        rand = 0
    return round(
        (seq.count(kmer) + seq.count(kmer_rc)) / (len(seq) / len(kmer)) + rand, 3
    )


def get_args():
    parser = ArgumentParser()
    parser.add_argument("vcfs", nargs="+", help="VCF files to parse")
    parser.add_argument(
        "-o",
        "--output",
        help="Output file",
        default="scatter_lengths_vs_composition.html",
    )
    parser.add_argument("--sample", help="Sample(s) to highlight", nargs="*")
    parser.add_argument(
        "--support", help="Minimum read support for a call", type=int, default=2
    )
    parser.add_argument(
        "-m", "--min-length", help="Minimum repeat length", type=int, default=20
    )
    parser.add_argument(
        "--copy_number", help="File with copy number information to add to the table")
    parser.add_argument(
        "-s", "--sampleinfo", help="Sampleinfo file to link samples to groups"
    )
    parser.add_argument("--showboth", help="show both haplotypes", action="store_true")
    parser.add_argument(
        "--haplotypes", help="Show haplotypes separately", action="store_true"
    )
    parser.add_argument(
        "--stddev", help="Show repeat length standard deviation", action="store_true"
    )
    parser.add_argument("--xline", help="Add a line at this x value", type=float)
    parser.add_argument("--yline", help="Add a line at this y value", type=float)
    parser.add_argument(
        "--arrow", help="Add an arrow at specific comma-separated samples"
    )
    parser.add_argument(
        "--overview",
        help="Output file for the overview table",
        default="repeat_length_overview.tsv",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
