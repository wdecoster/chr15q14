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
                        if len(alt) >= args.min_length and sup >= args.support:
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
                    # if not showing both alleles, show the longer one of the those with sufficient support
                    alts = [
                        (alt, sup, stddev)
                        for alt, sup, stddev in zip(
                            variant.ALT,
                            variant.format("SUP")[0],
                            variant.INFO.get("STDEV"),
                        )
                        if len(alt) >= args.min_length and sup >= args.support
                    ]
                    if alts:
                        alt, sup, stddev = max(alts, key=lambda x: len(x[0]))
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
            "stddev",
            "support",
            "sequence",
        ],
    )

    if args.sample:
        df["highlight"] = df["name"].isin(args.sample)
    else:
        df["highlight"] = False
    if args.groups:
        sampleinfo = pd.read_csv(
            args.groups,
            sep="\t",
            usecols=["individual", "cohort", "sex"],
        ).rename(columns={"cohort": "group", "individual": "name"})
        df = df.merge(sampleinfo, on="name", how="left").sort_values("group")
    with open(args.output, "w") as output:
        for locus in df["variant"].unique():
            title = (
                f"Repeat length {locus}"
                if df["variant"].nunique() > 1
                else "Repeat length"
            )
            fig = make_violin_plot(df.loc[df["variant"] == locus], title, args)

            fig.write_html(
                output,
                include_plotlyjs="cdn",
                full_html=True if output.tell() == 0 else False,
            )


def make_violin_plot(df, title, args):
    # color the plot based on the groups column if it exists
    # change the in-house controls label with a line break
    df["group"] = df["group"].replace("in-house non-aFTLD-U", "in-house<br>non-aFTLD-U")
    fig = px.strip(
        df,
        x="group",
        y="length",
        color_discrete_sequence=["black"],
        hover_data=[
            "name",
            "%CCTT",
            "%CCCTCT",
            "%AT",
            "support",
            "stddev",
            "sex",
            "haplotype" if "haplotype" in df.columns else None,
        ],
        title=title,
    )
    if args.line:
        fig.add_hline(y=args.line, line_dash="dash", line_color="grey")
    # make the dots smaller
    # fig.update_traces(spanmode="hard", marker=dict(size=4, color="black"))
    fig.update_traces(marker=dict(size=4, color="black"), jitter=1.0)
    fig.update_layout(plot_bgcolor="white")
    fig.update_yaxes(title_text="Consensus repeat length")
    fig.update_xaxes(title_text="")
    fig.update_layout(
        font=dict(size=18),
        legend=dict(title="Sample", itemsizing="constant"),
        margin=dict(l=0, r=0, t=50, b=0),
    )
    # show zerolines
    fig.update_xaxes(showline=True, linewidth=2, linecolor="black", mirror=True)
    upper_limit = 2500
    fig.update_yaxes(
        showline=True,
        linewidth=2,
        linecolor="black",
        mirror=True,
        range=[0, upper_limit],
    )
    # reduce the horizontal space
    fig.update_layout(width=400, height=800)
    # print a warning if any of the samples have a repeat length longer than the current y-axis limit
    if df["length"].max() > upper_limit:
        sys.stderr.write(
            "Warning: Some samples have repeat lengths longer than the current y-axis limit\n"
        )
    return fig


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
        default="lengths_violin.html",
    )
    parser.add_argument("--sample", help="Sample(s) to highlight", nargs="*")
    parser.add_argument(
        "--support", help="Minimum read support for a call", type=int, default=3
    )
    parser.add_argument(
        "-m", "--min-length", help="Minimum repeat length", type=int, default=20
    )
    parser.add_argument(
        "-g", "--groups", help="Sampleinfo file to link samples to groups"
    )
    parser.add_argument("--line", help="Height to show a vertical line on", type=int)
    parser.add_argument("--showboth", help="show both haplotypes", action="store_true")
    parser.add_argument(
        "--stddev", help="Show repeat length standard deviation", action="store_true"
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
