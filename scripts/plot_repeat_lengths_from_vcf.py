import pandas as pd
from argparse import ArgumentParser
import plotly.express as px
from cyvcf2 import VCF
from os.path import basename
import sys


def main():
    args = get_args()
    res = []
    alts = []
    for vcf in args.input:
        name = basename(vcf).replace(".vcf.gz", "").replace("_FCX", "")
        try:
            for variant in VCF(vcf):
                for phase in [0, 1]:
                    if variant.INFO.get("SEQS") is None:
                        sys.stderr.write(f"No SEQS in {name}\n")
                        continue
                    # probably lowly covered and/or homozygous
                    if len(variant.INFO.get("SEQS").split(",")) == 1 and phase == 1:
                        continue
                    for seq in variant.INFO.get("SEQS").split(",")[phase].split(":"):
                        res.append(
                            {
                                "name": name,
                                "phase": phase,
                                "length": len(seq),
                            }
                        )
                if variant.INFO.get("OUTLIERS") is not None:
                    for outlier in variant.INFO.get("OUTLIERS").split(","):
                        if (
                            not name == "rr_KUL_LEU1"
                        ):  # this sample has a weird outlier (artefact)
                            res.append(
                                {
                                    "name": name,
                                    "phase": "outlier",
                                    "length": len(outlier),
                                }
                            )
                # also append the length of the ALT allele
                for alt_allele in variant.ALT:
                    alts.append(
                        {
                            "name": name,
                            "alt_length": len(alt_allele),
                        }
                    )
        except Exception:
            sys.stderr.write(f"Problem parsing {name}\n\n\n")
            raise
    df = pd.DataFrame(res, columns=["name", "phase", "length"])
    alts = (
        pd.DataFrame(alts, columns=["name", "alt_length"])
        .groupby("name")
        .max()
        .reset_index(names="name")
    )
    if args.alphabetically:
        df = df.sort_values("name")
    else:
        # sort the names by the length of the alternative allele
        # put the largest alternative allele first
        df = df.merge(alts, on="name", how="left").sort_values(
            "alt_length", ascending=False
        )
    if args.sampleinfo:
        sampleinfo = pd.read_excel(
            args.sampleinfo, usecols=["Gentli_ID", "PathDx1"]
        ).rename(columns={"Gentli_ID": "name", "PathDx1": "Group"})
        df = df.merge(sampleinfo, on="name", how="left").fillna("unknown/misc")
    else:
        df["Group"] = "unknown"

    # only keep names for which the alt is long enough
    long_enough = alts.loc[alts["alt_length"] >= args.minlen, "name"]
    df = df[df["name"].isin(long_enough)]
    orders = df["name"].drop_duplicates()

    # all groups that are not aFTLD-U are in-house controls
    df.loc[df["Group"] != "aFTLD-U", "Group"] = "in-house control"

    fig = px.strip(
        df,
        y="name",
        x="length",
        color="Group",
        stripmode="overlay",
        color_discrete_map={
            "aFTLD-U": "red",
            "in-house control": "black",
        },
        orientation="h",
    )
    fig.update_traces(marker=dict(size=4))
    fig.update_layout(
        title=args.title,
        yaxis_title="Individuals",
        xaxis_title="Expansion length",
        plot_bgcolor="white",
        height=800,
        width=800,
        font=dict(size=20),
        legend=dict(yanchor="bottom", y=0.01, xanchor="right", x=0.99),
    )
    fig.update_xaxes(
        showgrid=True,
        gridcolor="lightgray",
        showline=True,
        linewidth=2,
        linecolor="black",
        mirror=True,
    )
    fig.update_yaxes(
        categoryorder="array",
        categoryarray=orders[::-1],
        showticklabels=False,
        showline=True,
        linewidth=2,
        linecolor="black",
        mirror=True,
    )
    if args.mark_hexamers:
        maxlengths = df.groupby(["name"]).max()
        # add a star symbol as annotation to 3 specific lines: the one of rr_EMC1999_019, rr_NA16_299 and rr_NA11_305
        for name in ["rr_EMC1999_019", "rr_NA16_299", "rr_NA11_305", "rr_ERA09_35"]:
            if name in maxlengths.index:
                fig.add_annotation(
                    x=maxlengths.loc[name, "length"] + 200,
                    y=name,
                    text="*",
                    showarrow=False,
                    font=dict(size=16),
                )
    with open(args.output, "w") as output:
        output.write(fig.to_html(include_plotlyjs="cdn"))


def get_args():
    parser = ArgumentParser(description="Plot repeat lengths from STRdust")
    parser.add_argument("-i", "--input", help="VCF files from STRdust", nargs="+")
    parser.add_argument(
        "-o", "--output", help="output plot", default="repeat_lengths.html"
    )
    parser.add_argument(
        "-m", "--minlen", help="minimal length to plot", default=20, type=int
    )
    parser.add_argument(
        "--alphabetically", help="Sort the samples alphabetically", action="store_true"
    )
    parser.add_argument("--sampleinfo", help="excel file with sample information")
    parser.add_argument(
        "--mark_hexamers",
        help="Mark the hexamers in the plot",
        action="store_true",
    )
    parser.add_argument(
        "--black_and_white",
        help="Use symbols in addition of colors",
        action="store_true",
    )
    parser.add_argument(
        "--title", help="Title of the plot", default="Repeat length per sequenced read"
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
