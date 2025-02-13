from argparse import ArgumentParser
import pandas as pd
from cyvcf2 import VCF
from os.path import basename
import sys
from plotly.subplots import make_subplots
import plotly.express as px


def main():
    args = get_args()
    args = get_args()
    res = []
    for vcf in args.vcfs:
        name = basename(vcf).replace(".vcf.gz", "").replace("_FCX", "")
        try:
            for variant in VCF(vcf):
                for alt, stdev in zip(variant.ALT, variant.INFO.get("STDEV")):
                    res.append(
                        {
                            "name": name,
                            "variant": f"{variant.CHROM}:{variant.POS}",
                            "length": len(alt),
                            "stdev": stdev,
                        }
                    )
        except Exception:
            sys.stderr.write(f"Problem parsing {name}\n\n\n")
            raise
    df = pd.DataFrame(res, columns=["name", "variant", "length", "stdev"]).merge(
        pd.read_csv(
            args.sample_info,
            sep="\t",
            usecols=["individual", "cohort", "haplotype"],
        ).rename(columns={"cohort": "group", "individual": "name"}),
        on="name",
        how="left",
    )
    df = df[df["length"] >= args.minlen]
    haplotype_alias = {"major": "HapA", "minor": "HapB", "none": "none"}
    df["haplotype"] = df["haplotype"].apply(lambda x: haplotype_alias[x])
    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=(
            "Repeat length vs. standard deviation",
            "Standard deviation per haplotype",
        ),
    )
    scatter = px.scatter(
        df,
        x="length",
        y="stdev",
        color="group",
        title="Repeat length vs. standard deviation",
        labels={
            "length": "Repeat length",
            "stdev": "Standard deviation of repeat length",
            "haplotype": "Haplotype",
        },
    )
    strip = px.strip(
        df,
        x="haplotype",
        y="stdev",
        color="group",
        title="Standard deviation of repeat lengths per group",
        labels={"stdev": "Standard deviation of repeat length", "group": "Group"},
        hover_data=["name"],
    )
    figures = [scatter, strip]
    for i, figure in enumerate(figures, start=1):
        for trace in figure.data:
            fig.add_trace(trace, row=1, col=i)

    names = set()
    fig.for_each_trace(
        lambda trace: (
            trace.update(showlegend=False)
            if (trace.name in names)
            else names.add(trace.name)
        )
    )

    fig.update_traces(
        marker=dict(
            size=10,
            opacity=0.7,
        )
    )
    fig.update_layout(
        plot_bgcolor="white",
        font=dict(size=28),
    )
    # put the legend on the bottom, centered between the two subplots
    fig.update_layout(
        legend=dict(
            yanchor="bottom",
            y=-0.15,
            xanchor="center",
            x=0.5,
            bordercolor="black",
            borderwidth=0.5,
            orientation="h",
        )
    )

    fig.update_xaxes(
        showline=True,
        linewidth=2,
        linecolor="black",
        mirror=True,
    )
    fig.update_yaxes(
        showline=True,
        linewidth=2,
        linecolor="black",
        mirror=True,
    )
    # change axis labels for each subplot
    fig.update_xaxes(title_text="Repeat length", row=1, col=1)
    fig.update_yaxes(title_text="Standard deviation of length", row=1, col=1)
    fig.update_xaxes(title_text="Haplotype", row=1, col=2)
    fig.update_yaxes(title_text="Standard deviation of length", row=1, col=2)
    # make the subplot titles larger, which turn out to be annotations
    fig.update_annotations(font_size=32)
    fig.write_html(args.output)


def get_args():
    parser = ArgumentParser("Plot somatic differences (stdev) per haplotype and group")
    parser.add_argument("vcfs", nargs="+", help="VCF files to parse")
    parser.add_argument(
        "-o", "--output", help="output plot", default="somatic_variation.html"
    )
    parser.add_argument(
        "-m",
        "--minlen",
        help="minimal (consensus) length to plot",
        default=50,
        type=int,
    )
    parser.add_argument("--sample_info", help="TSV file with sample information")
    return parser.parse_args()


if __name__ == "__main__":
    main()
