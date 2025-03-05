import plotly.express as px
import pandas as pd
from argparse import ArgumentParser


def main():
    args = get_args()
    df = pd.read_csv(args.input, sep="\t").drop_duplicates(subset=["name"])

    df = df[(df["group"].isin(["aFTLD-U", "in-house non-aFTLD-U", "1000G"]))]

    haplotype_alias = {
        "major": "HapA",
        "minor": "HapB",
        "none": "None"
    }
    df["Haplotype"] = df["haplotype"].map(haplotype_alias)

    fig = px.strip(
        df.sort_values("Haplotype", ascending=False),
        x="CT_dimer_count",
        y="Haplotype",
        color="group",
        hover_data=["name"],
        title="CT dimer count",
        labels={"CT_dimer_count": "Number of CT dimers"},
        color_discrete_map={
            "aFTLD-U": "red",
            "in-house non-aFTLD-U": "black",
            "1000G": "teal",
        },
    )
    fig.add_vline(x=190, line_width=2, line_dash="dash", line_color="black")

    fig.update_traces(marker=dict(size=6, opacity=0.7))
    fig.update_layout(
        plot_bgcolor="white",
        font=dict(size=20),
        height=300,
        width=1000,
        margin=dict(l=0, r=0, t=50, b=0),
        legend=dict(
            title="",
            yanchor="bottom",
            y=0,
            xanchor="right",
            x=1,
            bordercolor="black",
            borderwidth=0.5,
            font=dict(size=14),
            orientation="h",
        ),
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
        title="",
        tickangle=30
    )

    with open(args.output, "w") as output:
        output.write(fig.to_html(include_plotlyjs="cdn"))


def get_args():
    parser = ArgumentParser("Correlate age with variables related to the repeat")
    parser.add_argument("-i", "--input", help="Path to the input table", required=True)
    parser.add_argument("-o", "--output", default="ct-dimer-strip.html", help="Path to the output html file")
    return parser.parse_args()


if __name__ == "__main__":
    main()
