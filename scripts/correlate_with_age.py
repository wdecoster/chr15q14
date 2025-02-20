import plotly.express as px
import pandas as pd
from argparse import ArgumentParser
import statsmodels.api as sm


class OLS:
    def __init__(self, df, variable1, variable2):
        X = df[variable1]
        y = df[variable2]
        X = sm.add_constant(X)
        model = sm.OLS(y, X).fit()
        model.summary()

        # get the pvalue
        self.pval = model.pvalues[variable1]
        self.rsquared = model.rsquared

    def __str__(self):
        return f"R<sup>2</sup> = {self.rsquared:.2f}<br>p-value = {self.pval:.2f}"


def main():
    args = get_args()
    df = pd.read_csv(args.input, sep="\t")
    sampleinfo = pd.read_excel(
        args.sampleinfo, usecols=["Gentli_ID", "AGEATDEATH", "AGEATONSET"]
    ).rename(columns={"Gentli_ID": "name"})
    df = (
        df.merge(sampleinfo, on="name")
        .dropna(subset=["AGEATDEATH"])
        .drop_duplicates(subset=["name"])
    )
    if args.pat_only:
        df = df[df["group"] == "aFTLD-U"]
    # only keep the aFTLD-U and in-house non-aFTLD-U samples, which have an age at death and a haplotype, and removing problematic samples
    # problematic because of (respectively): no repeat, weird repeat, cerebellum sample, cerebellum sample
    df = df[
        (df["group"].isin(["aFTLD-U", "in-house non-aFTLD-U"]))
        & (df["AGEATDEATH"] > 0)
        & (df["haplotype"].isin(["major", "minor"]))
        & (~df["name"].isin(["rr_PIDN4463", "rr_UCL2783", "rr_R556", "rr_RZ257"]))
    ]
    df["repeat_length_CT"] = df["length"] * df["%CT"]

    annotation_location = {
            "aFTLD-U": (0.99, 0.70),
            "in-house non-aFTLD-U": (0.80, 0.70),
        }

    make_plot(df, args, annotation_location, variable="AGEATDEATH")
    make_plot(df, args, annotation_location, variable="AGEATONSET")

def make_plot(df, args, annotation_location, variable):
    fig = px.scatter(
        df,
        x="CT_dimer_count",
        y=variable,
        color="group",
        trendline="ols",
        hover_data=["name"],
        labels={
            "length": "Repeat length",
            "AGEATDEATH": "Age at death",
            "AGEATONSET": "Age at onset",
            "repeat_length_CT": "Repeat length * %CT",
            "CT_dimer_count": "CT dimer count",
        },
        color_discrete_map={"aFTLD-U": "red", "in-house non-aFTLD-U": "black"},
        title=f"Carriers of a chr15q14 expansion",
    )
    fig.update_traces(marker=dict(size=6, opacity=0.5))
    fig.update_layout(
        plot_bgcolor="white",
        font=dict(size=20),
        legend=dict(
            title="Group",
            yanchor="top",
            y=1,
            xanchor="right",
            x=1,
            bordercolor="black",
            borderwidth=1,
        ),
        width=800,
        height=800,
        margin=dict(l=0, r=0, t=50, b=0),
    )
    if args.pat_only:
        fig.update_layout(showlegend=False)
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

    # add OLS annotation for patients and controls
    # get the end of the ols line

    fig.add_annotation(
        x=annotation_location["aFTLD-U"][0],
        y=annotation_location["aFTLD-U"][1],
        xref="x domain",
        yref="y domain",
        text=OLS(
            df[df["group"].isin(["aFTLD-U"])],
            "CT_dimer_count",
            variable,
        ).__str__(),
        align="left",
        showarrow=False,
        font=dict(color="red", size=18),
    )
    if not args.pat_only:
        fig.add_annotation(
            x=annotation_location["in-house non-aFTLD-U"][0],
            y=annotation_location["in-house non-aFTLD-U"][1],
            xref="x domain",
            yref="y domain",
            text=OLS(
                df[df["group"].isin(["in-house non-aFTLD-U"])],
                "CT_dimer_count",
            ).__str__(),
            align="left",
            showarrow=False,
            font=dict(color="black", size=18),
        )
    print(
        fig.to_html(
            include_plotlyjs="cdn",
            config={
                "toImageButtonOptions": {
                    "format": "png",  # one of png, svg, jpeg, webp
                    "filename": "",
                    "height": 800,
                    "width": 800,
                    "scale": 5,  # Multiply title/legend/axis/canvas sizes by this factor
                }
            },
        )
    )


def get_args():
    parser = ArgumentParser("Correlate age with variables related to the repeat")
    parser.add_argument("input", help="Path to the input table")
    parser.add_argument("--sampleinfo", help="Path to Individuals.xlsx")
    parser.add_argument("--pat_only", help="Only include patients", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    main()
