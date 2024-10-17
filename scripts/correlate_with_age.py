import plotly.express as px
import pandas as pd
from argparse import ArgumentParser
import statsmodels.api as sm


class OLS:
    def __init__(self, df, variable):
        X = df[variable]
        y = df["AGEATDEATH"]
        X = sm.add_constant(X)
        model = sm.OLS(y, X).fit()
        model.summary()

        # get the pvalue
        self.pval = model.pvalues[variable]
        self.rsquared = model.rsquared

    def __str__(self):
        return f"R<sup>2</sup> = {self.rsquared:.2f}<br>p-value = {self.pval:.2f}"


def main():
    args = get_args()
    df = pd.read_csv(args.input, sep="\t")
    sampleinfo = pd.read_excel(
        args.sampleinfo, usecols=["Gentli_ID", "AGEATDEATH"]
    ).rename(columns={"Gentli_ID": "name"})
    df = (
        df.merge(sampleinfo, on="name")
        .dropna(subset=["AGEATDEATH"])
        .drop_duplicates(subset=["name"])
    )
    if args.pat_only:
        df = df[df["group"] == "aFTLD-U"]
    df = df[
        (df["group"].isin(["aFTLD-U", "in-house control"]))
        & (df["AGEATDEATH"] > 0)
        & (df["haplotype"].isin(["major", "minor"]))
        & (~df["name"].isin(["rr_PIDN4463", "rr_UCL2783"]))
    ]
    df["repeat_length_CT"] = df["length"] * df["%CT"]
    df["CT_dimer_count"] = df["CT_dimer_count"] * 2

    haplotype_alias = {
        "major": "haplotype A",
        "minor": "haplotype B",
    }

    annotation_location = {
        "major": {
            "aFTLD-U": (0.99, 0.70),
            "in-house control": (0.80, 0.70),
        },
        "minor": {
            "aFTLD-U": (0.99, 0.45),
            "in-house control": (0.30, 0.65),
        },
    }

    for haplotype in ["major", "minor"]:
        for variable in ["length", "repeat_length_CT", "CT_dimer_count"]:
            fig = px.scatter(
                df[df["haplotype"] == haplotype],
                x=variable,
                y="AGEATDEATH",
                color="group",
                trendline="ols",
                hover_data=["name"],
                labels={
                    "length": "Repeat length",
                    "AGEATDEATH": "Age at death",
                    "repeat_length_CT": "Repeat length * %CT",
                    "CT_dimer_count": "CT dimer count",
                },
                color_discrete_map={"aFTLD-U": "red", "in-house control": "black"},
                title=f"Carriers of {haplotype_alias[haplotype]}",
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

            # add OLS annotation for patients and controls
            # get the end of the ols line

            fig.add_annotation(
                x=annotation_location[haplotype]["aFTLD-U"][0],
                y=annotation_location[haplotype]["aFTLD-U"][1],
                xref="x domain",
                yref="y domain",
                text=OLS(
                    df[
                        (df["haplotype"] == haplotype) & (df["group"].isin(["aFTLD-U"]))
                    ],
                    variable,
                ).__str__(),
                align="left",
                showarrow=False,
                font=dict(color="red", size=18),
            )
            if not args.pat_only:
                fig.add_annotation(
                    x=annotation_location[haplotype]["in-house control"][0],
                    y=annotation_location[haplotype]["in-house control"][1],
                    xref="x domain",
                    yref="y domain",
                    text=OLS(
                        df[
                            (df["haplotype"] == haplotype)
                            & (df["group"].isin(["in-house control"]))
                        ],
                        variable,
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
