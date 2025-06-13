def generate_manhattan_plot(df, genes):
    significance_threshold = 5e-8

    fig = px.scatter(
        df,
        x="newPOS",
        y="-log10P",
        color="CHROM_str",
        title="Manhattan plot",
        hover_data=["CHROM", "POS"],
    )
    fig.add_hline(y=-np.log10(significance_threshold), line_color="red", line_width=0.6)

    ays = {"LOC100288798": -85, "SLCO2B1": -55, "TMEM19": -30, "CNTNAP4": -45}
    axs = {
        "CNTNAP4": 2,
        "TMEM19": 0,
        "SLCO2B1": -20,
        "PPP2R3B": 10,
        "LOC100288798": -15,
        "KCNN3": 10,
    }
    xanchors = {
        "TMEM19": "left",
        "SLCO2B1": "right",
        "PPP2R3B": "center",
        "CNTNAP4": "left",
        "AUTS2": "right",
        "LOC100288798": "right",
        "KCNN3": "left",
        "GOLGA8B": "right",
    }
    for _, row in genes.iterrows():
        if axs.get(row["Gene"], -10) == 0:
            xshift = 0
        elif axs.get(row["Gene"], -10) > 0:
            xshift = 3
        else:
            xshift = -3
        fig.add_annotation(
            x=row["newPOS"],
            y=-np.log10(row["P"]),
            text=row["Gene"],
            ax=axs.get(row["Gene"], -10),
            ay=ays.get(row["Gene"], -25),
            xanchor=xanchors.get(row["Gene"], "center"),
            yanchor="bottom",
            arrowhead=0,
            font=dict(size=24),
            yshift=3,
            xshift=xshift,
        )

    return fig


def prepare_data_for_plot(df):
    df["-log10P"] = -1 * np.log10(df["P"])
    df["newPOS"] = df.apply(
        lambda x: x["POS"]
        + sum([chrom_to_length[str(i)] for i in range(1, int(x["CHROM"]))]),
        axis=1,
    )
    df["CHROM_str"] = df["CHROM"].astype(str)
    return df


def filter_data(df, filters):
    if filters.get("p_value"):
        df = df[df["P"] <= filters["p_value"]]
    if filters.get("chromosome"):
        df = df[df["CHROM"].isin(filters["chromosome"])]
    return df