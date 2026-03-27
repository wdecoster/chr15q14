import logging
import sys
from itertools import cycle

import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots


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


def scatter_plot(df, motif, full=False):
    fig_scatter = go.Figure()
    colors = ["blue", "red", "green", "purple", "orange"]
    if len(df["sample"].unique()) > len(colors):
        sys.stderr.write(
            "Warning: more samples than colors in scatter plot, some samples will have the same color.\n"
        )
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
        title="Read length" if full else "non-ref length",
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


def plot_genotypes(genotypes):
    """
    Make a scatter plot showing for each sample the CT dimer count vs the CCCTCT hexamer count
    """
    fig_scatter = go.Figure(
        go.Scatter(
            x=[genotype.ct_dimer_count for genotype in genotypes],
            y=[genotype.hexamer_count for genotype in genotypes],
            mode="markers",
            text=[genotype.individual for genotype in genotypes],
            marker=dict(size=4, color="black"),
            hovertext=[genotype.individual for genotype in genotypes],
        )
    )
    fig_scatter.update_layout(
        plot_bgcolor="white",
        margin=dict(l=0, r=0, t=50, b=0),
        title="Genotypes from spanning PCR",
        font=dict(size=16),
        height=400,
        width=800,
    )
    fig_scatter.update_xaxes(
        showline=True,
        linewidth=1,
        linecolor="black",
        mirror=True,
        title="CT dimers",
        rangemode="tozero"
    )
    fig_scatter.update_yaxes(
        showline=True,
        linewidth=1,
        linecolor="black",
        mirror=True,
        title="CCCTCT hexamers",
    )
    fig_scatter.add_vline(x=190, line_width=2, line_dash="dash", line_color="black")
    return fig_scatter


def plot_violins(df, genotypes, min_ct_count=100):
    """
    Make a violin plot, with on the y-axis the samples and on the x-axis the CT count of the reads
    Samples on the y-axis should be sorted by the consensus ct_dimer_count
    """
    genotypes = [genotype for genotype in genotypes if genotype.ct_dimer_count > min_ct_count]
    sorted_genotypes = sorted(genotypes, key=lambda g: g.ct_dimer_count)
    sample_order = [genotype.individual for genotype in sorted_genotypes]
    df = df[df["sample"].isin(sample_order)]

    fig_violin = px.violin(df, x="sample", y="CT count", category_orders={"sample": sample_order},
                           points=False, title=f"Violin plot of CT counts by sample (minimum {min_ct_count} CTs in consensus sequence)")
    fig_violin.update_traces(spanmode="hard")
    # on the same plot, add the consensus length for each sample with a red dot
    fig_violin.add_trace(
        go.Scatter(
            x=sample_order,
            y=[genotype.ct_dimer_count for genotype in sorted_genotypes],
            mode="markers",
            marker=dict(color="red", size=6),
            name="Consensus Lengths",
        )
    )
    fig_violin.update_layout(plot_bgcolor="white", margin=dict(l=0, r=0, t=50, b=0))
    fig_violin.update_xaxes(showline=True, linewidth=1, linecolor="black", mirror=True, rangemode="nonnegative")
    fig_violin.update_yaxes(showline=True, linewidth=1, linecolor="black", mirror=True)
    fig_violin.add_hline(y=190, line_width=2, line_dash="dash", line_color="black")
    return fig_violin


def plot_truth_correlation(labels, cram_lengths, truth_lengths, cram_ct, truth_ct, both=False):
    if both:
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=("Expansion length (bp)", "CT dimer count"),
            horizontal_spacing=0.12,
        )
        fig.update_annotations(font=dict(size=18)) # increase subplot title font size

        # Scatter 1: expansion length
        fig.add_trace(
            go.Scatter(
                x=truth_lengths,
                y=cram_lengths,
                mode="markers",
                marker=dict(size=9, color="black"),
                text=labels,
                hovertext=labels,
                showlegend=False,
            ),
            row=1, col=1,
        )
        # identity line for length
        all_lengths = truth_lengths + cram_lengths
        len_min, len_max = min(all_lengths), max(all_lengths)
        fig.add_trace(
            go.Scatter(
                x=[len_min, len_max], y=[len_min, len_max],
                mode="lines", line=dict(dash="dash", color="grey", width=1),
                showlegend=False,
            ),
            row=1, col=1,
        )

        # Scatter 2: CT dimer count
        fig.add_trace(
            go.Scatter(
                x=truth_ct,
                y=cram_ct,
                mode="markers",
                marker=dict(size=9, color="black"),
                text=labels,
                hovertext=labels,
                showlegend=False,
            ),
            row=1, col=2,
        )
        # identity line for CT
        all_ct = truth_ct + cram_ct
        ct_min, ct_max = min(all_ct), max(all_ct)
        fig.add_trace(
            go.Scatter(
                x=[ct_min, ct_max], y=[ct_min, ct_max],
                mode="lines", line=dict(dash="dash", color="grey", width=1),
                showlegend=False,
            ),
            row=1, col=2,
        )

        fig.update_layout(
            plot_bgcolor="white",
            margin=dict(l=60, r=20, t=60, b=60),
            font=dict(size=18),
            height=500,
            width=1100,
        )
        for col in [1, 2]:
            fig.update_xaxes(
                showline=True, linewidth=1, linecolor="black", mirror=True,
                rangemode="tozero", row=1, col=col,
            )
            fig.update_yaxes(
                showline=True, linewidth=1, linecolor="black", mirror=True,
                rangemode="tozero", row=1, col=col,
            )
        fig.update_xaxes(title_text="ONT genome sequencing", row=1, col=1)
        fig.update_yaxes(title_text="Amplicon ", row=1, col=1)
        fig.update_xaxes(title_text="ONT genome sequencing", row=1, col=2)
        fig.update_yaxes(title_text="Amplicon sequencing", row=1, col=2)

    else:
        # this only does the CT dimer plot, without subplots, but otherwise the same style as above
        fig = go.Figure(
            go.Scatter(
                x=truth_ct,
                y=cram_ct,
                mode="markers",
                marker=dict(size=9, color="black"),
                text=labels,
                hovertext=labels,
                showlegend=False,
            )
        )
        # identity line for CT
        all_ct = truth_ct + cram_ct
        ct_min, ct_max = min(all_ct), max(all_ct)
        fig.add_trace(
            go.Scatter(
                x=[ct_min, ct_max], y=[ct_min, ct_max],
                mode="lines", line=dict(dash="dash", color="grey", width=1),
                showlegend=False,
            )
        )
        fig.update_layout(
            plot_bgcolor="white",
            margin=dict(l=60, r=20, t=60, b=60),
            font=dict(size=18),
            height=400,
            width=400,
            title=dict(text="CT dimer count", x=0.5, xanchor="center"),
        )
        fig.update_xaxes(
            showline=True, linewidth=1, linecolor="black", mirror=True,
            rangemode="tozero",
        )
        fig.update_yaxes(
            showline=True, linewidth=1, linecolor="black", mirror=True,
            rangemode="tozero",
        )
        fig.update_xaxes(title_text="ONT genome sequencing")
        fig.update_yaxes(title_text="Amplicon sequencing")
    with open("repeat_correlation.html", "w") as f:
        f.write(fig.to_html())
    logging.info("Wrote repeat_correlation.html")
