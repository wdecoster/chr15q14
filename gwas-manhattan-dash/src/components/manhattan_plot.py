from dash import dcc, html
import plotly.express as px
import pandas as pd
import numpy as np

def generate_manhattan_plot(df):
    """Generate a Manhattan plot from the filtered data"""
    significance_threshold = 5e-8
    
    # Check if required columns are present, add them to hover_data as appropriate
    hover_columns = []
    if "original_chrom" in df.columns and "original_pos" in df.columns:
        hover_columns = ["original_chrom", "original_pos", "P", "-log10P"]
    else:
        hover_columns = ["CHROM_str", "P", "-log10P"]
    
    fig = px.scatter(
        df,
        x="newPOS",
        y="-log10P",
        color="CHROM_str",
        title="Manhattan Plot",
        hover_data=hover_columns,
        labels={
            "original_chrom": "Chromosome",
            "original_pos": "Position",
            "P": "P-value",
            "-log10P": "-log10(P)"
        }
    )
    
    fig.add_hline(y=-np.log10(significance_threshold), line_color="red", line_width=0.6)

    # Add chromosome labels
    chrom_to_length = {
        "1": 248956422,
        "2": 242193529,
        "3": 198295559,
        "4": 190214555,
        "5": 181538259,
        "6": 170805979,
        "7": 159345973,
        "8": 145138636,
        "9": 138394717,
        "10": 133797422,
        "11": 135086622,
        "12": 133275309,
        "13": 114364328,
        "14": 107043718,
        "15": 101991189,
        "16": 90338345,
        "17": 83257441,
        "18": 80373285,
        "19": 58617616,
        "20": 64444167,
        "21": 46709983,
        "22": 50818468,
        "23": 156040895,
    }

    chromosome_middle = [
        sum([chrom_to_length[str(i)] for i in range(1, int(chrom))])
        + chrom_to_length[str(chrom)] / 2
        for chrom in range(1, 24)
    ]

    chroms = {
        ch: m
        for ch, m in zip(
            [f"chr{chrom}" for chrom in range(1, 23)] + ["chrX"], chromosome_middle
        )
        if ch not in ["chr18", "chr20", "chr22"]
    }

    fig.update_traces(marker=dict(size=4))
    fig.update_layout(
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=800,
        showlegend=False,
        margin=dict(l=0, r=0, t=30, b=0),
        xaxis_title="",
        xaxis=dict(
            tickmode="array",
            tickvals=list(chroms.values()),
            ticktext=list(chroms.keys()),
            ticklen=0,
            tickangle=-90,
            tickfont=dict(size=14),
        ),
        yaxis=dict(
            tickfont=dict(size=14),
            title=dict(text="-log10(P)", font=dict(size=16)),
        ),
    )
    
    return fig

def manhattan_plot_layout(df):
    return html.Div([
        html.Div(id='manhattan-plot-container', children=generate_manhattan_plot(df)),
    ])
