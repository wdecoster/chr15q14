import pandas as pd
import plotly.express as px
import numpy as np
from argparse import ArgumentParser
import sys
import os

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
    "Y": 57227415,
}


def downsample_manhattan_data(df, threshold=0.05, downsample_factor=100):
    """
    Downsample Manhattan plot data to make it manageable while keeping significant hits
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with P-values, already sorted by P-value
    threshold : float
        P-value threshold above which to downsample (default: 0.05)
    downsample_factor : int
        Keep every Nth point above threshold (default: 100)
    
    Returns:
    --------
    pandas.DataFrame
        Downsampled DataFrame
    """
    # Split data into significant and non-significant
    significant = df[df['P'] <= threshold].copy()
    non_significant = df[df['P'] > threshold].copy()
    
    print(f"Original data: {len(df):,} variants")
    print(f"Significant (P ≤ {threshold}): {len(significant):,} variants")
    print(f"Non-significant (P > {threshold}): {len(non_significant):,} variants")
    
    # Downsample the non-significant data
    if len(non_significant) > 0:
        # Keep every Nth point from the non-significant data
        downsampled_non_sig = non_significant.iloc[::downsample_factor].copy()
        print(f"Downsampled non-significant: {len(downsampled_non_sig):,} variants (every {downsample_factor}th point)")
        
        # Combine significant and downsampled non-significant
        final_df = pd.concat([significant, downsampled_non_sig], ignore_index=True)
    else:
        final_df = significant.copy()
    
    # Re-sort by P-value and reset index
    final_df = final_df.sort_values(by="P").reset_index(drop=True)
    print(f"Final dataset for plotting: {len(final_df):,} variants")
    
    return final_df


def main():
    args = get_args()
    df = load_data(args.input)
    
    # Save summary file BEFORE downsampling (using original full dataset)
    if args.save_summary:
        print(f"Saving summary file with {len(df):,} variants...")
        df[["CHROM", "POS", "P", "REF", "ALT", "OR"]].sort_values(
            by=["CHROM", "POS"]
        ).to_csv(args.save_summary, sep="\t", index=False)
    
    # Downsample the data for plotting
    df_downsampled = downsample_manhattan_data(
        df,
        threshold=args.downsample_threshold,
        downsample_factor=args.downsample_factor
    )
    
    # Only load genes if gene annotation is requested
    if not args.no_genes:
        if (
            os.path.basename(args.input)
            != "results.MAF0.01.VQSRpass.ctrlBatch.05.FUSBatch.05.ctrlHWE1E-8.ctrlMISS.05.caseMISS.05.forPublication.txt.gz"
        ):
            sys.stderr.write("Warning: Genes are hardcoded for the default input file\n")
        genes = load_genes(df_downsampled)  # Use downsampled data for gene loading
    else:
        genes = None
    
    # Create plot with downsampled data
    fig = plot(df_downsampled, genes)
    
    with open(args.output, "w") as f:
        f.write(fig.to_html(include_plotlyjs="cdn"))


def load_data(input_file):
    df = pd.read_csv(
        input_file,
        sep="\t",
        compression="gzip" if input_file.endswith('.gz') else None,
        usecols=["CHROM", "POS", "P", "REF", "ALT", "OR"],
    )
    df["-log10P"] = -1 * np.log10(df["P"])
    df["newPOS"] = df.apply(
        lambda x: x["POS"]
        + sum([chrom_to_length[str(i)] for i in range(1, int(x["CHROM"]))]),
        axis=1,
    )
    df["CHROM_str"] = df["CHROM"].astype(str)
    
    # Sort by P-value for downsampling
    df = df.sort_values(by="P").reset_index(drop=True)
    return df


def load_genes(df):
    genes = pd.DataFrame({"CHROM": ["15"], "POS": [34612468], "Gene": ["GOLGA8B"]})
    genes["newPOS"] = genes.apply(
        lambda x: x["POS"]
        + sum([chrom_to_length[str(i)] for i in range(1, int(x["CHROM"]))]),
        axis=1,
    )

    genes = (
        genes.set_index("newPOS")
        .join(df[["newPOS", "P"]].set_index("newPOS"), how="left", rsuffix="_r")
        .reset_index()
    )
    return genes


def plot(df, genes):
    significance_threshold = 5e-8

    fig = px.scatter(
        df,
        x="newPOS",
        y="-log10P",
        color="CHROM_str",
        hover_data=["CHROM", "POS"],
    )
    fig.add_hline(y=-np.log10(significance_threshold), line_color="red", line_width=0.6)

    # Only add gene annotations if genes are provided
    if genes is not None and not genes.empty:
        axs = {}
        ays = {}
        xanchors = {"GOLGA8B": "right"}
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

    # add the chromosomes as labels on the x-axis at the middle of the chromosome
    # however, the smaller chromosomes end up overlapping, so I will remove chr18, chr20 and chr22
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

    fig.update_traces(marker=dict(size=6))
    fig.update_layout(
        plot_bgcolor="white",
        paper_bgcolor="white",
        height=800,
        width=1000,
        showlegend=False,
        margin=dict(l=0, r=0, t=0, b=0),
        xaxis_title="",
        xaxis=dict(
            tickmode="array",
            tickvals=list(chroms.values()),
            ticktext=list(chroms.keys()),
            ticklen=0,
            tickangle=-90,
            tickfont=dict(size=20),
        ),
        yaxis=dict(
            tickfont=dict(size=20),
            title=dict(text="-log10(P)", font=dict(size=22)),
        ),
    )
    return fig


def get_args():
    parser = ArgumentParser()
    parser.add_argument("-i",
        "--input",
        type=str,
        help="Path to the input file",
        default="~/local/results.MAF0.01.VQSRpass.ctrlBatch.05.FUSBatch.05.ctrlHWE1E-8.ctrlMISS.05.caseMISS.05.forPublication.txt.gz",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Path to the output plot file",
        default="manhattan.html",
    )
    parser.add_argument("--save_summary", help="Save the summary file", default=False)
    parser.add_argument(
        "--no-genes",
        action="store_true",
        help="Don't annotate genes on the Manhattan plot"
    )
    parser.add_argument(
        "--downsample-threshold",
        type=float,
        default=0.05,
        help="P-value threshold for downsampling (default: 0.05)"
    )
    parser.add_argument(
        "--downsample-factor",
        type=int,
        default=100,
        help="Factor for downsampling (keep every Nth point, default: 100)"
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
