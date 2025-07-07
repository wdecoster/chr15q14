# This script generates a QQ plot for the given dataset
# The dataset is already filtered, and has column for CHROM, POS, P, REF, ALT and OR.


import numpy as np


def get_args():
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Generate a QQ plot from the input data.")
    parser.add_argument(
        "-i", "--input", type=str, required=True,
        help="Input file containing the data for the QQ plot."
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True,
        help="Output file to save the generated QQ plot.",
        default="qq_plot.html"
    )
    parser.add_argument(
        "--downsample-threshold", type=float, default=0.05,
        help="P-value threshold above which to downsample data (default: 0.05)"
    )
    parser.add_argument(
        "--downsample-factor", type=int, default=100,
        help="Keep every Nth point above the threshold (default: 100)"
    )
    return parser.parse_args()

def downsample_qq_data(df, threshold=0.05, downsample_factor=100):
    """
    Downsample QQ plot data to make it manageable while keeping significant hits
    
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
    import pandas as pd
    
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
    import pandas as pd
    import numpy as np
    import plotly.graph_objects as go
    from kaleido.scopes.plotly import PlotlyScope
    
    # Load and sort data
    print("Loading data...")
    df = pd.read_csv(args.input, sep="\t", compression='gzip' if args.input.endswith('.gz') else None)
    df = df.sort_values(by="P").reset_index(drop=True)
    df["expected"] = -np.log10(np.arange(1, len(df) + 1) / len(df))  # Expected values for QQ plot
    df["observed"] = -np.log10(df["P"])  # Observed values for QQ plot
    
    # Downsample the data
    df_downsampled = downsample_qq_data(
        df, 
        threshold=args.downsample_threshold, 
        downsample_factor=args.downsample_factor
    )
    
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=df_downsampled['expected'],
        y=df_downsampled['observed'],
        mode='markers',
        marker=dict(size=3, color='green'),
        name=f'P > {args.downsample_threshold} (downsampled)',
        hovertemplate='<b>Expected:</b> %{x:.3f}<br><b>Observed:</b> %{y:.3f}<extra></extra>'
    ))
    

    
    # Add diagonal reference line
    max_val = max(df_downsampled['expected'].max(), df_downsampled['observed'].max())
    fig.add_trace(go.Scatter(
        x=[0, max_val],
        y=[0, max_val],
        mode='lines',
        line=dict(color='grey', dash='dash', width=2),
        name='y = x',
        hoverinfo='skip'
    ))
    
    # Calculate lambda (genomic inflation factor)
    chi2_stats = -2 * np.log(df["P"])
    observed_median = np.median(chi2_stats)
    expected_median = -2 * np.log(0.5)  # This equals ~1.386
    lambda_gc = observed_median / expected_median
    
    fig.update_layout(
        title=f'QQ Plot<br><sub>λ<sub>GC</sub> = {lambda_gc:.3f}, {len(df):,} total variants, {len(df_downsampled):,} plotted</sub>',
        xaxis_title='Expected -log₁₀(p)',
        yaxis_title='Observed -log₁₀(p)',
        showlegend=False,
        height=800,
        width=800,
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis=dict(
            showgrid=True,
            gridcolor='lightgray',
            gridwidth=1,
            zeroline=True,
            zerolinecolor='gray',
            zerolinewidth=1
        ),
        yaxis=dict(
            showgrid=True,
            gridcolor='lightgray',
            gridwidth=1,
            zeroline=True,
            zerolinecolor='gray',
            zerolinewidth=1
        ),
    )
    
    # Save the plot
    print(f"Saving QQ plot to {args.output}...")
    fig.write_html(args.output)
    
    # Save PNG version
    png_output = args.output.replace(".html", ".png")
    print(f"Saving PNG version to {png_output}...")
    scope = PlotlyScope()
    with open(png_output, "wb") as f:
        f.write(scope.transform(fig, format="png"))
    
    print("QQ plot generation completed!")
    print(f"Original dataset: {len(df):,} variants")
    print(f"Plotted dataset: {len(df_downsampled):,} variants")
    print(f"Genomic inflation factor (λ_GC): {lambda_gc:.3f}")


if __name__ == "__main__":
    main()
