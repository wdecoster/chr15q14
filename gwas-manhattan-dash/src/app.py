from dash import Dash, dcc, html, Input, Output, State, ctx, dash_table
import pandas as pd
import dash_bootstrap_components as dbc
from components.filters import create_filter_layout
from components.manhattan_plot import generate_manhattan_plot
from utils.data_loader import load_data, filter_data, get_vqsr_values
import time
import gzip
from datetime import datetime
import argparse  # Add this import

# Add command line argument parsing
def parse_arguments():
    parser = argparse.ArgumentParser(description='GWAS Manhattan Plot Dashboard')
    parser.add_argument(
        '--data-file', 
        '-f',
        type=str,
        default="/home/wdecoster/local/CC.regenie.csv.gz",
        help='Path to the GWAS data file (default: /home/wdecoster/local/CC.regenie.csv.gz)'
    )
    parser.add_argument(
        '--pval-threshold',
        '-p',
        type=float,
        default=1e-2,
        help='P-value threshold for pre-filtering (default: 1e-2)'
    )
    parser.add_argument(
        '--port',
        type=int,
        default=8050,
        help='Port to run the app on (default: 8050)'
    )
    parser.add_argument(
        '--host',
        type=str,
        default='127.0.0.1',
        help='Host to run the app on (default: 127.0.0.1)'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Run in debug mode'
    )
    return parser.parse_args()

# Parse command line arguments
args = parse_arguments()

# Initialize the app with Bootstrap for styling
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Load the pre-filtered data for initial fast rendering
print(f"Loading pre-filtered data from: {args.data_file}")
prefiltered_data = load_data(args.data_file, pval_threshold=args.pval_threshold)
full_data = None  # Will be loaded on demand

# Get unique VQSR values for filter
vqsr_values = get_vqsr_values(prefiltered_data)

# Define the layout
app.layout = html.Div([
    html.H1("GWAS Manhattan Plot", className="mt-3 mb-4 text-center"),
    
    # Main content container - 2/3 plot, 1/3 filters
    html.Div([
        # Manhattan plot container (left side, 2/3 width)
        html.Div([
            # Variant count display
            html.Div(id='variant-count', className='variant-count-label'),
            
            # Button container for action button
            html.Div([
                # Single button to filter full dataset and download
                html.Button(
                    'Filter ALL variants and download (CHROM, POS, P, REF, ALT, OR)', 
                    id='filter-and-download-button',
                    className='filter-download-button',
                    n_clicks=0
                ),
                
                # Download component
                dcc.Download(id="download-dataframe-tsv")
            ], style={'marginBottom': '10px'}),
            
            # Add an informational note
            html.Div([
                html.P(
                    "Note: Manhattan plot shows pre-filtered data for performance. "
                    "Use the button above to apply current filters to ALL variants and download.",
                    style={
                        'fontSize': '12px', 
                        'color': '#666', 
                        'fontStyle': 'italic',
                        'marginBottom': '10px'
                    }
                )
            ]),
            
            dcc.Loading(
                id="loading",
                type="circle",
                children=dcc.Graph(id='manhattan-plot', className="manhattan-plot")
            ),
            
            dcc.Loading(
                id="table-loading",
                type="circle",
                children=html.Div(id='significant-variants-table')
            ),
            
            # Store for download status
            dcc.Store(id='download-status', data={})
        ], className="plot-container"),
        
        # Filters container (right side, 1/3 width)
        html.Div([
            create_filter_layout(vqsr_values)
        ], className="sidebar-container")
    ], className="main-container")
], className="app-container")

# Simplified callback for the plot (only uses pre-filtered data)
@app.callback(
    [Output('manhattan-plot', 'figure'),
     Output('variant-count', 'children'),
     Output('significant-variants-table', 'children')],
    [Input('apply-filters-button', 'n_clicks')],
    [State('vqsr-filter', 'value'),
     State('hwe-ctrl-filter', 'value'),
     State('batch-fepval-fus-filter', 'value'),
     State('ctrl-f-miss-filter', 'value'),
     State('case-f-miss-filter', 'value'),
     State('allele-freq-filter', 'value'),
     State('batch-pval-ctrl-filter', 'value')],
    prevent_initial_call=False
)
def update_plot(apply_clicks, vqsr_values, hwe_ctrl, batch_fepval_fus, 
                ctrl_f_miss, case_f_miss, allele_freq, batch_pval_ctrl):
    
    # Apply filters to the pre-filtered data (for plot only)
    try:
        filtered_plot_data = filter_data(
            prefiltered_data,
            vqsr_values=vqsr_values,
            hwe_ctrl_cutoff=float(hwe_ctrl) if hwe_ctrl else None,
            batch_fepval_fus_cutoff=float(batch_fepval_fus) if batch_fepval_fus else None,
            ctrl_f_miss_cutoff=float(ctrl_f_miss) if ctrl_f_miss else None,
            case_f_miss_cutoff=float(case_f_miss) if case_f_miss else None,
            allele_freq_cutoff=float(allele_freq) if allele_freq else None,
            batch_pval_ctrl_cutoff=float(batch_pval_ctrl) if batch_pval_ctrl else None,
        )
        
        # Generate the variant count label
        plot_count = len(filtered_plot_data)
        variant_count_text = f"Displaying {plot_count:,} variants (from pre-filtered data)"
        
        # Get significant peak variants from the plot data
        peak_variants = get_peak_variants(filtered_plot_data)
        
        # Create the table for significant variants
        table_title = "Significant Peaks (one variant per peak)"
        
        if len(peak_variants) > 0:
            # Format p-values properly
            display_data = peak_variants.copy()
            if 'P-value' in display_data.columns:
                display_data['P-value'] = display_data['P-value'].apply(
                    lambda x: f"{x:.2e}" if x < 1e-4 else f"{x:.6f}"
                )
            
            # Round other numeric columns
            numeric_columns = display_data.select_dtypes(include=['float64', 'int64']).columns
            for col in numeric_columns:
                if col != 'P-value':
                    if col in ['-log10(P)', 'OR']:
                        display_data[col] = display_data[col].round(4)
                    else:
                        display_data[col] = display_data[col].round(6)
            
            peak_table = html.Div([
                html.H3(table_title, className="mt-4"),
                dash_table.DataTable(
                    data=display_data.to_dict('records'),
                    columns=[{"name": i, "id": i} for i in display_data.columns],
                    style_table={'overflowX': 'auto'},
                    style_cell={
                        'textAlign': 'center',
                        'padding': '10px',
                        'minWidth': '100px',
                        'whiteSpace': 'normal',
                        'height': 'auto',
                    },
                    style_header={
                        'backgroundColor': 'rgb(230, 230, 230)',
                        'fontWeight': 'bold',
                        'textAlign': 'center'
                    },
                    style_data_conditional=[
                        {
                            'if': {'column_id': 'P-value'},
                            'fontWeight': 'bold',
                            'fontFamily': 'monospace',
                        }
                    ],
                    sort_action='native',
                    filter_action='native',
                    page_size=10,
                    cell_selectable=True,
                    row_selectable=False,
                    selected_cells=[],
                    export_format='xlsx',
                    export_headers='display',
                    css=[{
                        'selector': '.dash-cell div.dash-cell-value', 
                        'rule': 'display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;'
                    }]
                )
            ])
        else:
            peak_table = html.Div([
                html.H3(table_title, className="mt-4"),
                html.Div("No significant variants found", style={'color': 'gray', 'textAlign': 'center'})
            ])
        
        # Generate Manhattan plot
        manhattan_figure = generate_manhattan_plot(filtered_plot_data)
        
        return manhattan_figure, variant_count_text, peak_table
    
    except Exception as e:
        import traceback
        print(f"Error in update_plot: {e}")
        print(traceback.format_exc())
        empty_table = html.Div([
            html.H3("Significant Peaks (one variant per peak)", className="mt-4"),
            html.Div("Error generating table: " + str(e), style={'color': 'red'})
        ])
        return {}, f"Error: {str(e)}", empty_table

# New callback for filtering full dataset and direct download
@app.callback(
    Output("download-dataframe-tsv", "data"),
    Input("filter-and-download-button", "n_clicks"),
    [State('vqsr-filter', 'value'),
     State('hwe-ctrl-filter', 'value'),
     State('batch-fepval-fus-filter', 'value'),
     State('ctrl-f-miss-filter', 'value'),
     State('case-f-miss-filter', 'value'),
     State('allele-freq-filter', 'value'),
     State('batch-pval-ctrl-filter', 'value')],
    prevent_initial_call=True,
)
def filter_and_download_full_data(n_clicks, vqsr_values, hwe_ctrl, batch_fepval_fus, 
                                  ctrl_f_miss, case_f_miss, allele_freq, batch_pval_ctrl):
    if n_clicks is None or n_clicks == 0:
        return {}

    try:
        print("Loading full dataset with all columns for filtering and download...")
        start_time = time.time()

        # Load the full dataset with ALL required columns (including REF, ALT)
        full_data = load_full_data_with_all_columns(args.data_file)
        load_time = time.time() - start_time
        print(f"Full data loaded in {load_time:.2f} seconds ({len(full_data):,} variants)")

        # Apply filters to the full dataset
        print("Applying filters to full dataset...")
        filter_start = time.time()
        filtered_full_data = filter_data(
            full_data,
            vqsr_values=vqsr_values,
            hwe_ctrl_cutoff=float(hwe_ctrl) if hwe_ctrl else None,
            batch_fepval_fus_cutoff=float(batch_fepval_fus) if batch_fepval_fus else None,
            ctrl_f_miss_cutoff=float(ctrl_f_miss) if ctrl_f_miss else None,
            case_f_miss_cutoff=float(case_f_miss) if case_f_miss else None,
            allele_freq_cutoff=float(allele_freq) if allele_freq else None,
            batch_pval_ctrl_cutoff=float(batch_pval_ctrl) if batch_pval_ctrl else None,
        )
        filter_time = time.time() - filter_start
        print(f"Filtering completed in {filter_time:.2f} seconds ({len(filtered_full_data):,} variants remaining)")

        # Check which columns exist and create the final output dataframe
        final_output = filtered_full_data.copy()

        # Sort by chromosome and position
        final_output = final_output[
            ["CHROM", "POS", "P", "REF", "ALT", "OR"]
        ].sort_values(by=["CHROM", "POS"])

        print(f"Final output columns: {list(final_output.columns)}")

        # Create filename with timestamp
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"gwas_filtered_full_dataset_{timestamp}.tsv.gz"

        print(f"Preparing download of {len(final_output):,} variants with columns: {list(final_output.columns)}")

        # Convert DataFrame to compressed TSV using pandas directly
        import io
        buffer = io.BytesIO()

        # Write to gzipped buffer
        with gzip.open(buffer, 'wt', encoding='utf-8') as f:
            final_output.to_csv(f, sep='\t', index=False)

        # Get the compressed content
        buffer.seek(0)
        compressed_content = buffer.getvalue()

        print(f"Download ready: {filename} ({len(compressed_content)} bytes)")

        return dcc.send_bytes(
            compressed_content,
            filename=filename,
            type="application/gzip"
        )

    except Exception as e:
        print(f"Error in filter_and_download_full_data: {e}")
        import traceback
        traceback.print_exc()
        return {}

# Update the CSS to style the new button
app.index_string = '''
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
        <style>
            .app-container {
                margin: 0 auto;
                max-width: 1400px;
                padding: 20px;
            }
            
            .main-container {
                display: flex;
                flex-wrap: wrap;
            }
            
            .plot-container {
                width: 66%;
                padding-right: 20px;
            }
            
            .sidebar-container {
                width: 34%;
                padding: 15px;
                background-color: #f8f9fa;
                border-radius: 5px;
            }
            
            .filter-container {
                display: flex;
                flex-direction: column;
                gap: 15px;
            }
            
            .filter-group {
                margin-bottom: 10px;
            }
            
            .filter-label {
                display: block;
                margin-bottom: 5px;
                font-weight: bold;
            }
            
            .input-field {
                width: 100%;
                padding: 8px;
                border: 1px solid #ccc;
                border-radius: 4px;
            }
            
            .apply-button {
                margin-top: 15px;
                padding: 10px;
                background-color: #007bff;
                color: white;
                border: none;
                border-radius: 5px;
                cursor: pointer;
            }
            
            .apply-button:hover {
                background-color: #0069d9;
            }
            
            .variant-count-label {
                font-weight: bold;
                margin-bottom: 10px;
                padding: 8px;
                background-color: #f8f9fa;
                border-radius: 4px;
                display: inline-block;
            }
            
            .filter-download-button {
                margin-bottom: 10px;
                padding: 10px 15px;
                background-color: #28a745;
                color: white;
                border: none;
                border-radius: 5px;
                cursor: pointer;
                font-weight: bold;
            }
            
            .filter-download-button:hover {
                background-color: #218838;
            }
            
            @media (max-width: 992px) {
                .plot-container, .sidebar-container {
                    width: 100%;
                }
                .plot-container {
                    padding-right: 0;
                    margin-bottom: 20px;
                }
            }
        </style>
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
'''

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

def load_data(file_path, pval_threshold=1e-3):
    """
    Load GWAS data from a CSV file and apply initial filtering
    to improve performance
    
    Parameters:
    -----------
    file_path : str
        Path to the GWAS data file
    pval_threshold : float or None
        P-value threshold to pre-filter variants (default: 1e-3)
        If None, load all data without p-value filtering
    """
    import pandas as pd
    import numpy as np

    print(f"\n\nLoading data from {file_path}...")
    df = pd.read_csv(
        file_path,
        sep=",",
        compression="gzip",
        usecols=["chrom", "genpos", "pval", "VQSR", "HWE.ctrl", "batch.FEpval.FUS", 
                "ctrl_F_MISS", "case_F_MISS", "batch.pval.ctrl", "a1freq_cases", "a1freq_controls"]
    )

    # Pre-filter by p-value to improve performance, if threshold is provided
    print(f"Total variants before filtering: {len(df)}")
    if pval_threshold is not None:
        df = df[df["pval"] <= pval_threshold]
        print(f"Variants after p-value filtering (p ≤ {pval_threshold}): {len(df)}")
    else:
        print(f"Using all {len(df)} variants (no p-value pre-filtering)")

    # Keep original columns for hover
    df["original_chrom"] = df["chrom"]
    df["original_pos"] = df["genpos"]

    df["P"] = df["pval"]  # For compatibility with filter code
    df["-log10P"] = -1 * np.log10(df["pval"])

    # Calculate newPOS without using apply() to avoid potential issues with large dataframes
    df["CHROM_str"] = df["chrom"].astype(str)

    # Create a mapping of chromosome to cumulative length
    chrom_to_cumulative = {}
    for i in range(1, 24):
        chrom_to_cumulative[str(i)] = sum([chrom_to_length[str(j)] for j in range(1, i)])

    # Use vectorized operations instead of apply
    df["newPOS"] = df["genpos"] + df["chrom"].astype(str).map(chrom_to_cumulative)

    print("Data loaded successfully")
    return df.drop(columns=["chrom", "genpos", "pval"])  # Drop original columns but keep the copies

def get_peak_variants(df, significance_threshold=5e-8, distance_threshold=500_000):
    """
    Extract significant variants, keeping only the top variant in each peak region
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing GWAS results
    significance_threshold : float
        P-value threshold for significance (default: 5e-8)
    distance_threshold : int
        Distance (in bp) to consider variants in the same peak (default: 500kb)
    
    Returns:
    --------
    pandas.DataFrame
        DataFrame with one row per significant peak
    """
    # Filter to significant variants
    sig_variants = df[df["P"] <= significance_threshold].copy()

    if len(sig_variants) == 0:
        return pd.DataFrame()  # Return empty DataFrame if no significant variants

    # Sort by p-value (ascending)
    sig_variants = sig_variants.sort_values("P")

    # Check if original columns exist, otherwise use CHROM_str
    has_original_cols = "original_chrom" in df.columns and "original_pos" in df.columns
    chrom_col = "original_chrom" if has_original_cols else "CHROM_str"

    # Always use the original_pos for distance calculations if available,
    # otherwise fallback to newPOS (which might be less accurate for this purpose)
    pos_col = "original_pos" if has_original_cols else "newPOS"

    # Initialize list to store peak variants
    peak_variants = []

    # Track which variants have been assigned to a peak
    assigned_variants = set()

    for idx, variant in sig_variants.iterrows():
        if idx in assigned_variants:
            continue

        # Add this variant to the peak list
        peak_variants.append(variant)

        # Mark it as assigned
        assigned_variants.add(idx)

        # Find all variants within distance_threshold of this one
        chrom = variant[chrom_col]
        pos = variant[pos_col]

        # Find other variants on same chromosome within distance threshold
        same_chrom_variants = sig_variants[sig_variants[chrom_col] == chrom]
        nearby_variants = same_chrom_variants[same_chrom_variants[pos_col].between(pos - distance_threshold, pos + distance_threshold)]

        # Mark all these variants as assigned
        assigned_variants.update(nearby_variants.index)

    # Convert the list of peak variants to a DataFrame
    peak_df = pd.DataFrame(peak_variants)

    # Select only the columns we want for display
    # Always prioritize original chromosome and position for the table
    if has_original_cols:
        display_columns = [
            "original_chrom",
            "original_pos",
            "P",
            "-log10P",
            "VQSR",
            "HWE.ctrl",
            "batch.FEpval.FUS",
            "ctrl_F_MISS",
            "case_F_MISS",
            "a1freq_cases",
            "a1freq_controls",
            "batch.pval.ctrl",
        ]

        # Add OR column if it exists
        if "OR" in peak_df.columns:
            display_columns.append("OR")
            
        # Define column rename dictionary for original columns
        column_rename = {
            "original_chrom": "Chromosome",
            "original_pos": "Position",
            "P": "P-value",
            "-log10P": "-log10(P)",
        }

    else:
        display_columns = [
            "CHROM_str",
            "newPOS",
            "P",
            "-log10P",
            "VQSR",
            "HWE.ctrl",
            "batch.FEpval.FUS",
            "ctrl_F_MISS",
            "case_F_MISS",
            "a1freq_cases",
            "a1freq_controls",
            "batch.pval.ctrl",
        ]
        
        # Add OR column if it exists
        if "OR" in peak_df.columns:
            display_columns.append("OR")
            
        # Define column rename dictionary for fallback columns
        column_rename = {
            "CHROM_str": "Chromosome",
            "newPOS": "Position",
            "P": "P-value",
            "-log10P": "-log10(P)",
        }

    # Filter columns to only those that exist
    display_columns = [col for col in display_columns if col in peak_df.columns]

    # Return the DataFrame with selected columns and renamed
    # Only rename columns that exist in the dataframe
    rename_dict = {k: v for k, v in column_rename.items() if k in peak_df.columns}
    return peak_df[display_columns].rename(columns=rename_dict)

def load_full_data_with_all_columns(file_path):
    """
    Load GWAS data with all columns needed for final output
    """
    import pandas as pd
    import numpy as np

    print(f"Loading full dataset with all columns from {file_path}...")

    # Required columns for filtering
    filter_columns = ["chrom", "genpos", "pval", "VQSR", "HWE.ctrl", "batch.FEpval.FUS", 
                     "ctrl_F_MISS", "case_F_MISS", "batch.pval.ctrl", "a1freq_cases", "a1freq_controls"]

    # Required columns for final output
    output_columns = ["allele0", "allele1"]

    # Combine all required columns
    all_columns = filter_columns + output_columns

    # Try to load with OR column first, then beta column
    try:
        df = pd.read_csv(
            file_path,
            sep=",",
            compression="gzip",
            usecols=all_columns + ["OR"]
        )
        print("Loaded data with OR column")
        has_beta = False

    except ValueError as e:
        if "OR" in str(e):
            try:
                df = pd.read_csv(
                    file_path,
                    sep=",",
                    compression="gzip",
                    usecols=all_columns + ["beta"]
                )
                print("Loaded data with beta column")
                has_beta = True

            except ValueError as e2:
                if "beta" in str(e2):
                    print("Warning: Neither OR nor beta column found")
                    df = pd.read_csv(
                        file_path,
                        sep=",",
                        compression="gzip",
                        usecols=all_columns
                    )
                    has_beta = False
                else:
                    raise e2
        else:
            raise e

    print(f"Total variants loaded: {len(df)}")

    # Convert beta to OR if needed
    if has_beta:
        print("Converting beta to OR (OR = e^beta)")
        df["OR"] = np.exp(df["beta"])
        df = df.drop(columns=["beta"])

    # Keep original columns for output

    df["P"] = df["pval"]  # For compatibility with filter code

    print("Full data loaded successfully with all required columns")

    # Drop original columns but keep the copies
    return df.rename(columns={"allele0": "REF", 
                              "allele1": "ALT", 
                              "chrom": "CHROM", 
                              "genpos": "POS",
                              "pval": "P"})

if __name__ == "__main__":
    print(f"Starting GWAS Manhattan Plot Dashboard")
    print(f"Data file: {args.data_file}")
    print(f"Pre-filter p-value threshold: {args.pval_threshold}")
    print(f"Running on http://{args.host}:{args.port}")
    
    app.run(
        debug=args.debug,
        host=args.host,
        port=args.port
    )
