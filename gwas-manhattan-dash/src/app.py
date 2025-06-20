from dash import Dash, dcc, html, Input, Output, State, ctx, dash_table
import pandas as pd
import dash_bootstrap_components as dbc
from components.filters import create_filter_layout
from components.manhattan_plot import generate_manhattan_plot
from utils.data_loader import load_data, filter_data, get_vqsr_values
import time

# Initialize the app with Bootstrap for styling
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Load the pre-filtered data for initial fast rendering
prefiltered_data = load_data("/home/wdecoster/local/CC.regenie.csv.gz", pval_threshold=1e-2)
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
            
            # Button to toggle between full and pre-filtered data
            html.Button(
                'Apply filters to ALL variants (slow)', 
                id='use-full-data-button',
                className='use-full-data-button',
                n_clicks=0
            ),
            
            dcc.Loading(
                id="loading",
                type="circle",
                children=dcc.Graph(id='manhattan-plot', className="manhattan-plot")
            ),
            
            # Add this section for the significant peaks table
            html.H3("Significant Peaks (one variant per peak)", className="mt-4"),
            dcc.Loading(
                id="table-loading",
                type="circle",
                children=html.Div(id='significant-variants-table')
            ),
            
            # Store to keep track of whether we're using full data
            dcc.Store(id='using-full-data', data=False)
        ], className="plot-container"),
        
        # Filters container (right side, 1/3 width)
        html.Div([
            create_filter_layout(vqsr_values)
        ], className="sidebar-container")
    ], className="main-container")
], className="app-container")

# Callback to update the plot when filters are applied
@app.callback(
    [Output('manhattan-plot', 'figure'),
     Output('variant-count', 'children'),
     Output('use-full-data-button', 'children'),
     Output('using-full-data', 'data'),
     Output('significant-variants-table', 'children')],
    [Input('apply-filters-button', 'n_clicks'),
     Input('use-full-data-button', 'n_clicks')],
    [State('vqsr-filter', 'value'),
     State('hwe-ctrl-filter', 'value'),
     State('batch-fepval-fus-filter', 'value'),
     State('ctrl-f-miss-filter', 'value'),
     State('case-f-miss-filter', 'value'),
     State('allele-freq-filter', 'value'),  # Updated ID
     State('batch-pval-ctrl-filter', 'value'),
     State('using-full-data', 'data')],
    prevent_initial_call=False
)
def update_plot(apply_clicks, full_data_clicks, vqsr_values, hwe_ctrl, 
                batch_fepval_fus, ctrl_f_miss, case_f_miss, allele_freq, batch_pval_ctrl, using_full_data):
    global full_data
    
    # Initialize button text based on current state
    button_text = "Using ALL variants - Click to switch back" if using_full_data else "Apply filters to ALL variants (slow)"
    
    # Determine which button was clicked
    triggered_id = ctx.triggered_id if ctx.triggered_id else 'apply-filters-button'
    
    # Toggle full data mode if the full-data button was clicked
    if triggered_id == 'use-full-data-button' and full_data_clicks > 0:
        using_full_data = not using_full_data
    
    # Load full data if needed and not already loaded
    if using_full_data and full_data is None:
        start_time = time.time()
        full_data = load_data("/home/wdecoster/local/CC.regenie.csv.gz", pval_threshold=None)
        load_time = time.time() - start_time
        print(f"Full data loaded in {load_time:.2f} seconds")
    
    # Choose which dataset to use
    data_to_use = full_data if using_full_data else prefiltered_data
    
    # Apply filters to the data
    try:
        filtered_data = filter_data(
            data_to_use,
            vqsr_values=vqsr_values,
            hwe_ctrl_cutoff=float(hwe_ctrl) if hwe_ctrl else None,
            batch_fepval_fus_cutoff=float(batch_fepval_fus) if batch_fepval_fus else None,
            ctrl_f_miss_cutoff=float(ctrl_f_miss) if ctrl_f_miss else None,
            case_f_miss_cutoff=float(case_f_miss) if case_f_miss else None,
            allele_freq_cutoff=float(allele_freq) if allele_freq else None,  # Updated parameter name
            batch_pval_ctrl_cutoff=float(batch_pval_ctrl) if batch_pval_ctrl else None,
        )
                
        # Generate the variant count label
        dataset_type = "all" if using_full_data else "pre-filtered"
        variant_count_text = f"Displaying {len(filtered_data):,} variants (from {dataset_type} data)"
        
        # Update the button text based on the current state
        button_text = "Using ALL variants - Click to switch back" if using_full_data else "Apply filters to ALL variants (slow)"
        
        # Get significant peak variants
        peak_variants = get_peak_variants(filtered_data)
        
        # Create the table for significant variants
        if len(peak_variants) > 0:
            peak_table = dash_table.DataTable(
                data=peak_variants.round(10).to_dict('records'),
                columns=[{"name": i, "id": i} for i in peak_variants.columns],
                style_table={'overflowX': 'auto'},
                style_cell={
                    'textAlign': 'center',
                    'padding': '10px',
                    'minWidth': '100px'
                },
                style_header={
                    'backgroundColor': 'rgb(230, 230, 230)',
                    'fontWeight': 'bold',
                    'textAlign': 'center'
                },
                style_data_conditional=[
                    {
                        'if': {'column_id': 'P-value'},
                        'fontWeight': 'bold'
                    }
                ],
                sort_action='native',
                filter_action='native',
                page_size=10,
                # Add these properties to enable selection and copying:
                cell_selectable=True,
                row_selectable=False,
                selected_cells=[],
                # Enable copy to clipboard
                export_format='xlsx',  # Add export button (optional)
                export_headers='display',
                # Enable right-click menu with copy option
                css=[{
                    'selector': '.dash-cell div.dash-cell-value', 
                    'rule': 'display: inline; white-space: inherit; overflow: inherit; text-overflow: inherit;'
                }]
            )
        else:
            peak_table = html.Div("No significant variants found", style={'color': 'gray', 'textAlign': 'center'})
        
        # Generate and return the plot, variant count, button text, and updated full data state
        return generate_manhattan_plot(filtered_data), variant_count_text, button_text, using_full_data, peak_table
    
    except Exception as e:
        import traceback
        print(f"Error: {e}")
        print(traceback.format_exc())
        empty_table = html.Div("Error generating table: " + str(e), style={'color': 'red'})
        return {}, f"Error: {str(e)}", button_text, using_full_data, empty_table

# Add CSS for layout positioning
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
            
            .use-full-data-button {
                margin-left: 10px;
                margin-bottom: 10px;
                padding: 8px;
                background-color: #28a745;
                color: white;
                border: none;
                border-radius: 4px;
                cursor: pointer;
            }
            
            .use-full-data-button:hover {
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

def get_peak_variants(df, significance_threshold=5e-8, distance_threshold=500000):
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
        nearby_variants = same_chrom_variants[
            (same_chrom_variants[pos_col] >= pos - distance_threshold) &
            (same_chrom_variants[pos_col] <= pos + distance_threshold)
        ]
        
        # Mark all these variants as assigned
        assigned_variants.update(nearby_variants.index)
    
    # Convert the list of peak variants to a DataFrame
    peak_df = pd.DataFrame(peak_variants)
    
    # Select only the columns we want for display
    # Always prioritize original chromosome and position for the table
    if has_original_cols:
        display_columns = [
            "original_chrom", "original_pos", "P", "-log10P", 
            "VQSR", "HWE.ctrl", "batch.FEpval.FUS", "ctrl_F_MISS", "case_F_MISS", 
            "a1freq_cases", "a1freq_controls", "batch.pval.ctrl"  # Added both allele frequency columns
        ]
        column_rename = {
            "original_chrom": "Chromosome",
            "original_pos": "Position",
            "P": "P-value",
            "-log10P": "-log10(P)",
        }
    else:
        display_columns = [
            "CHROM_str", "newPOS", "P", "-log10P", 
            "VQSR", "HWE.ctrl", "batch.FEpval.FUS", "ctrl_F_MISS", "case_F_MISS", 
            "a1freq_cases", "a1freq_controls", "batch.pval.ctrl"  # Added both allele frequency columns
        ]
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

if __name__ == "__main__":
    app.run(debug=True)
