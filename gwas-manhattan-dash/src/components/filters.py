from dash import html, dcc

def create_filter_layout(vqsr_values=None):
    """Create the filter layout for the dashboard"""
    
    # Create VQSR checkbox options based on data values
    if vqsr_values:
        vqsr_options = [{'label': str(value), 'value': str(value)} for value in vqsr_values]
        # Set only PASS as checked by default if available
        default_vqsr = ["PASS"] if "PASS" in vqsr_values else []
    else:
        # Fallback options if no values provided
        vqsr_options = [{'label': 'PASS', 'value': 'PASS'}]
        default_vqsr = ['PASS']
    
    return html.Div([
        html.H3("Filters", className="mb-3"),
        
        # VQSR Filter (Categorical, Checkboxes)
        html.Div([
            html.Label("VQSR Filter:", className="filter-label"),
            dcc.Checklist(
                id='vqsr-filter',
                options=vqsr_options,
                value=default_vqsr,
                labelStyle={'display': 'block'}
            )
        ], className="filter-group"),
        
        # HWE.ctrl Filter
        html.Div([
            html.Label("HWE.ctrl Filter (remove below):", className="filter-label"),
            dcc.Input(
                id='hwe-ctrl-filter',
                type='number',
                placeholder='e.g. 1E-8',
                value=1e-6,
                step=1e-9,
                className="input-field"
            )
        ], className="filter-group"),
        
        # batch.FEpval.FUS Filter
        html.Div([
            html.Label("batch.FEpval.FUS Filter (remove below):", className="filter-label"),
            dcc.Input(
                id='batch-fepval-fus-filter',
                type='number',
                placeholder='e.g. 0.05',
                value=0.01,
                step=0.01,
                className="input-field"
            )
        ], className="filter-group"),
        
        # ctrl_F_MISS Filter
        html.Div([
            html.Label("ctrl_F_MISS Filter (remove above):", className="filter-label"),
            dcc.Input(
                id='ctrl-f-miss-filter',
                type='number',
                placeholder='e.g. 0.05',
                value=0.1,
                step=0.01,
                min=0,
                max=1,
                className="input-field"
            )
        ], className="filter-group"),
        
        # case_F_MISS Filter (new)
        html.Div([
            html.Label("case_F_MISS Filter (remove above):", className="filter-label"),
            dcc.Input(
                id='case-f-miss-filter',
                type='number',
                placeholder='e.g. 0.05',
                value=0.1,  # Same default as ctrl_F_MISS
                step=0.01,
                min=0,
                max=1,
                className="input-field"
            )
        ], className="filter-group"),
        
        # batch.pval.ctrl Filter
        html.Div([
            html.Label("batch.pval.ctrl Filter (remove below):", className="filter-label"),
            dcc.Input(
                id='batch-pval-ctrl-filter',
                type='number',
                placeholder='e.g. 0.05',
                value=0.01,
                step=0.01,
                className="input-field"
            )
        ], className="filter-group"),
        
        # Apply button
        html.Button(
            'Apply Filters', 
            id='apply-filters-button', 
            className="apply-button",
            n_clicks=0
        ),
        
    ], className="filter-container")
