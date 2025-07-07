import pandas as pd
import numpy as np

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

    # Keep original columns for hover and table (rename them for clarity)
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
    
    # Don't drop original_chrom and original_pos since we need them
    return df.drop(columns=["chrom", "genpos", "pval"])

def filter_data(df, vqsr_values=None, hwe_ctrl_cutoff=None, 
                batch_fepval_fus_cutoff=None, ctrl_f_miss_cutoff=None,
                case_f_miss_cutoff=None, batch_pval_ctrl_cutoff=None,
                allele_freq_cutoff=None):
    """Filter the data based on the specified criteria"""
    filtered_df = df.copy()
    
    # Apply VQSR filter (categorical)
    if vqsr_values and len(vqsr_values) > 0:
        filtered_df = filtered_df[filtered_df["VQSR"].isin(vqsr_values)]
    
    # Apply HWE.ctrl filter (numeric, below cutoff)
    if hwe_ctrl_cutoff is not None:
        filtered_df = filtered_df[filtered_df["HWE.ctrl"] > hwe_ctrl_cutoff]
    
    # Apply batch.FEpval.FUS filter (numeric, below cutoff)
    if batch_fepval_fus_cutoff is not None:
        filtered_df = filtered_df[filtered_df["batch.FEpval.FUS"] > batch_fepval_fus_cutoff]
    
    # Apply ctrl_F_MISS filter (numeric, above cutoff)
    if ctrl_f_miss_cutoff is not None:
        filtered_df = filtered_df[filtered_df["ctrl_F_MISS"] < ctrl_f_miss_cutoff]
    
    # Apply case_F_MISS filter (numeric, above cutoff)
    if case_f_miss_cutoff is not None:
        filtered_df = filtered_df[filtered_df["case_F_MISS"] < case_f_miss_cutoff]
    
    # Apply batch.pval.ctrl filter (numeric, below cutoff)
    if batch_pval_ctrl_cutoff is not None:
        filtered_df = filtered_df[filtered_df["batch.pval.ctrl"] > batch_pval_ctrl_cutoff]
    
    # Apply allele frequency filter - keep variants with frequency above cutoff in EITHER cases OR controls
    if allele_freq_cutoff is not None:
        filtered_df = filtered_df[
            (filtered_df["a1freq_cases"] > allele_freq_cutoff) | 
            (filtered_df["a1freq_controls"] > allele_freq_cutoff)
        ]
    
    return filtered_df

def get_vqsr_values(df):
    """Get unique VQSR values from the dataframe"""
    return sorted(df["VQSR"].unique().tolist())