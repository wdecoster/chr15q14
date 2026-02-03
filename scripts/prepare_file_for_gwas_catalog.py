import pandas as pd
import sys

def get_args():
    from argparse import ArgumentParser
    parser= ArgumentParser(
        description="Prepare GWAS summary statistics file for submission to GWAS Catalog"
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="Path to the input GWAS summary statistics file (TSV format)",
        default="/home/wdecoster/Downloads/FTLD-FUS_allCTRL.CC.regenie.csv.gz",
    )
    return parser.parse_args()

def load_data(input_file):
    """Load GWAS data from a TSV file into a pandas DataFrame."""
    df = pd.read_csv(
        input_file,
        sep=",",
        usecols=[
            "chrom",
            "genpos",
            "pval",
            "allele0",
            "allele1",
            "beta",
            "se",
            "VQSR",
            "HWE.ctrl",
            "batch.FEpval.FUS",
            "ctrl_F_MISS",
            "case_F_MISS",
            "batch.pval.ctrl",
            "a1freq_cases",
            "a1freq_controls",
            "RSID",
        ],
    ).rename(
        columns={
            "chrom": "chromosome",
            "genpos": "base_pair_location",
            "allele0": "effect_allele",
            "allele1": "other_allele",
            "beta": "beta",
            "se": "standard_error",
            "a1freq_controls": "effect_allele_frequency",
            "pval": "p_value",
            "RSID": "rs_id",
        }
    )
    return df

def filter_data(
    df,
    vqsr_values=["PASS"],
    hwe_ctrl_cutoff=1e-6,
    batch_fepval_fus_cutoff=0.01,
    ctrl_f_miss_cutoff=0.1,
    case_f_miss_cutoff=0.1,
    batch_pval_ctrl_cutoff=0.01,
    allele_freq_cutoff=0.01,
):
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
        filtered_df = filtered_df[
            filtered_df["batch.FEpval.FUS"] > batch_fepval_fus_cutoff
        ]

    # Apply ctrl_F_MISS filter (numeric, above cutoff)
    if ctrl_f_miss_cutoff is not None:
        filtered_df = filtered_df[filtered_df["ctrl_F_MISS"] < ctrl_f_miss_cutoff]

    # Apply case_F_MISS filter (numeric, above cutoff)
    if case_f_miss_cutoff is not None:
        filtered_df = filtered_df[filtered_df["case_F_MISS"] < case_f_miss_cutoff]

    # Apply batch.pval.ctrl filter (numeric, below cutoff)
    if batch_pval_ctrl_cutoff is not None:
        filtered_df = filtered_df[
            filtered_df["batch.pval.ctrl"] > batch_pval_ctrl_cutoff
        ]

    # Apply allele frequency filter - keep variants with frequency above cutoff in EITHER cases OR controls
    if allele_freq_cutoff is not None:
        filtered_df = filtered_df[
            (filtered_df["a1freq_cases"] > allele_freq_cutoff)
            | (filtered_df["effect_allele_frequency"] > allele_freq_cutoff)
        ]

    before = len(filtered_df)
    filtered_df = filtered_df[filtered_df["effect_allele"] != "*"]
    after = len(filtered_df)
    sys.stderr.write(f"Filtered out {before - after} rows with effect_allele='*'\n")
    filtered_df = filtered_df[filtered_df["other_allele"] != "*"]
    after2 = len(filtered_df)
    sys.stderr.write(f"Filtered out {after - after2} rows with other_allele='*'\n")
    return filtered_df[
        [
            "chromosome",
            "base_pair_location",
            "effect_allele",
            "other_allele",
            "beta",
            "standard_error",
            "effect_allele_frequency",
            "p_value",
            "rs_id",
        ]
    ]

def main():
    args = get_args()
    df = load_data(args.input)
    df = filter_data(df)
    print(df.to_csv(sep="\t", index=False))

if __name__ == "__main__":
    main()
