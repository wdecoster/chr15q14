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
        default="/home/wdecoster/Downloads/results.MAF0.01.VQSRpass.ctrlBatch.05.FUSBatch.05.ctrlHWE1E-8.ctrlMISS.05.caseMISS.05.forPublication.txt.gz",
    )
    return parser.parse_args()

def load_data(input_file):
    """Load GWAS data from a TSV file into a pandas DataFrame."""
    df = pd.read_csv(
        input_file,
        sep="\t",
        usecols=[
            "CHROM",
            "POS",
            "A1",
            "OMITTED",
            "OR",
            "LOG(OR)_SE",
            "A1_CTRL_FREQ",
            "P",
            "RSID"
        ],
    ).rename(columns={
        "CHROM": "chromosome",
        "POS": "base_pair_location",
        "A1": "effect_allele",
        "OMITTED": "other_allele",
        "OR": "odds_ratio",
        "LOG(OR)_SE": "standard_error",
        "A1_CTRL_FREQ": "effect_allele_frequency",
        "P": "p_value",
        "RSID": "rs_id"

    })[[
        "chromosome",
        "base_pair_location",
        "effect_allele",
        "other_allele",
        "odds_ratio",
        "standard_error",
        "effect_allele_frequency",
        "p_value",
        "rs_id",
    ]]
    return df

def filter_data(df):
    before = len(df)
    df = df[df["effect_allele"] != "*"]
    after = len(df)
    sys.stderr.write(f"Filtered out {before - after} rows with effect_allele='*'\n")
    df = df[df["other_allele"] != "*"]
    after2 = len(df)
    sys.stderr.write(f"Filtered out {after - after2} rows with other_allele='*'\n")
    return df

def main():
    args = get_args()
    df = load_data(args.input)
    df = filter_data(df)
    print(df.to_csv(sep="\t", index=False))

if __name__ == "__main__":
    main()
