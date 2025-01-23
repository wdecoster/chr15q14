from argparse import ArgumentParser
import pandas as pd

def main():
    args = get_args()
    df = pd.read_csv(args.input, sep="\t")
    df = df[df["haplotype"] != "none"]
    haplotype_alias = {'major': 'A', 'minor': 'B'}
    df["haplotype"] = df["haplotype"].apply(lambda x: haplotype_alias[x])
    df = df.sort_values(by=["haplotype", "group", "CT_dimer_count"], ascending=[True, True, False]).reset_index(drop=True)
    df["number"] = df.index + 1

    # introduce anonymous aliases for the individuals, with PAT1 etc for patients and CON1 etc for controls
    # but keep 1000G samples with their original names
    # add a column with the alias, including a number to distinguish between individuals with the same group
    group_alias = {
        "1000G": "1000G",
        "in-house non-aFTLD-U": "CON",
        "aFTLD-U": "aFTLD-U",
    }
    df["alias"] = df.apply(lambda x: group_alias[x["group"]] + "_" + str(x["number"]) if x["group"] != "1000G" else x["name"], axis=1)
    df["classifier1"] = df.apply(lambda x: "pathogenic" if x["length"] >= 450 and x["%CT"] >= 0.8 else "", axis=1)
    df["classifier2"] = df["CT_dimer_count"].apply(lambda x: "pathogenic" if x >= 190 else "")
    columns = ["name", "alias", "group", "haplotype", "copy number", "source tissue", "length", "%CT", "%CCCTCT", "CT_dimer_count", "classifier1", "classifier2"]
    df[columns].to_excel(args.output, index=False)


def get_args():
    parser = ArgumentParser("Make a summary table for publication of haplotype carriers (controls and patients)")
    parser.add_argument("-i", "--input", help="Input overview file generated in the workflow")
    parser.add_argument("-o", "--output", help="Output file name", default="haplotype_carriers.xlsx")
    return parser.parse_args()


if __name__ == "__main__":
    main()
