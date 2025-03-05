import pandas as pd
from argparse import ArgumentParser
import plotly.express as px


def main():
    args = get_args()
    df = pd.read_csv(args.input, sep="\t", header=None, names=["sample", "copy_number"])
    if args.sampleinfo:
        sampleinfo = pd.read_excel(
            args.sampleinfo, usecols=["Gentli_ID", "PathDx1"]
        ).rename(columns={"Gentli_ID": "sample", "PathDx1": "Group"})
        df = df.merge(sampleinfo, on="sample", how="left").fillna("unknown/misc")
        # if a group has fewer than 3 samples, merge it with the unknown group
        counts = df["Group"].value_counts()
        for group in counts[counts <= 5].index:
            df.loc[df["Group"] == group, "Group"] = "unknown/misc"
    else:
        df["Group"] = "unknown"
    fig = px.strip(
        df,
        x="Group",
        y="copy_number",
        title="Copy number",
        hover_data=["sample"],
        color="Group",
    )
    html = fig.to_html()
    with open(args.output, "w") as f:
        f.write(html)


def get_args():
    parser = ArgumentParser(
        "Plot copy numbers from multiple samples for a single region"
    )
    parser.add_argument(
        "-i", "--input", help="file with copy numbers for multiple samples"
    )
    parser.add_argument(
        "-o", "--output", help="output file name", default="copy_numbers.html"
    )
    parser.add_argument("--sampleinfo", help="excel file with sample information")
    return parser.parse_args()


if __name__ == "__main__":
    main()
