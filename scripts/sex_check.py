import pandas as pd
from argparse import ArgumentParser
import plotly.express as px
import numpy as np


def main():
    args = get_args()
    df = pd.read_csv(args.cramino, sep="\t", usecols=["identifier", "name", "chrX", "chrY"])
    # add a little bit of random noise to avoid overlapping points
    df["chrX"] = df["chrX"] + np.random.normal(-0.02, 0.02, len(df))
    df["chrY"] = df["chrY"] + np.random.normal(-0.02, 0.02, len(df))
    sampleinfo = pd.read_excel(args.sampleinfo, usecols=["Gentli_ID", "SEXatBirth"]).rename(
        columns={"Gentli_ID": "name"}
    )
    df = df.merge(sampleinfo, on="name", how="left").fillna("unknown")
    fig = px.scatter(df, x="chrX", y="chrY", color="SEXatBirth", hover_data=["name", "identifier"])
    # change opacity to 0.6
    fig.update_traces(opacity=0.8, marker_size=3)
    fig.update_layout(scattermode="group", scattergap=0.75)
    html = fig.to_html()
    with open(args.output, "w") as f:
        f.write(html)


def get_args():
    parser = ArgumentParser("Plot copy numbers from multiple samples for a single region")
    parser.add_argument("-c", "--cramino", help="summary file of cramino")
    parser.add_argument("--sampleinfo", help="excel file with sample information")
    parser.add_argument("-o", "--output", help="output file name", default="sex_check.html")
    return parser.parse_args()


if __name__ == "__main__":
    main()
