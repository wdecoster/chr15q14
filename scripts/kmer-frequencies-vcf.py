from argparse import ArgumentParser
from cyvcf2 import VCF
from collections import Counter
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from os.path import basename
import sys
import numpy as np


def main():
    df = pd.DataFrame()
    args = get_args()
    if args.somatic:
        for vcff in args.vcf:
            kmers = Counter()
            vcf = VCF(vcff)
            name = basename(vcff).replace(".vcf.gz", "").replace("_FCX", "")
            for v in vcf:
                allele_used = []
                # this code is just adding the kmers from both alleles together, which could be suboptimal in the case of two different expansions on both haplotypes
                # however, this is rather rare, and most typically those are incorrectly phased alleles due to large somatic differences
                # I checked this by making an astronaut plot for the samples that triggered the warning below, and only one sample is problematic
                for allele in [0, 1]:
                    genotype = v.genotypes[0][allele]
                    if genotype > 0 and len(v.ALT[genotype - 1]) > args.minlength:
                        for counts in [
                            count_kmers(seq, k=args.kmer)
                            for seq in v.INFO.get("SEQS").split(",")[allele].split(":")
                            if len(seq) > args.minlength
                        ]:
                            kmers.update(counts)
                        allele_used.append(allele)
                if len(allele_used) == 2:
                    sys.stderr.write(f"Warning: more than one allele used for {name}\n")
                pruned_kmers = prune_counts(kmers)
                if pruned_kmers:
                    pruned_kmers.update({"individual": name})
                    pruned_kmers_df = pd.DataFrame.from_dict(
                        pruned_kmers, orient="index"
                    ).transpose()
                    df = pd.concat([df, pruned_kmers_df], ignore_index=True)
                break
    else:
        for vcff in args.vcf:
            kmers = Counter()
            vcf = VCF(vcff)
            name = basename(vcff).replace(".vcf.gz", "").replace("_FCX", "")
            for v in vcf:
                allele_used = []
                for allele in [0, 1]:
                    genotype = v.genotypes[0][allele]
                    if genotype > 0 and len(v.ALT[genotype - 1]) > args.minlength:
                        counts = count_kmers(v.ALT[genotype - 1], k=args.kmer)
                        kmers.update(counts)
                        allele_used.append(allele)
                if len(allele_used) == 2:
                    sys.stderr.write(f"Warning: more than one allele used for {name}\n")
                pruned_kmers = prune_counts(kmers)
                if pruned_kmers:
                    pruned_kmers.update({"individual": name})
                    pruned_kmers_df = pd.DataFrame.from_dict(
                        pruned_kmers, orient="index"
                    ).transpose()
                    df = pd.concat([df, pruned_kmers_df], ignore_index=True)
                break
    if args.sampleinfo:
        sampleinfo = pd.read_table(
            args.sampleinfo, usecols=["individual", "cohort", "haplotype"]
        ).rename(columns={"cohort": "Group"})
        df = df.merge(sampleinfo, on="individual", how="left")
        df["Group"] = df["Group"].fillna("unknown/misc")
        # add a 'aFTLD-U' column which is 1 if the sample is aFTLD-U and 0 otherwise
        df["aFTLD-U"] = df["Group"].apply(lambda x: 1 if x == "aFTLD-U" else 0)
        # the haplotype has to be encoded as 1 for major, 0.5 for minor and 0 for the rest
        df["Haplotype"] = 0
        df.loc[df["haplotype"] == "minor", "Haplotype"] = 0.5
        df.loc[df["haplotype"] == "major", "Haplotype"] = 1

        # drop the Group column
        df = df.drop(columns=["Group", "haplotype"])

    df = df.fillna(0).set_index("individual")
    df["CTCTCTCTCTCT"] = df["CTCTCTCTCTCT"].apply(lambda x: x if x > 0.01 else 0)
    df.to_csv(args.counts if args.counts else f"kmer{args.kmer}-counts.tsv", sep="\t")
    plot_heatmap(df, args=args)


def count_kmers(seq, k=4):
    kmers = Counter()
    for i in range(len(seq) - k + 1):
        kmers[seq[i : i + k]] += 1
    return kmers


def get_rotations(kmer):
    """
    Rotate a kmer to get all equivalent representations
    """
    e = len(kmer)
    rotations = [kmer[i:e] + kmer[:i] for i in range(e)]
    return sorted(rotations)[0], rotations


def prune_counts(kmers):
    """
    For all rotations of a kmer, keep only the lexicographical first
    Return the number as a fraction of the total kmers
    And only return those kmers that are above 1%
    """
    pruned = dict()
    for key in kmers:
        first, rotations = get_rotations(key)
        if first in pruned.keys():
            continue
        else:
            pruned[first] = sum([kmers[r] for r in rotations])
    total_kmers = sum(pruned.values())
    return {k: v / total_kmers for k, v in pruned.items()}


def plot_heatmap(df, args, max_missing=0.1):
    if args.nosort:
        df = df.transpose()
    else:
        df = df.sort_values(
            by=[
                "CT" * int(args.kmer / 2),
                "CCTT" * int(args.kmer / 4),
                "CTTT" * int(args.kmer / 4),
                "CCCT" * int(args.kmer / 4),
            ],
            ascending=False,
        ).transpose()
    # only keep rows that are not < 0.01 for too many samples
    mask1 = (df < 0.01).sum(axis=1) < ((1 - max_missing) * len(df.columns))
    # but keep also rows that are above 0.2 for at least one sample
    mask2 = (df > 0.2).sum(axis=1) > 0

    df = df.loc[mask1 | mask2, :]

    kmers = [
        c
        for c in df.index
        if c not in ["aFTLD-U", "Haplotype"]
    ]

    # change the order of the kmers to show CT dimer first and CCTT tetramer second
    df = df.reindex(
        ["CTCTCTCTCTCT", "CCTTCCTTCCTT",]
        + [k for k in kmers if k not in ["CTCTCTCTCTCT", "CCTTCCTTCCTT"]]
        + ["aFTLD-U", "Haplotype"]
    )

    rename_dict = {
        # "aFTLD-U": "Phenotype",
        # "Haplotype": "Haplotype",
        "CTCTCTCTCTCT": "(CT)<sub>6</sub>",
        "CTTTCTTTCTTT": "(CTTT)<sub>3</sub>",
        "CCTTCCTTCCTT": "(CCTT)<sub>3</sub>",
        "CCCTCTCCCTCT": "(CCCTCT)<sub>2</sub>",
        "CCCTCCCTCCCT": "(CCCT)<sub>3</sub>",
        "CCCCCCCCCCCC": "C<sub>12</sub>",
        "ATATATATATAT": "(AT)<sub>6</sub>",
        "GTGTGTGTGTGT": "(GT)<sub>6</sub>",
    }

    df = df.transpose()

    from plotly.subplots import make_subplots
    fig = make_subplots(rows=1, cols=3, column_widths=[0.1, 0.1, 0.8], shared_yaxes=True, horizontal_spacing=0.01)

    kmer_df = df[
        [k for k in df.columns if k not in ["aFTLD-U", "Haplotype"]]
    ]
    info_df = df[["aFTLD-U", "Haplotype"]]

    fig.add_trace(
        go.Heatmap(
            z=kmer_df.values,
            y=kmer_df.index,
            x=[rename_dict.get(c, c) for c in kmer_df.columns],
            # labels=dict(
            #     x="",
            #     y="Individuals with an expanded repeat allele",
            #     color="Fraction",
            # ),
            colorscale=[(0, "white"), (1, "black")],
            colorbar=dict(
                title="Fraction",
                nticks=5,
                tickmode="auto",
                len=0.75
            ),
        ),
        row=1,
        col=3,
    )
    fig.add_trace(go.Heatmap(
        z=info_df.loc[:, "aFTLD-U"].values.reshape(-1, 1),
        y=info_df.index,
        x=["Phenotype"],
        colorscale=[(0, "black"), (1, "red")],
        showscale=False,
        ),
        row=1, col=1,
        )
    fig.add_trace(go.Heatmap(
        z=info_df.loc[:, "Haplotype"].values.reshape(-1, 1),
        y=info_df.index,
        x=["Haplotype"],
        colorscale=[(0, "blue"), (0.5, "orange"), (1, "red")],
        showscale=False,
        ),
        row=1, col=2,
        )

    title = (
        f"Repeat composition (individual reads)"
        if args.somatic
        else "Repeat composition"
    )

    fig.update_layout(
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
        height=800,
        width=600,
        font=dict(size=18),
        title=title,
        xaxis={"dtick": 1},
    )

    fig.update_xaxes(tickfont_size=20, tickangle=-90, showline=True, linewidth=2, linecolor="black")
    fig.update_yaxes(showticklabels=False, title="Individuals with an expanded repeat", showline=True, linewidth=2, linecolor="black", row=1, col=1)
    plotname = args.output if args.output else f"kmer{args.kmer}-heatmap.html"
    with open(plotname, "w") as output:
        output.write(fig.to_html(include_plotlyjs="cdn"))
    if args.pdf:
        fig.write_image(plotname.replace('.html', '.pdf'))


def get_args():
    parser = ArgumentParser(description="sum insertions")
    parser.add_argument(
        "vcf", help="vcf files from STRdust, assuming a single target", nargs="+"
    )
    parser.add_argument("-k", "--kmer", help="kmer to count", default=4, type=int)
    parser.add_argument(
        "--nosort", help="sort kmers before plotting", action="store_true"
    )
    parser.add_argument("-o", "--output", help="output html plot")
    parser.add_argument("-c", "--counts", help="output kmer counts table")
    parser.add_argument(
        "-m",
        "--minlength",
        help="minimal length of an entry to be used",
        type=int,
        default=100,
    )
    parser.add_argument(
        "--somatic",
        help="plot kmers from individual sequences, if not use consensus sequence",
        action="store_true",
    )
    parser.add_argument("--sampleinfo", help="excel file with sample information")
    parser.add_argument("--pdf", help="Additionally create pdf output", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    main()
