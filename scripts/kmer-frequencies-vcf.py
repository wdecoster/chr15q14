from argparse import ArgumentParser
from cyvcf2 import VCF
from collections import Counter
import pandas as pd
import plotly.express as px
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

        df["spacer"] = 1
        df["spacer2"] = 1
        df["spacer3"] = 1
        # drop the Group column
        df = df.drop(columns=["Group", "haplotype"])

    df = df.fillna(0).set_index("individual")
    df["CTCTCTCTCTCT"] = df["CTCTCTCTCTCT"].apply(lambda x: x if x > 0.01 else 0)
    df.to_csv(args.counts if args.counts else f"kmer{args.kmer}-counts.tsv", sep="\t")
    if args.nosort:
        plot_heatmap(df.transpose(), k=args.kmer, outputfile=args.output, args=args)
    else:
        plot_heatmap(
            df.sort_values(
                by=[
                    "CT" * int(args.kmer / 2),
                    "CCTT" * int(args.kmer / 4),
                    "CTTT" * int(args.kmer / 4),
                    "CCCT" * int(args.kmer / 4),
                ],
            ).transpose(),
            k=args.kmer,
            outputfile=args.output,
            args=args,
        )


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


def plot_heatmap(df, k, outputfile, args, max_missing=0.1):
    color_scale = [(0, "white"), (1, "black")]
    # only keep rows that are not < 0.01 for too many samples
    mask1 = (df < 0.01).sum(axis=1) < ((1 - max_missing) * len(df.columns))
    # but keep also rows that are above 0.2 for at least one sample
    mask2 = (df > 0.2).sum(axis=1) > 0

    df = df.loc[mask1 | mask2, :]

    # change the row with label "spacer" to 0
    # I set it to 1 to make sure it is not removed by the mask
    # this is quite ridiculous but hey bear with me
    df.loc["spacer", :] = 0
    df.loc["spacer2", :] = 0
    df.loc["spacer3", :] = 0

    kmers = [
        c
        for c in df.index
        if c not in ["spacer2", "spacer", "spacer3", "aFTLD-U", "Haplotype"]
    ]

    # change the order of the kmers to show CT dimer first and CCTT tetramer second
    df = df.reindex(
        [
            "Haplotype",
            "spacer",
            "aFTLD-U",
            "spacer2",
            "spacer3",
            "CTCTCTCTCTCT",
            "CCTTCCTTCCTT",
        ]
        + [k for k in kmers if k not in ["CTCTCTCTCTCT", "CCTTCCTTCCTT"]]
    )

    rename_dict = {
        "spacer": "",
        "spacer2": " ",
        "spacer3": "  ",
        "aFTLD-U": "Phenotype",
        "Haplotype": "Haplotype",
        "CTCTCTCTCTCT": "(CT)<sub>6</sub>",
        "CTTTCTTTCTTT": "(CTTT)<sub>3</sub>",
        "CCTTCCTTCCTT": "(CCTT)<sub>3</sub>",
        "CCCTCTCCCTCT": "(CCCTCT)<sub>2</sub>",
        "CCCTCCCTCCCT": "(CCCT)<sub>3</sub>",
        "CCCCCCCCCCCC": "C<sub>12</sub>",
        "ATATATATATAT": "(AT)<sub>6</sub>",
        "GTGTGTGTGTGT": "(GT)<sub>6</sub>",
    }

    fig = px.imshow(
        df.transpose(),
        x=[rename_dict.get(c, c) for c in df.index],
        labels=dict(
            x="",
            y="Individuals with an expanded repeat allele",
            color="Fraction",
        ),
        color_continuous_scale=color_scale,
        aspect="auto",
    )
    fig.update_xaxes(
        tickfont_size=20,
        tickangle=-90,
    )
    fig.update_yaxes(showticklabels=False)

    title = (
        f"Repeat composition (individual reads)"
        if args.somatic
        else "Repeat consensus sequence composition"
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
    fig.update_xaxes(showline=True, linewidth=2, linecolor="black", mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor="black", mirror=True)
    plotname = outputfile if outputfile else f"kmer{k}-heatmap.html"
    with open(plotname, "w") as output:
        output.write(fig.to_html(include_plotlyjs="cdn"))


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
    return parser.parse_args()


if __name__ == "__main__":
    main()
