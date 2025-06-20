#!/usr/bin/env python3
from argparse import ArgumentParser
import sys
import pandas as pd
from cyvcf2 import VCF
import scipy.stats as stats
import numpy as np


def main():
    args = parse_arguments()

    (
        positive_indices,
        negative_indices,
        con_hapA_indices,
        patient_indices,
        control_indices,
    ) = define_groups(args)
    gwas_tested = parse_gwas(args.gwas) if args.gwas else []

    # Print header
    varinfo = "SVLEN" if args.type == "sv" else "VARIANT"
    print(
        f"CHROM\tPOS\t{varinfo}\tPOS_COUNT\tNEG_COUNT\tCON_HAPA_COUNT\tPAT_COUNT\tCTRL_COUNT\t"
        f"HAPA_FISHER_P\tHAPA_ODDS\tPAT_FISHER_P\tPAT_ODDS\tSEEN_IN_GWAS"
    )

    # Iterate through variants
    variants_kept = 0
    pos_total = len(positive_indices)
    neg_total = len(negative_indices)
    pat_total = len(patient_indices)
    ctrl_total = len(control_indices)

    for variant in VCF(args.vcf):
        # Get genotypes (0 = ref, 1 = het, 2 = unknown, 3 = hom_alt)
        gt_types = variant.gt_types

        # Count carriers in each group
        pos_carriers = sum(1 for i in positive_indices if gt_types[i] in [1, 3])
        neg_carriers = sum(1 for i in negative_indices if gt_types[i] in [1, 3])
        third_carriers = sum(1 for i in con_hapA_indices if gt_types[i] in [1, 3])
        pat_carriers = sum(1 for i in patient_indices if gt_types[i] in [1, 3])
        ctrl_carriers = sum(1 for i in control_indices if gt_types[i] in [1, 3])

        hap_p_value, hap_odds_ratio = calculate_fisher_exact(
            pos_carriers, pos_total, neg_carriers, neg_total
        )
        pat_p_value, pat_odds_ratio = calculate_fisher_exact(
            pat_carriers, pat_total, ctrl_carriers, ctrl_total
        )

        varinfo = (
            variant.INFO.get("SVLEN", "NA")
            if args.type == "sv"
            else get_snv_varinfo(variant)
        )
        seen_in_gwas = variant.POS in gwas_tested

        print(
            f"{variant.CHROM}\t{variant.POS}\t{varinfo}\t{pos_carriers}\t{neg_carriers}\t"
            f"{third_carriers}\t{pat_carriers}\t{ctrl_carriers}\t"
            f"{hap_p_value:.2e}\t{hap_odds_ratio:.2f}\t{pat_p_value:.2e}\t{pat_odds_ratio:.2f}\t"
            f"{seen_in_gwas}"
        )

    print(f"Kept {variants_kept} variants", file=sys.stderr)


def get_snv_varinfo(variant):
    change = variant.REF + ">" + ",".join(str(a) for a in variant.ALT)
    csq = variant.INFO.get("CSQ")
    effect = ",".join(set([ann.split("|")[1] for ann in csq.split(",")]))
    gene = ",".join(
        set([ann.split("|")[3] for ann in csq.split(",") if ann.split("|")[3]])
    )
    return change + " " + effect + " " + gene


def define_groups(args):
    # Load sample information
    sample_info = pd.read_csv(args.sampleinfo, sep="\t")

    # Define positive group
    positive_group = sample_info[
        (sample_info["haplotype"] == "major")
        & (sample_info["cohort"] == "aFTLD-U")
        & (sample_info["collection"] == "normal")
        & (sample_info["inclusion"] != "no")
    ]["individual"].tolist()

    # Define negative group
    negative_group = sample_info[
        (sample_info["haplotype"] == "none")
        & (sample_info["cohort"] == "in-house non-aFTLD-U")
        & (sample_info["collection"] == "normal")
        & (sample_info["inclusion"] != "no")
    ]["individual"].tolist()

    # Define third group - in-house non-aFTLD-U major haplotype carriers
    con_hapA_group = sample_info[
        (sample_info["haplotype"] == "major")
        & (sample_info["cohort"] == "in-house non-aFTLD-U")
        & (sample_info["collection"] == "normal")
        & (sample_info["inclusion"] != "no")
    ]["individual"].tolist()

    # Define patient group:
    patient_group = sample_info[
        (sample_info["cohort"] == "aFTLD-U")
        & (sample_info["collection"] == "normal")
        & (sample_info["inclusion"] != "no")
    ]["individual"].tolist()
    control_group = sample_info[
        (sample_info["cohort"] == "in-house non-aFTLD-U")
        & (sample_info["collection"] == "normal")
        & (sample_info["inclusion"] != "no")
    ]["individual"].tolist()

    for group in [positive_group, negative_group, con_hapA_group, patient_group, control_group]:
        print(f"Identified {len(group)} samples in the group for {group}", file=sys.stderr)


    # Open VCF file
    vcf = VCF(args.vcf)

    # Get indices of all groups in the VCF
    samples = vcf.samples
    positive_indices = [samples.index(s) for s in positive_group if s in samples]
    negative_indices = [samples.index(s) for s in negative_group if s in samples]
    con_hapA_indices = [samples.index(s) for s in con_hapA_group if s in samples]
    patient_indices = [samples.index(s) for s in patient_group if s in samples]
    control_indices = [samples.index(s) for s in control_group if s in samples]

    # Count missing samples
    missing_pos = len(positive_group) - len(positive_indices)
    missing_neg = len(negative_group) - len(negative_indices)
    missing_con_hapA = len(con_hapA_group) - len(con_hapA_indices)
    missing_patients = len(patient_group) - len(patient_indices)
    missing_controls = len(control_group) - len(control_indices)

    if missing_pos > 0:
        print(
            f"Warning: {missing_pos} samples from positive group not found in VCF",
            file=sys.stderr,
        )
    if missing_neg > 0:
        print(
            f"Warning: {missing_neg} samples from negative group not found in VCF",
            file=sys.stderr,
        )
    if missing_con_hapA > 0:
        print(
            f"Warning: {missing_con_hapA} samples from con hapA group not found in VCF",
            file=sys.stderr,
        )
    if missing_patients > 0:
        print(
            f"Warning: {missing_patients} samples from patient group not found in VCF",
            file=sys.stderr,
        )
    if missing_controls > 0:
        print(
            f"Warning: {missing_controls} samples from control group not found in VCF",
            file=sys.stderr,
        )
    return (
        positive_indices,
        negative_indices,
        con_hapA_indices,
        patient_indices,
        control_indices,
    )


def parse_gwas(gwas_file):
    gwas = pd.read_csv(gwas_file, sep="\t", usecols=["CHROM", "POS"])
    return gwas.loc[
        (gwas["CHROM"] == 15) & (gwas["POS"].between(34362469, 534862469)), "POS"
    ].tolist()


def parse_arguments():
    parser = ArgumentParser(
        description="Filter VCF for variants shared by haplotype carriers but not controls"
    )
    parser.add_argument(
        "-i", "--vcf", required=True, help="Input VCF file (can be gzipped)"
    )
    parser.add_argument(
        "--sampleinfo", required=True, help="TSV file with sample information"
    )
    parser.add_argument("--gwas", help="GWAS tsv file with tested variants")
    parser.add_argument(
        "--type",
        choices=["sv", "snv"],
        default="sv",
        help="Type of variant as input (default: sv)",
    )
    return parser.parse_args()


def calculate_fisher_exact(carriers1, total1, carriers2, total2):
    """
    Calculate Fisher's exact test for comparing two groups.

    Args:
        carriers1: Number of carriers in first group
        total1: Total number of samples in first group
        carriers2: Number of carriers in second group
        total2: Total number of samples in second group

    Returns:
        tuple: (p_value, odds_ratio)
    """
    # Calculate non-carriers
    non_carriers1 = total1 - carriers1
    non_carriers2 = total2 - carriers2

    # Handle edge cases
    if total1 == 0 or total2 == 0:
        return np.nan, np.nan

    # Create contingency table
    contingency_table = [[carriers1, non_carriers1], [carriers2, non_carriers2]]

    # Calculate Fisher's exact test
    odds_ratio, p_value = stats.fisher_exact(contingency_table)

    # Handle potential numerical issues with odds ratio
    if np.isinf(odds_ratio):
        # If one cell is 0, add a pseudocount of 0.5 to all cells
        if carriers1 == 0 or carriers2 == 0 or non_carriers1 == 0 or non_carriers2 == 0:
            adjusted_table = [
                [carriers1 + 0.5, non_carriers1 + 0.5],
                [carriers2 + 0.5, non_carriers2 + 0.5],
            ]
            odds_ratio, _ = stats.fisher_exact(adjusted_table)

    return p_value, odds_ratio


if __name__ == "__main__":
    main()
