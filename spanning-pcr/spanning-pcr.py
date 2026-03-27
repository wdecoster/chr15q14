from argparse import ArgumentParser
import logging
import sys
import plotting
import genotype_expansion as gt
import extract_from_bam as extract


def main():
    args = get_args()
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="[%(levelname)s] %(message)s",
        stream=sys.stderr,
    )

    # Check if we should load from pickle
    if args.load_pickle:
        logging.info(f"Loading data from pickle: {args.load_pickle}")
        import pickle
        with open(args.load_pickle, 'rb') as f:
            data = pickle.load(f)
            df = data['df']
            genotypes = data['genotypes']
        logging.info(f"Loaded {len(df)} records and {len(genotypes)} genotypes from pickle")
    else:
        # Process data from CRAMs
        logging.info(f"Processing {len(args.crams)} CRAM file(s)")
        df = extract.get_data(args)
        logging.info(f"Loaded {len(df)} total records across {df['sample'].nunique()} sample(s)")
        genotypes = gt.genotype_samples(df, args)

        if args.pickle:
            logging.info(f"Saving processed data to pickle: {args.pickle}")
            import pickle
            with open(args.pickle, 'wb') as f:
                pickle.dump({'df': df, 'genotypes': genotypes}, f)

    if args.tsv:
        df.to_csv(args.tsv, sep="\t", index=False)
        logging.info(f"Wrote read DataFrame to {args.tsv}")

    # Display genotypes
    print(gt.Genotype.write_header())
    for genotype in sorted(genotypes, key=lambda x: x.ct_dimer_count):
        print(genotype)

    # Generate plots if requested
    if not args.noplot:
        logging.info(f"Writing plots to {args.output}")
        df_filtered = df[df["length"].between(args.minlength, args.maxlength)]
        with open(args.output, "w") as out:
            out.write(plotting.plot_genotypes(genotypes).to_html())
            out.write(plotting.plot_violins(df_filtered, genotypes).to_html())
            # out.write(plotting.ridges_plot(df_filtered).to_html())
            # out.write(plotting.scatter_plot(df_filtered, motif="CCCTCT count", full=args.full).to_html())
            out.write(plotting.scatter_plot(df_filtered, motif="CT count", full=args.full).to_html())
            out.write(plotting.scatter_motifs(df_filtered).to_html())
        logging.info("Done")

    if args.truth:
        import pandas as pd
        truth = pd.read_csv(args.truth, sep="\t")
        truth_dict = {row["name"]: row["sequence"] for _, row in truth.iterrows()}
        cram_lengths, truth_lengths, cram_ct, truth_ct, labels = [], [], [], [], []
        for g in genotypes:
            if g.individual in truth_dict:
                seq = truth_dict[g.individual]
                cram_lengths.append(g.length)
                truth_lengths.append(len(seq))
                cram_ct.append(g.ct_dimer_count)
                truth_ct.append(gt.Genotype.count_ct_dimer(seq))
                labels.append(g.individual)
        if labels:
            plotting.plot_truth_correlation(labels, cram_lengths, truth_lengths, cram_ct, truth_ct)
        else:
            logging.warning("No common individuals between CRAMs and truth file; skipping correlation plots.")


def get_args():
    parser = ArgumentParser()
    parser.add_argument("--crams", nargs="+", help="Input CRAM file(s)")
    parser.add_argument(
        "--minlength", type=int, default=50, help="Minimum fragment length"
    )
    parser.add_argument(
        "--maxlength", type=int, default=1200, help="Maximum fragment length"
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="Use entire read sequence, not just non-reference bases",
    )
    parser.add_argument(
        "--skip-softclip",
        action="store_true",
        help="When extracting non-ref bases, skip softclipped regions (consider only insertions)",
    )
    parser.add_argument("-o", "--output", default="read_lengths.html", help="Output file")
    parser.add_argument(
        "--threads", type=int, default=4, help="Number of threads to use"
    )
    parser.add_argument("--noplot", action="store_true", help="Don't make plots")
    parser.add_argument(
        "--pickle", 
        help="Path to save processed data as pickle for faster reuse"
    )
    parser.add_argument(
        "--load-pickle", 
        help="Path to load previously processed data in pickle format"
    )
    parser.add_argument(
        "--region",
        default="chr15:34419288-34419527",
        help="Genomic region to fetch reads from (format: chrom:start-end)",
    )
    parser.add_argument(
        "--downsample",
        choices=["longest", "random"],
        default="longest",
        help="Strategy to select up to 100 reads for consensus",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging to stderr",
    )
    parser.add_argument(
        "--truth",
        help="TSV file with columns 'individual' and 'sequence' for truth expansion sequences",
    )
    parser.add_argument(
        "--tsv",
        help="Path to dump the entire per-read DataFrame as a TSV file, for debugging",
    )
    parser.add_argument("--groups", help="Path to TSV file with columns 'sample' and 'group' to assign samples to groups for plotting")
    args = parser.parse_args()
    # verify input arguments
    if args.region != "chr15:34419288-34419527" and not args.skip_softclip:
        sys.exit(
            "Custom region specified without --skip-softclip; which is currently required."
        )
    if (args.load_pickle and args.crams) or (not args.load_pickle and not args.crams):
        sys.exit("Error: either --crams or --load-pickle must be specified, but not both")

    return args


if __name__ == "__main__":
    main()
