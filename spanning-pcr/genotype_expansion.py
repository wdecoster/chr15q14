import os
import logging
import random
import extract_from_bam as extract


class Genotype(object):
    def __init__(
        self,
        individual,
        num_reads,
        reads_initially,
        seq,
        smallest_length_used,
        filter="PASS",
        group = None,
    ):
        self.individual = individual
        self.num_reads = num_reads
        self.reads_initially = reads_initially
        self.seq = seq if seq else "-"
        self.ct_dimer_count = self.count_ct_dimer(seq) if seq else 0
        self.hexamer_count = seq.count("CCCTCT") if seq else 0
        self.length = len(seq) if seq else 0
        self.smallest_length_used = smallest_length_used if smallest_length_used else 0
        self.filter = filter
        self.group = group

    def __repr__(self):
        return f"{self.individual}\t{self.num_reads}\t{self.reads_initially}\t{self.length}\t{self.ct_dimer_count}\t{self.hexamer_count}\t{self.seq}\t{self.filter}"

    @staticmethod
    def write_header():
        return "Individual\tNum reads expanded\tNum reads initially\tLength\tCT dimer count\tCCCTCT hexamer count\tConsensus sequence\tFilter"

    @staticmethod
    def count_ct_dimer(seq):
        _seq = seq
        for motif in ["CCCTCT", "CCCCT", "CCTTT", "CCTT", "CCCT", "CTTT"]:
            _seq = _seq.replace(motif, "-")
        return _seq.count("CT")


def genotype_samples(df, args):
    import concurrent.futures
    import tqdm

    def next_sample(df):
        for sample in df["sample"].unique():
            yield df[df["sample"] == sample]

    num_samples = len(df["sample"].unique())
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        genotypes = list(
            tqdm.tqdm(
                executor.map(genotype_sample, next_sample(df), [args] * num_samples),
                total=num_samples,
                desc="Genotyping",
            )
        )
    return genotypes


def genotype_sample(df_sample, args, min_reads=100):
    reads_before = len(df_sample)
    name = df_sample["sample"].iloc[0]
    df_sample = df_sample[
        (df_sample["flag"] == extract.ReadFlag.PASS_QC)
        & df_sample["length"].between(args.minlength, args.maxlength)
    ]

    num_reads = len(df_sample)

    # Set filter flag based on coverage
    filter_flag = "PASS" if num_reads >= min_reads else "LOWCOV"
    logging.info(
        f"Sample {name}: {reads_before} reads initially, {num_reads} after length filter "
        f"[{args.minlength}-{args.maxlength}], filter={filter_flag}"
    )

    # Always attempt consensus with available reads
    if num_reads > 0:
        import pypoars

        num_seqs = min(num_reads, min_reads)
        if args.downsample == "random":
            selected_reads = df_sample.sample(
                n=num_seqs, random_state=random.randint(0, 2**32 - 1)
            )
        else:
            selected_reads = df_sample.sort_values(by="length", ascending=False).head(
                num_seqs
            )
        seqs = selected_reads["seq"].to_list()
        logging.info(
            f"  Running POA consensus on {num_seqs} sequence(s) for {name} "
            f"using {args.downsample} downsampling"
        )
        consensus = pypoars.poa_consensus( # pyright: ignore[reportAttributeAccessIssue]
            seqs
        )
        return Genotype(
            name,
            num_reads,
            reads_before,
            consensus,
            selected_reads["length"].min(),
            filter=filter_flag,
        )
    else:
        logging.info(f"  No reads for {name}, skipping consensus")
        return Genotype(name, num_reads, reads_before, None, None, filter=filter_flag)