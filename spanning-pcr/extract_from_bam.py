import logging
import os
from typing import NamedTuple, Optional
from enum import IntEnum
import sys

class ReadFlag(IntEnum):
    PASS_QC = 0
    ACCIDENTAL_2D = 1
    NO_SEQ = 2


class ReadFeatures(NamedTuple):
    sample: str
    read_id: str
    strand: str
    length: int
    seq: Optional[str]
    ccctct_count: int
    ct_count: int
    flag: ReadFlag


def load_one_cram(cram, args):
    """Load reads from a single CRAM file. Runs in a worker process."""
    import pysam  # pyright: ignore[reportMissingImports]

    name = (
        os.path.basename(cram)
        .replace("masked_rm_map-sminimap2-", "")
        .replace("_hg38s.cram", "")
        .split(".")[0]
        .replace("_v7", "")
    )
    use_softclip = not args.skip_softclip
    chrom, coords = args.region.split(":")
    start, end = (int(x) for x in coords.split("-"))
    records = []
    with pysam.AlignmentFile(cram) as f:
        for r in f.fetch(chrom, start, end):
            records.append(parse_read(r, args.full, name, use_softclip))
    return records


def get_data(args):
    import concurrent.futures
    import pandas as pd
    import tqdm

    logging.info(
        f"Extracting reads from region {args.region} in {len(args.crams)} CRAM file(s)"
    )
    logging.info(f"  Using full read sequences: {args.full}")
    logging.info(
        f"  Using softclipped regions as non-ref bases: {not args.skip_softclip}"
    )

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {
            executor.submit(load_one_cram, cram, args): cram for cram in args.crams
        }
        records = []
        for future in tqdm.tqdm(
            concurrent.futures.as_completed(futures),
            total=len(futures),
            desc="Loading CRAMs",
        ):
            cram_records = future.result()
            logging.info(f"  {len(cram_records)} reads fetched from {futures[future]}")
            records.extend(cram_records)

    df = pd.DataFrame(
        records,
        columns=[
            "sample",
            "read_id",
            "strand",
            "length",
            "seq",
            "CCCTCT count",
            "CT count",
            "flag",
        ],
    )
    df["flag"] = df["flag"].astype("category")
    return df


def parse_read(read, full, name, use_softclip=True):
    """
    Extracts features of each read as a ReadFeatures named tuple.
    Filters out reads that are likely to be accidental 2D reads (ONT artefact)
     based on the presence of an SA tag indicating a chimeric alignment on the opposite strand.
    Or reads that have no sequence record.
    """
    if is_accidental_2d_read(read):
        return ReadFeatures(
            name,
            read.query_name,
            get_strand(read),
            0,
            read.query_sequence,  # include sequence for debugging purposes, even though it will be flagged as ACCIDENTAL_2D and not used
            0,
            0,
            ReadFlag.ACCIDENTAL_2D,
        )
    seq = (
        read.query_sequence if full else non_ref_bases(read, use_softclip=use_softclip)
    )
    if not seq:
        return ReadFeatures(
            name, read.query_name, get_strand(read), 0, None, 0, 0, ReadFlag.NO_SEQ
        )
    fragment_length = read.query_length if full else len(seq)
    return ReadFeatures(
        sample=name,
        read_id=read.query_name,
        strand=get_strand(read),
        length=fragment_length,
        seq=seq,
        ccctct_count=get_ccctct(seq),
        ct_count=count_ct_by_subtracting_motifs(seq),
        flag=ReadFlag.PASS_QC,
    )


def is_accidental_2d_read(read):
    """Check if a read alignment might be an accidental 2D read (ONT artefact).

    An accidental 2D read means that right after the template strand, the complement
    strand was also sequenced. The read will align in two pieces of similar length to
    the reference genome, with the second piece on the opposite strand.
    """
    read_strand = "-" if read.is_reverse else "+"

    if not read.has_tag("SA"):
        return False

    sa_tag = read.get_tag("SA")
    if not isinstance(sa_tag, str):
        logging.warning(f"Unexpected type of SA auxiliary tag: {type(sa_tag)}")
        return False

    sa_entries = [e for e in sa_tag.split(";") if e]

    if len(sa_entries) != 1:
        return False

    fields = sa_entries[0].split(",")
    if len(fields) < 4:
        logging.warning("Malformed SA tag entry: insufficient fields")
        return False

    sa_strand = fields[2]
    if read_strand == sa_strand:
        return False

    start = read.reference_start
    end = read.reference_end
    sa_start = int(fields[1])
    sa_end = sa_start + cigar_string_to_rlen(fields[3])

    if max(start, sa_start) < min(end, sa_end):
        logging.debug(
            f"Identified read {read.query_name} as accidental 2D read, discarding"
        )
        return True

    return False


def cigar_string_to_rlen(cigar_string):
    """Calculate reference length consumed by a CIGAR string."""
    import re

    # Operations that consume reference: M(0), D(2), N(3), =(7), X(8)
    ref_consuming = {"M", "D", "N", "=", "X"}
    rlen = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar_string):
        if op in ref_consuming:
            rlen += int(length)
    return rlen


def non_ref_bases(read, minlength=50, use_softclip=True):
    """
    This function slices out the bases from a read that do not match the reference genome, by parsing the CIGAR string.
    Only cigar operations longer than minlength are considered, and only insertions and softclips are returned.
    The script iterates over the cigar string, while moving cursors in both the read sequence and the reference.
    Softclips are only extracted if the reference position falls within the target interval.
    """
    if not read.query_sequence:
        return ""
    non_ref = ""
    read_position = 0
    reference_position = read.reference_start
    softclip_start = 34419380  # start of interval in which softclips are considered
    softclip_end = 34419494  # end of interval in which softclips are considered
    operations_to_consider = [1, 4] if use_softclip else [1]
    for operation, length in read.cigartuples:
        if operation in [0, 7, 8]:
            # operation 0 is match (M), 7 is match with sequence (=), 8 is alignment match (X)
            read_position += length
            reference_position += length
        elif operation == 3:
            # operation 3 is refskip (N), this shouldn't happen
            sys.stderr.write("Warning: unexpected refskip cigar operation in read\n")
            reference_position += length
        elif operation == 2:
            # operation 2 is deletion (D). This script does not care about repeat contractions.
            reference_position += length
        elif operation in [1, 4]:
            # operation 4 is softclip (S), operation 1 is insertion (I)
            # read_position must always advance for both, since both consume query sequence
            if operation in operations_to_consider and length >= minlength:
                # For softclips, only extract if reference position is within the target interval
                if operation == 4 and not (
                    softclip_start < reference_position < softclip_end
                ):
                    read_position += length
                    continue
                non_ref += read.query_sequence[read_position : read_position + length]
            read_position += length
    # attempt to trim off the reference sequences that may be caugth in the softclipped region
    # however, it is not likely that perfect matches are found for every read
    # therefore, using a fuzzy search allowing for a few mismatches (~5%) is used
    from fuzzysearch import find_near_matches  # pyright: ignore[reportMissingImports]

    right_seq = "GAGACGGAGTTTCTCTCTTGTTGCCCAGGCTGGAGTGCATGTTGCTGTGCACTTTGAGGGCAGGAACTG"
    matches = find_near_matches(right_seq, non_ref, max_l_dist=4)
    if len(matches) > 1:
        return None
    if len(matches) == 1:
        non_ref = non_ref[: matches[0].start]

    left_seq = "TTATCAGGGCCTCTCTTCGCAGGCAGTGGGGCCTCATCCACAACCCTGGAAAAGAACTGGAAAGCGTTGCTCAGCCAGGTACGGAGGGCAGGGCCATGTGGGACTCCCGTCTCCAGGCCCCCTCTCCCCAGCTCCCG"
    matches = find_near_matches(left_seq, non_ref, max_l_dist=7)
    if len(matches) > 1:
        return None
    if len(matches) == 1:
        non_ref = non_ref[matches[0].end :]
    return non_ref


def get_strand(read):
    return "-" if read.is_reverse else "+"


def get_ccctct(seq):
    return seq.count("CCCTCT") if seq else 0


def count_ct_by_subtracting_motifs(seq):
    if not seq:
        return 0
    for motif in ["CCCTCT", "CCCCT", "CCTTT", "CCTT", "CCCT", "CTTT"]:
        seq = seq.replace(motif, "-")
    return seq.count("CT")
