# chr15q14

This repository contains the snakemake workflow and scripts for the analysis of the GOLGA8A repeat expansion using long read sequencing.

The workflow, `workflow/chr15q14.smk` uses a tab separate file to describe the whole cohort. This file contains sample identifiers and can as such not be publicly shared.
The workflow uses conda environment (from yml files in `envs/`) to install the necessary software and run reproducible analysis.

## Analysis of spanning PCR data generated on flongle

The `scripts/analyze_flongle_data.py` script can be used to analze results from the spanning PCR assay to genotype shorter alleles of the expansion.
It uses Rust code (`pypoars`) to generate a consensus sequence of the PCR products of expansions. See below on installation and usage.
Further improvements can be made on how the rust code is packaged (installable), please reach out if the code below is not working as expected.

```bash
mamba create -n analyze_spanning pysam plotly pandas maturin tqdm fuzzysearch -y
mamba activate analyze_spanning
cd /path/to/repository/chr15q14/pypoars/
maturin develop --release
python ~/chr15q14/scripts/analyze_flongle_data.py --crams refmasked/*.cram --output flongle_analysis.html > genotypes.tsv
```

```text
USAGE
analyze_flongle_data.py [-h] [--crams CRAMS [CRAMS ...]] [--minlength MINLENGTH] [--maxlength MAXLENGTH] [--full] [--output OUTPUT] [--threads THREADS] [--noplot]

options:
  -h, --help            show this help message and exit
  --crams CRAMS         Input CRAM file(s) (minimum 1)
  --minlength MINLENGTH Minimum fragment length [default 50]
  --maxlength MAXLENGTH Maximum fragment length [default 1200]
  --full                Use entire read sequence, not just non-reference bases
  --output OUTPUT       Output file
  --threads THREADS     Number of threads to use
  --noplot              Don't make plots
```

Note that this script expects cram files to be generated based on an alignment to a reference genome in which the GOLGA8B gene is hardmasked, e.g. using bedtools and a bed file with the coordinates of the GOLGA8B gene (chr15 34525095 34583651)

```bash
bedtools maskfasta -fi ~/GRCh38.fa -bed mask.bed -fo GRCh38_with_GOLGA8B_masked.fa
```
