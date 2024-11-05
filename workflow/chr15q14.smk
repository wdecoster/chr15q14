import os
import glob
import pandas as pd
from os.path import basename

outdir = "/results/rr/study/hg38s/study252-P200_analysis/workflow_results/chr15q14/"
ref = "/home/wdecoster/GRCh38_recommended/GRCh38.fa"

# load the dataframe from the excel sheet
crams = pd.read_excel("/home/wdecoster/chr15q14/full_cohort_for_paper.xlsx")
# remove individuals that have "inclusion" set to "no"
# currently, also remove samples that are not_brain (corresponding to individuals for which multiple tissues were sequenced or no brain available) are dropped
# but this still has to be fixed later, to get those samples back in for specific analyses
crams = crams.loc[
    (crams["inclusion"] == "yes")
    & (~crams["collection"].isin(["not_brain", "other_brain"]))
]

local_cram_path = "/home/wdecoster/p200/study169-p200/good_crams/"
crams["path-to-cram"] = crams["path-to-cram"].astype(str)
crams["path-to-cram"] = crams.apply(
    lambda x: (
        x["path-to-cram"]
        if x["path-to-cram"].startswith("https://")
        else local_cram_path + x["individual"] + ".cram"
    ),
    axis=1,
)

# check that all paths that do not start with "https://" exist
for path in crams["path-to-cram"]:
    if not path.startswith("https://"):
        assert os.path.exists(path), f"Path does not exist: {path}"

# dump the dataframe to a tsv file, as it is used in the workflow
crams["name"] = crams["individual"]
crams.to_csv("/home/wdecoster/chr15q14/full_cohort_for_paper.tsv", sep="\t", index=False)

def get_duplicates():
    other_tissues = crams.loc[crams["collection"] == "other_tissue_duplicate", "individual"].tolist()
    fcxs = crams.loc[crams["collection"] == "normal", "individual"].tolist()
    samples = []
    for other in other_tissues:
        stripped = other.replace("_CELL_LINE", "").replace("_CER", "").replace("_OCX", "").replace("_tissue_CAU", "")
        if stripped in fcxs:
            samples.append(stripped)
        else:
            print(f"Sample {stripped} not found in the cohort")
    return samples + other_tissues

duplicates = get_duplicates()

def get_representative_cohort():
    """return all individuals from the normal collection, 
    but exclude those controls with the major haplotype that have been intentionally included, 
    as listed in privacy.intentionally_included_controls"""
    try:
        from privacy import intentionally_included_controls
        from privacy import asian_case
        return crams.loc[
            (crams["collection"] == "normal")
            & (~crams["individual"].isin(intentionally_included_controls+ asian_case)),
            "individual",
        ].tolist()
    except ImportError:
        return crams.loc[crams["collection"] == "normal", "individual"].tolist()

representative_cohort = get_representative_cohort()

def get_cram(wildcards):
    return crams.loc[crams["individual"] == wildcards.id, "path-to-cram"].values[0]

def fix_names_relatives(wildcards):
    """
    This very specific function prepares the names of the relatives for the axis labels in the aSTRonaut plot.
    This requires a dictionary with the names of the relatives, which is imported from a separate file.
    This separate file is not included in the repository, as it contains personal information (identifiers).
    As a fallback, the function will return the original identifiers.
    """
    try:
        from privacy import name_changes_relatives
        relatives = crams.loc[crams["collection"] == "relatives", "individual"].tolist()
        return ','.join([name_changes_relatives[r] for r in relatives])
    except ImportError:
        return ','.join(crams.loc[crams["collection"] == "relatives", "individual"].tolist())

def fix_names_duplicates(wildcards):
    """
    This very specific function prepares the names of the individuals for which multiple samples were sequenced for the axis labels in the aSTRonaut plot.
    This requires a dictionary with the names, which is imported from a separate file.
    The dictionary contains the individual names as keys and the new names as values.
    This is replaced in the identifiers, to preserve the sample type information.
    This separate file is not included in the repository, as it contains personal information (identifiers).
    As a fallback, the function will return the original identifiers.
    """
    try:
        from privacy import name_changes_duplicates
        out = []
        for d in duplicates:
            for k,v in name_changes_duplicates.items():
                if k == d:
                    out.append(f"{v}_FCX")
                elif k in d:
                    out.append(d.replace(k, v).replace("_CELL_LINE", "_LCL").replace("_tissue_CAU", "_CAU"))
        return ','.join(out)
    except ImportError:
        return ','.join(crams.loc[crams["collection"] == "other_tissue_duplicate", "individual"].tolist())


coords = {
    "golga8a_unphased": ["chr15:34419425-34419451", "--unphased"],
    "golga8a_phased": ["chr15:34419425-34419451", ""],
    "inbetween": ["chr15:34480576-34480608", "--unphased"],
    "golga8b": ["chr15:34565656-34565682", "--unphased"],
    "mga": ["chr15:41656320-41656381", ""],
    "linc02177": ["chr16:9393980-9394893", ""],
    "tata_golga8a": ["chr15:34454183-34454242", "--unphased"],
    "golga8a_utr": ["chr15:34437426-34437720", "--unphased"],
}

def get_coords(wildcards):
    try:
        return coords[wildcards.target][0]
    except KeyError:
        raise ValueError(f"Unknown target: {wildcards.target}")

def get_method(wildcards):
    try:
        return coords[wildcards.target][1]
    except KeyError:
        raise ValueError(f"Unknown target: {wildcards.target}")



targets = ["golga8a_phased", "golga8a_unphased"]  # select loci from the keys in coords dictionary
assert all([t in coords for t in targets]), "Not all targets are present in the coords dictionary"

rule all:
    input:
        strdust=expand(
            os.path.join(outdir, "strdust", "{target}/{id}.vcf.gz"),
            id=crams["individual"],
            target=targets,
        ),
        somatic_astronaut=expand(
            os.path.join(outdir, "plots/{target}/somatic_astronaut/{id}.html"),
            id=crams["individual"],
            target=targets,
        ),
        length=expand(
            os.path.join(outdir, "plots", "{target}/length_plot.html"),
            target=targets,
        ),
        lengthdelT=expand(os.path.join(outdir, "plots", "{target}/length_plot_delT.html"), target=targets),
        length625=expand(os.path.join(outdir, "plots", "{target}/length_plot_625.html"), target=targets),
        length_strip = expand(
            os.path.join(outdir, "plots", "{target}/length_plot_violin.html"),
            target=targets,
        ),
        astronaut_all=expand(os.path.join(outdir, "plots/{target}/aSTRonaut_all.html"), target=targets),
        astronaut_delT=expand(os.path.join(outdir, "plots/{target}/aSTRonaut_delT.html"), target=targets),
        astronaut_625=expand(os.path.join(outdir, "plots/{target}/aSTRonaut_625.html"), target=targets),
        astronaut_relatives=expand(os.path.join(outdir, "plots/{target}/aSTRonaut_relatives.html"), target=targets),
        astronaut_multiple_samples=expand(os.path.join(outdir, "plots/{target}/aSTRonaut_multiple_samples.html"), target=targets),
        kmer_plot=expand(os.path.join(outdir, "plots/{target}/kmer_plot.html"), target=targets),
        ct_vs_length=expand(os.path.join(outdir, "plots/{target}/ct_vs_length.html"), target=targets),
        combined_inquistr=os.path.join(outdir, "inquistr/representative_cohort.inq"),
        ct_stretch=expand(os.path.join(outdir, "analysis_overview-ct-stretch_{target}.tsv"), target=targets),
        corr_with_age = expand(os.path.join(outdir, "plots/{target}/correlations-with-age.html"), target=targets),
        corr_with_age_only_patients = expand(os.path.join(outdir, "plots/{target}/correlations-with-age_pat_only.html"), target=targets),
        copy_number_plot=os.path.join(outdir, "plots/copy_number.html"),




rule strdust:
    output:
        os.path.join(outdir, "strdust", "{target}/{id}.vcf.gz"),
    log:
        os.path.join(outdir, "logs/strdust", "{target}/{id}.log"),
    params:
        cram=get_cram, # using this as a param to avoid checking for existence of the cram file, as remote files are not considered present, and presence of local files is already checked
        ref=ref,
        locus=get_coords,
        binary="/home/wdecoster/repositories/STRdust/target/release/STRdust",
        method=get_method,
        support=2,
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/tabix.yml")
    shell:
        """RUST_LOG=debug {params.binary} \
        -r {params.locus} \
        {params.method} \
        --support {params.support} \
        --find-outliers \
        --somatic \
        {params.ref} {params.cram} 2> {log} | bgzip > {output} 2>> {log}"""


rule somatic_astronaut:
    input:
        os.path.join(outdir, "strdust/{target}/{id}.vcf.gz"),
    output:
        os.path.join(outdir, "plots/{target}/somatic_astronaut/{id}.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/somatic_astronaut_{id}.log"),
    params:
        script="/home/wdecoster/pathSTR-1000G/scripts/aSTRonaut.py",
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """
        python {params.script} \
        --motifs CT,CCTT,CTTT,CCCT,CCCTCT,CCCCT,CCCCCC \
        --somatic -m 100 --publication \
        {input} -o {output} 2> {log}
        """


rule length_plot:
    input:
        expand(
            os.path.join(outdir, "strdust", "{{target}}/{id}.vcf.gz"),
            id=crams.loc[crams["collection"] == "normal", "individual"], 
        ), # this corresponds to our own cohort
    output:
        os.path.join(outdir, "plots", "{target}/length_plot.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/length_plot.log"),
    params:
        script=os.path.join(
            os.path.dirname(workflow.basedir),
            "scripts/plot_repeat_lengths_from_vcf.py",
        ),
        minlen=100,
        sampleinfo="/home/wdecoster/cohorts/Individuals.xlsx",
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """python {params.script} -i {input} -o {output} --title "Repeat length per read" --mark_hexamers --minlen {params.minlen} --sampleinfo {params.sampleinfo} 2> {log}"""


rule length_plot_strip:
    input:
        expand(
            os.path.join(outdir, "strdust/{{target}}/{id}.vcf.gz"),
            id=crams.loc[crams["collection"] == "normal", "individual"], # this corresponds to our own cohort
        ),
    output:
        os.path.join(outdir, "plots", "{target}/length_plot_violin.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/length_plot_violin.log"),
    params:
        script=os.path.join(
            os.path.dirname(workflow.basedir),
            "scripts/repeat_length_violin.py",
        ),
        sample_info = "/home/wdecoster/fus-analysis/full_cohort_for_paper.tsv",
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        "python {params.script} -o {output} -g {params.sample_info} {input} 2> {log}"


rule length_plot_delT:
    input:
        expand(
            os.path.join(outdir, "strdust", "{{target}}/{id}.vcf.gz"),
            id=crams.loc[(crams["collection"] == "normal") & (crams["haplotype"] == 'major'), "individual"], # this corresponds to our own cohort, major haplotype carriers
        ),
    output:
        os.path.join(outdir, "plots", "{target}/length_plot_delT.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/length_plot_delT.log"),
    params:
        script=os.path.join(
            os.path.dirname(workflow.basedir),
            "scripts/plot_repeat_lengths_from_vcf.py",
        ),
        minlen=0,
        sampleinfo="/home/wdecoster/cohorts/Individuals.xlsx",
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """python {params.script} -i {input} --title "Repeat length per read for haplotype A" --mark_hexamers -o {output} --minlen {params.minlen} --sampleinfo {params.sampleinfo} 2> {log}"""


rule length_plot_625:
    input:
        expand(
            os.path.join(outdir, "strdust", "{{target}}/{id}.vcf.gz"),
            id=crams.loc[(crams["collection"] == "normal") & (crams["haplotype"] == 'minor'), "individual"], # this corresponds to our own cohort, minor haplotype carriers
        ),
    output:
        os.path.join(outdir, "plots", "{target}/length_plot_625.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/length_plot_625.log"),
    params:
        script=os.path.join(
            os.path.dirname(workflow.basedir),
            "scripts/plot_repeat_lengths_from_vcf.py",
        ),
        minlen=0,
        sampleinfo="/home/wdecoster/cohorts/Individuals.xlsx",
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """python {params.script} -i {input} -o {output} --title "Repeat length per read for haplotype B" --mark_hexamers --minlen {params.minlen} --sampleinfo {params.sampleinfo} 2> {log}"""


# rule length_plot_duplicates:
#     input:
#         expand(
#             os.path.join(outdir, "strdust/{{target}}/{id}.vcf.gz"),
#             id=duplicates,
#         ),
#     output:
#         os.path.join(outdir, "plots", "{target}/length_plot_duplicates.html"),
#     log:
#         os.path.join(outdir, "logs/workflows/{target}/length_plot_duplicates.log"),
#     params:
#         script=os.path.join(
#             os.path.dirname(workflow.basedir),
#             "analysis/plot_repeat_lengths_from_vcf.py",
#         ),
#         minlen=0,
#         sampleinfo="/home/wdecoster/cohorts/Individuals.xlsx",
#     conda:
#         os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
#     shell:
#         "python {params.script} -i {input} -o {output} --minlen {params.minlen} --alphabetically 2> {log}"


rule kmer_plot:
    input:
        expand(
            os.path.join(outdir, "strdust", "{{target}}/{id}.vcf.gz"),
            id=crams.loc[crams["collection"] == "normal", "individual"],
        ),
    output:
        html=os.path.join(outdir, "plots", "{target}/kmer_plot.html"),
        counts=os.path.join(outdir, "plots", "{target}/kmer_counts.tsv"),
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    params:
        script=os.path.join(
            os.path.dirname(workflow.basedir), "scripts/kmer-frequencies-vcf.py"
        ),
        kmer=12,
        sampleinfo="/home/wdecoster/p200/fus-analysis/full_cohort_for_paper.tsv",
    log:
        os.path.join(outdir, "logs/{target}/kmer_plot.log"),
    shell:
        """
        python {params.script} \
        -k {params.kmer} \
        --output {output.html} \
        --counts {output.counts} \
        --sampleinfo {params.sampleinfo} \
        {input} 2> {log}"""

# rule kmer_plot_duplicates:
#     input:
#         expand(
#             os.path.join(outdir, "strdust/{{target}}/{id}.vcf.gz"),
#             id=duplicates,
#         ),
#     output:
#         html=os.path.join(outdir, "plots", "{target}/kmer_plot_duplicates.html"),
#         counts=os.path.join(outdir, "plots", "{target}/kmer_counts_duplicates.tsv"),
#     conda:
#         os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
#     params:
#         script=os.path.join(
#             os.path.dirname(workflow.basedir), "analysis/kmer-frequencies-vcf.py"
#         ),
#         kmer=12,
#         sampleinfo="/home/wdecoster/cohorts/Individuals.xlsx",
#     log:
#         os.path.join(outdir, "logs/{target}/kmer_plot_duplicates.log"),
#     shell:
#         """
#         python {params.script} \
#         -k {params.kmer} \
#         --output {output.html} \
#         --counts {output.counts} \
#         --sampleinfo {params.sampleinfo} \
#         {input} 2> {log}"""


rule astronaut_all:
    input:
        expand(
            os.path.join(outdir, "strdust/{{target}}/{id}.vcf.gz"),
            id=crams.loc[crams["collection"] == "normal", "individual"], # this corresponds to our own cohort
        ),
    output:
        os.path.join(outdir, "plots/{target}/aSTRonaut_all.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/aSTRonaut_all.log"),
    params:
        script="/home/wdecoster/pathSTR-1000G/scripts/aSTRonaut.py",
        sample_info="/home/wdecoster/chr15q14/full_cohort_for_paper.tsv",
        dotsize=8
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """
        python {params.script} \
        --motifs CT,CCTT,CTTT,CCCT,CCCTCT,CCCCT,CCCCCC \
        -m 100 \
        -o {output} \
        --size {params.dotsize} \
        --hide-labels --longest_only --publication \
        --title "Repeat composition sequence" \
        --sampleinfo {params.sample_info} \
        {input} 2> {log}
        """


rule astronaut_delT:
    input:
        expand(
            os.path.join(outdir, "strdust/{{target}}/{id}.vcf.gz"),
            id=crams.loc[(crams["collection"] == "normal") & (crams["haplotype"] == 'major'), "individual"], # this corresponds to our own cohort, major haplotype carriers
        ),
    output:
        os.path.join(outdir, "plots/{target}/aSTRonaut_delT.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/aSTRonaut_delT.log"),
    params:
        script="/home/wdecoster/pathSTR-1000G/scripts/aSTRonaut.py",
        sample_info="/home/wdecoster/chr15q14/full_cohort_for_paper.tsv",
        dotsize = 8,
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """
        python {params.script} \
        --motifs CT,CCTT,CTTT,CCCT,CCCTCT,CCCCT,CCCCCC \
        -m 100 \
        -o {output} \
        --size {params.dotsize} \
        --hide-labels --longest_only --publication \
        --title "Repeat composition sequence in carriers of haplotype A" \
        --sampleinfo {params.sample_info} \
        {input} 2> {log}
        """


rule astronaut_625:
    input:
        expand(
            os.path.join(outdir, "strdust/{{target}}/{id}.vcf.gz"),
            id=crams.loc[(crams["collection"] == "normal") & (crams["haplotype"] == 'minor'), "individual"], # this corresponds to our own cohort, minor haplotype carriers
        ),
    output:
        os.path.join(outdir, "plots/{target}/aSTRonaut_625.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/aSTRonaut_625.log"),
    params:
        script="/home/wdecoster/pathSTR-1000G/scripts/aSTRonaut.py",
        sample_info="/home/wdecoster/chr15q14/full_cohort_for_paper.tsv",
        dotsize = 8,
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """
        python {params.script} \
        --motifs CT,CCTT,CTTT,CCCT,CCCTCT,CCCCT,CCCCCC \
        -m 100 \
        -o {output} \
        --size {params.dotsize} \
        --hide-labels --longest_only --publication \
        --title "Repeat composition sequence in carriers of haplotype B" \
        --sampleinfo {params.sample_info} \
        {input} 2> {log}
        """

rule astronaut_relatives:
    input:
        expand(
            os.path.join(outdir, "strdust/{{target}}/{id}.vcf.gz"),
            id=crams.loc[crams["collection"] == "relatives", "individual"], # this corresponds to the individuals with relatives in the cohort
        ),
    output:
        os.path.join(outdir, "plots/{target}/aSTRonaut_relatives.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/aSTRonaut_relatives.log"),
    params:
        script="/home/wdecoster/pathSTR-1000G/scripts/aSTRonaut.py",
        sample_info="/home/wdecoster/chr15q14/full_cohort_for_paper.tsv",
        relative_names=fix_names_relatives,
        dotsize = 8,
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """
        python {params.script} \
        --motifs CT,CCTT,CTTT,CCCT,CCCTCT,CCCCT,CCCCCC \
        --names {params.relative_names} \
        --label_size 20 \
        --alphabetic \
        --hide_allele_label \
        --size {params.dotsize} \
        {input} -m 100 -o {output} --longest_only --publication --title "Repeat composition sequence in relatives" --sampleinfo {params.sample_info} 2> {log}
        """

rule astronaut_multiple_samples:
    input:
        expand(
            os.path.join(outdir, "strdust/{{target}}/{id}.vcf.gz"),
            id=duplicates, # this corresponds to the individuals with multiple samples in the cohort
        ),
    output:
        os.path.join(outdir, "plots/{target}/aSTRonaut_multiple_samples.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/aSTRonaut_multiple_samples.log"),
    params:
        script="/home/wdecoster/pathSTR-1000G/scripts/aSTRonaut.py",
        sample_info="/home/wdecoster/chr15q14/full_cohort_for_paper.tsv",
        duplicate_names=fix_names_duplicates
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """
        python {params.script} \
        --motifs CT,CCTT,CTTT,CCCT,CCCTCT,CCCCT,CCCCCC \
        --names {params.duplicate_names} \
        --label_size 20 \
        --alphabetic \
        --hide_allele_label \
        {input} -m 100 -o {output} --longest_only --publication --title "Repeat composition sequence in individuals with multiple samples" --sampleinfo {params.sample_info} 2> {log}
        """

rule ct_vs_length:
    input:
        expand(
            os.path.join(outdir, "strdust/{{target}}/{id}.vcf.gz"),
            id=crams.loc[crams["collection"].isin(["normal", "1000G", "owen", "kristel"]), "individual"], # this corresponds to the full cohort
        ),
    output:
        plot = os.path.join(outdir, "plots/{target}/ct_vs_length.html"),
        overview = os.path.join(outdir, "analysis_overview_{target}.tsv"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/ct_vs_length.log"),
    params:
        script=os.path.join(
            os.path.dirname(workflow.basedir),
            "scripts/plot_ct_vs_repeat_length.py",
        ),
        sample_info = "/home/wdecoster/chr15q14/full_cohort_for_paper.tsv",
        xline = 0.8,
        yline = 400,
        arrows = "rr_UCL2783,rr_CW05_004"
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """
        python {params.script} \
        --output {output.plot} \
        --overview {output.overview} \
        --groups {params.sample_info} \
        --haplotypes \
        --xline {params.xline} \
        --yline {params.yline} \
        --arrow {params.arrows} \
        {input} 2> {log} 
        """

rule ct_stretch:
    input:
        os.path.join(outdir, "analysis_overview_{target}.tsv"),
    output:
        os.path.join(outdir, "analysis_overview-ct-stretch_{target}.tsv"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/ct_stretch.log"),
    params:
        script=os.path.join(
            os.path.dirname(workflow.basedir),
            "scripts/ct-stretch.py",
        ),
    shell:
        """
        python {params.script} {input} > {output} 2> {log}
        """

rule correlate_with_age:
    input:
        os.path.join(outdir, "analysis_overview-ct-stretch_{target}.tsv"),
    output:
        os.path.join(outdir, "plots/{target}/correlations-with-age.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/ct_stretch_age.log"),
    params:
        sample_info = "/home/wdecoster/cohorts/Individuals.xlsx",
        script = os.path.join(
            os.path.dirname(workflow.basedir),
            "scripts/correlate_with_age.py",
        ),
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """
        python {params.script} {input} --sampleinfo {params.sample_info} > {output} 2> {log}
        """


rule correlate_with_age_only_patients:
    input:
        os.path.join(outdir, "analysis_overview-ct-stretch_{target}.tsv"),
    output:
        os.path.join(outdir, "plots/{target}/correlations-with-age_pat_only.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/ct_stretch_age.log"),
    params:
        sample_info = "/home/wdecoster/cohorts/Individuals.xlsx",
        script = os.path.join(
            os.path.dirname(workflow.basedir),
            "scripts/correlate_with_age.py",
        ),
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """
        python {params.script} {input} --pat_only --sampleinfo {params.sample_info} > {output} 2> {log}
        """




rule make_combined_inquistr_file:
    input:
        expand("/results/rr/study/hg38s/study252-P200_analysis/workflow_results/inquistr/inquistr-call/{id}.inq.gz", id=representative_cohort)
    output:
        os.path.join(outdir, "inquistr/representative_cohort.inq")
    log:
        os.path.join(outdir, "logs/workflows/combine_inquistr.log")
    params:
        inquistr = "/home/wdecoster/repositories/inquiSTR/target/release/inquiSTR"
    shell:
        "{params.inquistr} combine {input} > {output} 2> {log}"


rule mosdepth:
    output:
        os.path.join(outdir, "mosdepth/{id}.regions.bed.gz"),
    log:
        os.path.join(outdir, "logs/workflows/{id}-mosdepth.log"),
    params:
        cram=get_cram, # using this as a param to avoid checking for existence of the cram file, as remote files are not considered present, and presence of local files is already checked
        chrom="chr15",
        window=100,
        ref=ref,
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/mosdepth.yml")
    threads: 5
    shell:
        """
        output=$(dirname {output}) && \
        mosdepth --mapq 10 \
        --chrom {params.chrom} \
        --by {params.window} \
        -t {threads} \
        --no-per-base \
        --fast-mode \
        --fasta {params.ref} \
        ${{output}}/{wildcards.id} {params.cram} 2> {log}
        """


rule call_copy_number:
    input:
        os.path.join(outdir, "mosdepth/{id}.regions.bed.gz"),
    output:
        os.path.join(outdir, "mosdepth/copy_number_{id}.tsv"),
    log:
        os.path.join(outdir, "logs/workflows/calculate_copy_number_{id}.log"),
    params:
        call="chr15:34438297-34524132",
        norm="chr15:54,033,377-56,279,876",
        script=os.path.join(
            os.path.dirname(workflow.basedir), "scripts/calculate_copy_number.py"
        ),
    shell:
        "python {params.script} -i {input} -c {params.call} --control {input} -n {params.norm} > {output} 2> {log}"


rule cat_copy_number_files:
    input:
        expand(
            os.path.join(outdir, "mosdepth/copy_number_{id}.tsv"),
            id=crams.loc[crams["collection"].isin(["normal", "1000G", "owen", "kristel"]), "individual"], # this corresponds to the full cohort
        ),
    output:
        os.path.join(outdir, "mosdepth/copy_number.tsv"),
    log:
        os.path.join(outdir, "logs/workflows/cat_copy_number_files.log"),
    shell:
        "cat {input} | sort -k2,2n > {output} 2> {log}"


rule plot_copy_number:
    input:
        os.path.join(outdir, "mosdepth/copy_number.tsv"),
    output:
        os.path.join(outdir, "plots/copy_number.html"),
    log:
        os.path.join(outdir, "logs/workflows/plot_copy_number.log"),
    params:
        sampleinfo="/home/wdecoster/cohorts/Individuals.xlsx",
        script=os.path.join(
            os.path.dirname(workflow.basedir), "scripts/plot_copy_numbers.py"
        ),
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        "python {params.script} -i {input} -o {output} --sampleinfo {params.sampleinfo} 2> {log}"


# rule join_lengths_and_kmers:
#     input:
#         golga=os.path.join(outdir, "{taget}/length_summary.tsv"),
#         copy_number=os.path.join(outdir, "mosdepth/copy_number.tsv"),
#         counts=os.path.join(outdir, "plots", "{target}/kmer_counts.tsv"),
#     output:
#         os.path.join(outdir, "summary/{target}.tsv"),
#     log:
#         os.path.join(outdir, "logs/workflows/{target}/oin_lengths_and_kmers.log"),
#     params:
#         script=os.path.join(
#             os.path.dirname(workflow.basedir), "analysis/make_fus_overview.py"
#         ),
#     conda:
#         os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
#     shell:
#         "python {params.script} -l {input.golga} -c {input.copy_number} -k {input.counts} > {output} 2> {log}"


# rule make_overview:
#     input:
#         tata=os.path.join(outdir, "{tata}/length_summary.tsv"),
#         golga=os.path.join(outdir, "{golga8a}/length_summary.tsv"),
#         mga=os.path.join(outdir, "{mga}/length_summary.tsv"),
#         linc02177=os.path.join(outdir, "{linc02177}/length_summary.tsv"),
#         copy_number=os.path.join(outdir, "mosdepth/copy_number.tsv"),
#     output:
#         os.path.join(outdir, "overview_table.tsv"),
#     log:
#         os.path.join(outdir, "logs/workflows/make_overview.log"),
#     params:
#         sampleinfo="/home/wdecoster/cohorts/Individuals.xlsx",
#         script=os.path.join(
#             os.path.dirname(workflow.basedir), "scripts/make_overview.py"
#         ),
#     conda:
#         os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
#     shell:
#         """
#         python {params.script} \
#         --tata {input.tata} \
#         --golga {input.golga} \
#         --mga {input.mga} \
#         --linc {input.linc02177} \
#         --cn {input.copy_number} > {output} 2> {log}
#         """
