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



targets = ["golga8a_unphased"]  # select loci from the keys in coords dictionary
# note that working with phased data is problematic for unphased samples, particularly 1000G
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
            os.path.join(outdir, "plots", "{target}/somatic_length_plot.html"),
            target=targets,
        ),
        #lengthdelT=expand(os.path.join(outdir, "plots", "{target}/length_plot_delT.html"), target=targets),
        #length625=expand(os.path.join(outdir, "plots", "{target}/length_plot_625.html"), target=targets),
        length_strip = expand(
            os.path.join(outdir, "plots", "{target}/length_plot_violin.html"),
            target=targets,
        ),
        astronaut_all=expand(os.path.join(outdir, "plots/{target}/aSTRonaut_all.html"), target=targets),
        astronaut_hapA=expand(os.path.join(outdir, "plots/{target}/aSTRonaut_hapA.html"), target=targets),
        astronaut_hapB=expand(os.path.join(outdir, "plots/{target}/aSTRonaut_hapB.html"), target=targets),
        astronaut_relatives=expand(os.path.join(outdir, "plots/{target}/aSTRonaut_relatives.html"), target=targets),
        astronaut_multiple_samples=expand(os.path.join(outdir, "plots/{target}/aSTRonaut_multiple_samples.html"), target=targets),
        kmer_plot=expand(os.path.join(outdir, "plots/{target}/kmer_plot.html"), target=targets),
        ct_vs_length=expand(os.path.join(outdir, "plots/{target}/ct_vs_length.html"), target=targets),
        ct_dimer_strip=expand(os.path.join(outdir, "plots/{target}/ct_dimer_strip.html"), target=targets),
        combined_inquistr=os.path.join(outdir, "inquistr/representative_cohort.inq"),
        corr_with_age = expand(os.path.join(outdir, "plots/{target}/correlations-with-age.html"), target=targets),
        corr_with_age_only_patients = expand(os.path.join(outdir, "plots/{target}/correlations-with-age_pat_only.html"), target=targets),
        copy_number_plot=os.path.join(outdir, "plots/copy_number.html"),
        overview = expand(os.path.join(outdir, "analysis_overview_{target}.tsv"), target=targets),
        table_carriers = expand(os.path.join(outdir, "tables/haplotype_carriers_{target}.xlsx"), target=targets),
        somatic_variation = expand(os.path.join(outdir, "plots/{target}/somatic_variation_plot.html"), target=targets),
        sex_check = os.path.join(outdir, "sex_check.html"),
        precision_recall = expand(os.path.join(outdir, "plots/{target}/precision_recall.txt"), target=targets),


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
        --motifs CT,CCTT,CTTT,CCCT,CCCTCT,CCCCT,CCTTT,CCCCCC \
        --somatic -m 100 --publication \
        {input} -o {output} 2> {log}
        """


rule somatic_length_plot:
    input:
        expand(
            os.path.join(outdir, "strdust", "{{target}}/{id}.vcf.gz"),
            id=crams.loc[crams["collection"].isin(["normal", "1000G", "owen", "kristel", "mayo"]), "individual"],
        ), # this corresponds to the full cohort, without the relatives
    output:
        os.path.join(outdir, "plots", "{target}/somatic_length_plot.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/length_plot.log"),
    params:
        script=os.path.join(
            os.path.dirname(workflow.basedir),
            "scripts/plot_repeat_lengths_from_vcf.py",
        ),
        minlen=100,
        sampleinfo="/home/wdecoster/chr15q14/full_cohort_for_paper.tsv",
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """python {params.script} -i {input} -o {output} --title "Repeat length per read" --minlen {params.minlen} --sampleinfo {params.sampleinfo} 2> {log}"""


rule length_plot_strip:
    input:
        expand(
            os.path.join(outdir, "strdust/{{target}}/{id}.vcf.gz"),
            id=crams.loc[crams["collection"] == "normal", "individual"], 
        ), # this corresponds to our own cohort
    output:
        os.path.join(outdir, "plots", "{target}/length_plot_violin.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/length_plot_violin.log"),
    params:
        script=os.path.join(
            os.path.dirname(workflow.basedir),
            "scripts/repeat_length_violin.py",
        ),
        sample_info="/home/wdecoster/chr15q14/full_cohort_for_paper.tsv",
        show_line = 100,
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        "python {params.script} --line {params.show_line} -o {output} -g {params.sample_info} {input} 2> {log}"


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


rule somatic_variation_plot:
    input:
        expand(
            os.path.join(outdir, "strdust", "{{target}}/{id}.vcf.gz"),
            id=crams.loc[crams["collection"] == "normal", "individual"],
        ),
    output:
        os.path.join(outdir, "plots/{target}/somatic_variation_plot.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/somatic_variation_plot.log"),
    params:
        script=os.path.join(
            os.path.dirname(workflow.basedir),
            "scripts/somatic_variation_plot.py",
        ),
        sample_info="/home/wdecoster/chr15q14/full_cohort_for_paper.tsv",
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        "python {params.script} -o {output} --sample_info {params.sample_info} {input} 2> {log}"

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
        minlength=100,
        sampleinfo="/home/wdecoster/chr15q14/full_cohort_for_paper.tsv",
    log:
        os.path.join(outdir, "logs/{target}/kmer_plot.log"),
    shell:
        """
        python {params.script} \
        -k {params.kmer} \
        --output {output.html} \
        --counts {output.counts} \
        --minlength {params.minlength} \
        --sampleinfo {params.sampleinfo} \
        {input} 2> {log}"""


rule ct_vs_length:
    """
    Create a scatter plot of the CT repeat length versus the total repeat length.
    Additionally, create an overview table
    """
    input:
        vcfs = expand(
            os.path.join(outdir, "strdust/{{target}}/{id}.vcf.gz"),
            id=crams.loc[crams["collection"].isin(["normal", "1000G", "owen", "kristel", "mayo"]), "individual"], # this corresponds to the full cohort, without the relatives
        ),
        copy_number = os.path.join(outdir, "mosdepth/copy_number.tsv"),
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
        yline = 450,
        arrows = "rr_UCL2783,rr_CW05_004"
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """
        python {params.script} \
        --output {output.plot} \
        --overview {output.overview} \
        --sampleinfo {params.sample_info} \
        --copy_number {input.copy_number} \
        --haplotypes \
        --xline {params.xline} \
        --yline {params.yline} \
        --arrow {params.arrows} \
        {input.vcfs} 2> {log} 
        """

rule ct_dimer_strip:
    input:
        os.path.join(outdir, "analysis_overview_{target}.tsv"),
    output:
        os.path.join(outdir, "plots/{target}/ct_dimer_strip.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/ct_dimer_strip.log"),
    params:
        script=os.path.join(
            os.path.dirname(workflow.basedir),
            "scripts/plot_ct_dimer_strip.py",
        ),
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """
        python {params.script} -i {input} -o {output} 2> {log}
        """

rule correlate_with_age:
    input:
        os.path.join(outdir, "analysis_overview_{target}.tsv")
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
        os.path.join(outdir, "analysis_overview_{target}.tsv")
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

rule table_of_carriers:
    """
    This rule uses the overview file(s) to create a table of haplotype carriers (major and minor) and formats it for publication in the supplementary data of the paper.
    """
    input:
        overview = os.path.join(outdir, "analysis_overview_{target}.tsv"),
    output:
        os.path.join(outdir, "tables/haplotype_carriers_{target}.xlsx"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/haplotype_carriers.log"),
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    params:
        script=os.path.join(
            os.path.dirname(workflow.basedir),
            "scripts/table_of_carriers.py",
        ),
    shell:
        """
        python {params.script} -i {input} -o {output} 2> {log}
        """


rule prepare_astronaut_tables:
    """
    This rule slices the analysis_overview table into tables for aSTRonaut, to create the plots for the different groups of individuals.
    """
    input:
        overview = os.path.join(outdir, "analysis_overview_{target}.tsv"),
    output:
        all = os.path.join(outdir, "temp/analysis_overview_{target}_all.tsv"),
        hapA = os.path.join(outdir, "temp/analysis_overview_{target}_hapA.tsv"),
        hapB = os.path.join(outdir, "temp/analysis_overview_{target}_hapB.tsv"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/prepare_astronaut_tables.log"),
    params:
        script=os.path.join(
            os.path.dirname(workflow.basedir),
            "scripts/prepare_astronaut_tables.py",
        ),
        all = crams.loc[crams["collection"] == "normal", "individual"].tolist(), # this corresponds to our own cohort
        hapA = crams.loc[(crams["collection"] == "normal") & (crams["haplotype"] == 'major'), "individual"].tolist(), # this corresponds to our own cohort, major haplotype carriers
        hapB = crams.loc[(crams["collection"] == "normal") & (crams["haplotype"] == 'minor'), "individual"].tolist(), # this corresponds to our own cohort, minor haplotype carriers 
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """
        python {params.script} -i {input} \
        --all {params.all} --all_out {output.all} \
        --hapA {params.hapA} --hapA_out {output.hapA} \
        --hapB {params.hapB} --hapB_out {output.hapB} \
        2> {log}
        """

rule astronaut_all:
    input:
        os.path.join(outdir, "temp/analysis_overview_{target}_all.tsv"),
    output:
        os.path.join(outdir, "plots/{target}/aSTRonaut_all.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/aSTRonaut_all.log"),
    params:
        script="/home/wdecoster/pathSTR-1000G/scripts/aSTRonaut.py",
        sample_info="/home/wdecoster/chr15q14/full_cohort_for_paper.tsv",
        dotsize=8,
        minsize = 100,
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """
        python {params.script} \
        --motifs CT,CCTT,CTTT,CCCT,CCCTCT,CCCCT,CCTTT,CCCCCC \
        -m {params.minsize} \
        -o {output} \
        --size {params.dotsize} \
        --hide-labels --longest_only --publication \
        --title "Repeat composition sequence" \
        --sampleinfo {params.sample_info} \
        --height 1200 \
        --table {input} 2> {log}
        """


rule astronaut_hapA:
    input:
        os.path.join(outdir, "temp/analysis_overview_{target}_hapA.tsv"),
    output:
        os.path.join(outdir, "plots/{target}/aSTRonaut_hapA.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/aSTRonaut_hapA.log"),
    params:
        script="/home/wdecoster/pathSTR-1000G/scripts/aSTRonaut.py",
        sample_info="/home/wdecoster/chr15q14/full_cohort_for_paper.tsv",
        dotsize = 8,
        minsize = 100,
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """
        python {params.script} \
        --motifs CT,CCTT,CTTT,CCCT,CCCTCT,CCCCT,CCTTT,CCCCCC \
        -m {params.minsize} \
        -o {output} \
        --size {params.dotsize} \
        --hide-labels --longest_only --publication \
        --title "Repeat composition sequence in carriers of haplotype A" \
        --sampleinfo {params.sample_info} \
        --table {input} 2> {log}
        """


rule astronaut_hapB:
    input:
        os.path.join(outdir, "temp/analysis_overview_{target}_hapB.tsv"),        
    output:
        os.path.join(outdir, "plots/{target}/aSTRonaut_hapB.html"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/aSTRonaut_hapB.log"),
    params:
        script="/home/wdecoster/pathSTR-1000G/scripts/aSTRonaut.py",
        sample_info="/home/wdecoster/chr15q14/full_cohort_for_paper.tsv",
        dotsize = 8,
        minsize = 100,
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """
        python {params.script} \
        --motifs CT,CCTT,CTTT,CCCT,CCCTCT,CCCCT,CCTTT,CCCCCC \
        -m {params.minsize} \
        -o {output} \
        --size {params.dotsize} \
        --hide-labels --longest_only --publication \
        --title "Repeat composition sequence in carriers of haplotype B" \
        --sampleinfo {params.sample_info} \
        --table {input} 2> {log}
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
        minsize = 100,
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """
        python {params.script} \
        --motifs CT,CCTT,CTTT,CCCT,CCCTCT,CCCCT,CCTTT,CCCCCC \
        --names {params.relative_names} \
        --label_size 20 \
        --alphabetic \
        --hide_allele_label \
        --size {params.dotsize} \
        {input} \
        --legend_corner topright \
        -m {params.minsize} \
        -o {output} --longest_only --publication --title "Repeat composition sequence in relatives" --sampleinfo {params.sample_info} 2> {log}
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
        duplicate_names=fix_names_duplicates,
        minsize =100,
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        """
        python {params.script} \
        --motifs CT,CCTT,CTTT,CCCT,CCCTCT,CCCCT,CCTTT,CCCCCC \
        --names {params.duplicate_names} \
        --label_size 20 \
        --alphabetic \
        --hide_allele_label \
        {input} -m {params.minsize} -o {output} --longest_only --publication --title "Repeat composition sequence in individuals with multiple samples" --sampleinfo {params.sample_info} 2> {log}
        """

rule precision_recall:
    input:
        overview = os.path.join(outdir, "analysis_overview_{target}.tsv")
    output:
        os.path.join(outdir, "plots/{target}/precision_recall.txt"),
    log:
        os.path.join(outdir, "logs/workflows/{target}/precision_recall.log"),
    params:
        script=os.path.join(
            os.path.dirname(workflow.basedir),
            "scripts/precision_recall.py",
        ),
        cutoff_CT_dimer = 190,
        cutoff_double_length = 450,
        cutoff_double_ct = 0.8,
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/precision_recall.yml")
    shell:
        """
        python {params.script} \
        --data {input.overview} \
        --cutoff_CT_dimer {params.cutoff_CT_dimer} \
        --cutoff_double {params.cutoff_double_length} {params.cutoff_double_ct} > {output} 2> {log}
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
            id=crams.loc[crams["collection"].isin(["normal", "1000G", "owen", "kristel", "mayo"]), "individual"], # this corresponds to the full cohort
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
    
rule cramino:
    input:
        cram=get_cram,
    output:
        os.path.join(outdir, "cramino/{id}.cramino"),
    log:
        os.path.join(outdir, "logs/workflows/{id}-cramino.log"),
    params:
        executable="/home/wdecoster/repositories/cramino/target/release/cramino",
        ref=ref,
    threads:
        4
    shell:
        """
        {params.executable} --karyotype --reference {params.ref} --threads {threads} {input.cram} > {output} 2> {log}
        """

rule gather_cramino:
    input:
        expand(os.path.join(outdir, "cramino/{id}.cramino"), id=crams.loc[crams["collection"] == "normal", "individual"])
    output:
        os.path.join(outdir, "cramino/cramino-all.txt"),
    log:
        os.path.join(outdir, "logs/workflows/gather_cramino.log"),
    params:
        script = os.path.join(os.path.dirname(workflow.basedir), "scripts/gather_cramino.py"),
    shell:
        "python {params.script} -i {input} -o {output} 2> {log}"

rule sex_check:
    input:
        os.path.join(outdir, "cramino/cramino-all.txt"),
    output:
        os.path.join(outdir, "sex_check.html"),
    log:
        os.path.join(outdir, "logs/sex_check.log"),
    params:
        script = os.path.join(os.path.dirname(workflow.basedir), "scripts/sex_check.py"),
        cohort = "/home/wdecoster/cohorts/Individuals.xlsx",
    conda:
        os.path.join(os.path.dirname(workflow.basedir), "envs/pandas_cyvcf2_plotly.yml")
    shell:
        "python {params.script} -c {input} -o {output} --sampleinfo {params.cohort} 2> {log}"