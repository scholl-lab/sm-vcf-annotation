import os
import glob
import functools
import json

# ----------------------------------------------------------------------------------- #
# Script Description:
# This script integrates a complete annotation pipeline for large multisample VCFs:
#   1. Generate equally sized intervals using GATK ScatterIntervalsByNs and SplitIntervals.
#   2. Scatter the input VCF into these intervals, excluding non-canonical contigs.
#   3. Annotate each interval-chunked VCF with:
#       - snpEff for variant annotation
#       - SnpSift varType to determine variant types
#       - SnpSift dbnsfp for additional functional annotations from dbNSFP
#       - Optional extra VCF-based annotations
#   4. Concatenate all annotated interval files into a single final annotated VCF.
#
# Intermediate files can be cleaned up if desired.
#
# The workflow uses a "config.yaml" for parameters.
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Read the configuration file
configfile: "config.yaml"
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Define temporary directory using an environment variable (usually set by the cluster scheduler)
SCRATCH_DIR = os.environ.get('TMPDIR', '/tmp')
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Extract user-defined parameters from the configuration file
OUTPUT_FOLDER = config["output_folder"]
VCF_INPUT_FOLDER = config["vcf_input_folder"]
SNPEFF_ANNOTATION_DB = config["snpeff_annotation_db"]
SNPEFF_ADDITIONAL_FLAGS = config["snpeff_additional_flags"]
SNPSIFT_DB_LOCATION = config["snpsift_db_location"]
ANNOTATION_SUBFOLDER = config["annotation_subfolder"]
LOG_SUBFOLDER = config["log_subfolder"]
CONDA_ENVIRONMENT_ANNOTATION = config["conda_environment_annotation"]
SNPSIFT_DBNSFP_FIELDS = config["snpsift_dbnsfp_fields"]
EXTRA_VCF_ANNOTATIONS = config.get("extra_vcf_annotations", [])

REFERENCE_GENOME = config["reference_genome"]
REFERENCE_DICT = config["reference_dict"]
SCATTER_COUNT = int(config.get("scatter_count", 100))
GATK_ENV = config.get("gatk_env", CONDA_ENVIRONMENT_ANNOTATION)  # Environment with GATK installed

# Added canonical contigs parameter
CANONICAL_CONTIGS = config["canonical_contigs"]
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# DIRECTORY SETUP:
# Define directories
prefix_results = functools.partial(os.path.join, OUTPUT_FOLDER)
annotation_dir = prefix_results(ANNOTATION_SUBFOLDER)
log_dir = prefix_results(LOG_SUBFOLDER)
intervals_dir = prefix_results("intervals")

os.makedirs(intervals_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Helper Functions
def get_vcf_files():
    """Retrieve all VCF files in the specified input folder."""
    return glob.glob(os.path.join(VCF_INPUT_FOLDER, "*.vcf.gz"))

def get_mem_from_threads(wildcards, threads):
    """Calculate the amount of memory to allocate based on the number of threads."""
    return 4096 * threads

def format_extra_annotations(annotations):
    """Format extra annotations into a space-separated string."""
    formatted_annotations = []
    for annotation in annotations:
        formatted_annotations.append(
            f"{annotation['vcf_file']},{annotation['info_field']},{annotation['annotation_prefix']}"
        )
    return ' '.join(formatted_annotations)

formatted_extra_annotations = format_extra_annotations(EXTRA_VCF_ANNOTATIONS)

# Added function to parse dict and extract non-canonical contigs
def get_noncanonical_contigs(dict_file, canonical):
    """Parse the reference dictionary to find non-canonical contigs."""
    noncanonical = []
    with open(dict_file, 'r') as f:
        for line in f:
            if line.startswith('@SQ'):
                fields = line.strip().split('\t')
                for field in fields:
                    if field.startswith('SN:'):
                        contig = field.replace('SN:', '')
                        if contig not in canonical:
                            noncanonical.append(contig)
    return noncanonical

noncanonical_contigs = get_noncanonical_contigs(REFERENCE_DICT, CANONICAL_CONTIGS)
exclude_intervals = ' '.join(["-XL " + c for c in noncanonical_contigs])
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# SAMPLE AND INTERVAL SETUP:
# Determine the input samples and interval IDs
samples = [os.path.basename(x).replace('.vcf.gz', '') for x in get_vcf_files()]
interval_ids = [str(i).zfill(4) for i in range(0, SCATTER_COUNT)]

# Set intervals list path (always genome-style scattering)
INTERVALS_LIST = os.path.join(intervals_dir, "intervals.interval_list")
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# SNAKEMAKE RULES:
# The workflow is defined below. Each rule corresponds to a step in the pipeline.
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# RULE: all
# The 'all' rule specifies the final targets. It ensures that the pipeline produces:
#   - A single fully annotated VCF per sample.
rule all:
    input:
        expand(os.path.join(annotation_dir, "{sample}.annotated.vcf.gz"), sample=samples)
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# RULE: scatter_intervals_by_ns
# Uses GATK ScatterIntervalsByNs to identify suitable intervals for scattering the genome.
# This creates a primary interval list without subdividing it into smaller chunks yet.
rule scatter_intervals_by_ns:
    input:
        reference=REFERENCE_GENOME
    output:
        intervals_list=INTERVALS_LIST
    log:
        os.path.join(log_dir, "scatter_intervals_by_ns.log")
    conda:
        GATK_ENV
    shell:
        """
        set -e
        echo "Starting ScatterIntervalsByNs at: $(date)" > {log}
        gatk ScatterIntervalsByNs \
            -R {input.reference} \
            -O {output.intervals_list} &>> {log}
        echo "Finished ScatterIntervalsByNs at: $(date)" >> {log}
        """

# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# RULE: split_intervals
# Splits the primary intervals generated above into SCATTER_COUNT approximately equal parts.
# This ensures intervals are roughly similar in size for balanced parallelization.
rule split_intervals:
    input:
        intervals_list=INTERVALS_LIST,
        reference=REFERENCE_GENOME
    output:
        interval_files=expand(os.path.join(intervals_dir, "{interval_id}-scattered.interval_list"),
                              interval_id=interval_ids)
    params:
        scatter_count=SCATTER_COUNT,
        exclude_intervals=exclude_intervals
    log:
        os.path.join(log_dir, "split_intervals.log")
    conda:
        GATK_ENV
    shell:
        """
        set -e
        echo "Starting SplitIntervals at: $(date)" > {log}
        gatk SplitIntervals \
            -R {input.reference} \
            -L {input.intervals_list} \
            {params.exclude_intervals} \
            --scatter-count {params.scatter_count} \
            -O {intervals_dir} &>> {log}
        echo "Finished SplitIntervals at: $(date)" >> {log}
        """
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# RULE: scatter_vcf
# Extracts variants from the original large VCF for each scattered interval.
# This produces one VCF per interval per sample.
rule scatter_vcf:
    input:
        vcf_file=lambda w: os.path.join(VCF_INPUT_FOLDER, f"{w.sample}.vcf.gz"),
        interval_file=lambda w: os.path.join(intervals_dir, f"{w.interval_id}-scattered.interval_list")
    output:
        os.path.join(annotation_dir, "{sample}.{interval_id}.vcf.gz")
    log:
        os.path.join(log_dir, "scatter_vcf.{sample}.{interval_id}.log")
    conda:
        GATK_ENV
    threads: 4
    resources:
        mem_mb=get_mem_from_threads,
        time='72:00:00',
        tmpdir=SCRATCH_DIR
    shell:
        r'''
        echo "Starting scatter_vcf at: $(date)" >> {log}
        gatk --java-options "-Xmx8g" SelectVariants \
          -R {REFERENCE_GENOME} \
          -V {input.vcf_file} \
          -L {input.interval_file} \
          -O {output} 2>> {log}
        echo "Finished scatter_vcf at: $(date)" >> {log}
        '''
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# RULE: snpeff_annotation
# Annotates the variants in each interval-chunked VCF using snpEff.
# Adds gene, transcript, and other functional annotations.
rule snpeff_annotation:
    input:
        vcf_file=os.path.join(annotation_dir, "{sample}.{interval_id}.vcf.gz")
    output:
        ann_vcf=os.path.join(annotation_dir, "{sample}.{interval_id}.ann.vcf.gz")
    log:
        os.path.join(log_dir, "snpeff_annotation.{sample}.{interval_id}.log")
    conda:
        CONDA_ENVIRONMENT_ANNOTATION
    threads: 4
    resources:
        mem_mb=get_mem_from_threads,
        time='72:00:00',
        tmpdir=SCRATCH_DIR
    shell:
        '''
        echo "Starting snpeff_annotation at: $(date)" >> {log}
        snpEff -Xms4000m -Xmx8g {SNPEFF_ANNOTATION_DB} {SNPEFF_ADDITIONAL_FLAGS} -stats {output.ann_vcf}.html {input.vcf_file} | bgzip -c > {output.ann_vcf} 2>> {log}
        bcftools index --threads {threads} -t {output.ann_vcf} 2>> {log}
        echo "Finished snpeff_annotation at: $(date)" >> {log}
        '''
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# RULE: snpsift_variant_type
# Uses SnpSift varType to determine variant types (e.g., SNP, INS, DEL).
rule snpsift_variant_type:
    input:
        ann_vcf=os.path.join(annotation_dir, "{sample}.{interval_id}.ann.vcf.gz")
    output:
        ann_vartype_vcf=os.path.join(annotation_dir, "{sample}.{interval_id}.ann.vartype.vcf.gz")
    log:
        os.path.join(log_dir, "snpsift_variant_type.{sample}.{interval_id}.log")
    conda:
        CONDA_ENVIRONMENT_ANNOTATION
    threads: 4
    resources:
        mem_mb=get_mem_from_threads,
        time='72:00:00',
        tmpdir=SCRATCH_DIR
    shell:
        '''
        echo "Starting snpsift_variant_type at: $(date)" >> {log}
        # Check if the VCF has any variants (non-header lines).
        variant_count=$(bcftools view -H {input.ann_vcf} | wc -l)
        if [ "$variant_count" -eq 0 ]; then
            echo "No variants found in {input.ann_vcf}, copying input to output." >> {log}
            cp {input.ann_vcf} {output.ann_vartype_vcf}
            cp {input.ann_vcf}.tbi {output.ann_vartype_vcf}.tbi
        else
            echo "Variants found, running SnpSift varType." >> {log}
            SnpSift -Xms4000m -Xmx8g varType {input.ann_vcf} | bgzip -c > {output.ann_vartype_vcf} 2>> {log}
            bcftools index --threads {threads} -t {output.ann_vartype_vcf} 2>> {log}
        fi
        echo "Finished snpsift_variant_type at: $(date)" >> {log}
        '''
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# RULE: snpsift_annotation_dbnsfp
# Uses SnpSift dbnsfp to add functional scores and predictions from dbNSFP.
rule snpsift_annotation_dbnsfp:
    input:
        ann_vartype_vcf=os.path.join(annotation_dir, "{sample}.{interval_id}.ann.vartype.vcf.gz")
    output:
        ann_dbnsfp_vcf=os.path.join(annotation_dir, "{sample}.{interval_id}.ann.dbnsfp.vcf.gz")
    log:
        os.path.join(log_dir, "snpsift_annotation_dbnsfp.{sample}.{interval_id}.log")
    conda:
        CONDA_ENVIRONMENT_ANNOTATION
    threads: 4
    resources:
        mem_mb=get_mem_from_threads,
        time='72:00:00',
        tmpdir=SCRATCH_DIR
    shell:
        '''
        echo "Starting snpsift_annotation_dbnsfp at: $(date)" >> {log}
        SnpSift -Xms4000m -Xmx8g dbnsfp -f {SNPSIFT_DBNSFP_FIELDS} -db {SNPSIFT_DB_LOCATION} {input.ann_vartype_vcf} | bgzip -c > {output.ann_dbnsfp_vcf} 2>> {log}
        bcftools index --threads {threads} -t {output.ann_dbnsfp_vcf} 2>> {log}
        echo "Finished snpsift_annotation_dbnsfp at: $(date)" >> {log}
        '''
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# RULE: extra_annotations
# Optionally adds other annotations from provided VCF files using SnpSift annotate.
rule extra_annotations:
    input:
        vcf=os.path.join(annotation_dir, "{sample}.{interval_id}.ann.dbnsfp.vcf.gz")
    output:
        vcf=os.path.join(annotation_dir, "{sample}.{interval_id}.annotated.vcf.gz")
    log:
        os.path.join(log_dir, "extra_annotations.{sample}.{interval_id}.log")
    conda:
        CONDA_ENVIRONMENT_ANNOTATION
    threads: 4
    resources:
        mem_mb=get_mem_from_threads,
        time='72:00:00',
        tmpdir=SCRATCH_DIR
    params:
        annotations=formatted_extra_annotations
    shell:
        '''
        echo "Starting extra_annotations at: $(date)" >> {log}
        final_vcf={input.vcf}
        for vcf_annotation in {params.annotations}; do
            IFS=',' read -r vcf_file info_field annotation_prefix <<< "$vcf_annotation"
            output_vcf=$(echo $final_vcf | sed 's/.vcf.gz/.tmp.vcf.gz/')
            SnpSift -Xms4000m -Xmx8g annotate -info $info_field -name $annotation_prefix $vcf_file $final_vcf | bgzip -c > $output_vcf
            final_vcf=$output_vcf
        done
        mv $final_vcf {output.vcf}
        bcftools index --threads {threads} -t {output.vcf} 2>> {log}
        echo "Finished extra_annotations at: $(date)" >> {log}
        '''
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# RULE: concatenate_annotated_vcfs
# After all intervals are annotated, this rule concatenates them into a single VCF per sample.
rule concatenate_annotated_vcfs:
    input:
        lambda w: expand(os.path.join(annotation_dir, "{sample}.{interval_id}.annotated.vcf.gz"),
                         sample=[w.sample],
                         interval_id=interval_ids)
    output:
        concatenated_vcf=os.path.join(annotation_dir, "{sample}.annotated.vcf.gz")
    log:
        os.path.join(log_dir, "concatenate_annotated_vcfs.{sample}.log")
    conda:
        CONDA_ENVIRONMENT_ANNOTATION
    threads: 4
    resources:
        mem_mb=get_mem_from_threads,
        time='72:00:00',
        tmpdir=SCRATCH_DIR
    shell:
        '''
        echo "Starting concatenate_annotated_vcfs at: $(date)" >> {log}
        bcftools concat -Oz --threads {threads} {input} -o {output.concatenated_vcf} 2>> {log}
        bcftools index --threads {threads} -t {output.concatenated_vcf} 2>> {log}
        echo "Finished concatenate_annotated_vcfs at: $(date)" >> {log}
        '''
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# RULE: clean_intermediate_files
# (Optional) Cleans up intermediate files after the final concatenation is done.
# This frees up disk space by removing per-interval annotated VCFs.
rule clean_intermediate_files:
    input:
        os.path.join(annotation_dir, "{sample}.{interval_id}.ann.vcf.gz")
    shell:
        '''
        rm {input} {input}.tbi
        '''
# ----------------------------------------------------------------------------------- #
