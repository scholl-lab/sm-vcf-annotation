import os
import glob
import functools
import json

# ----------------------------------------------------------------------------------- #
# Script Description:
# This script integrates:
#   1. Generation of equally sized intervals using GATK ScatterIntervalsByNs and SplitIntervals.
#   2. Scattering the input VCF into these intervals.
#   3. Annotating each interval-chunked VCF with snpEff, SnpSift varType, dbNSFP, and optional extra annotations.
#   4. Concatenating all annotated interval files into a single final annotated VCF.
#
# Intermediate files can be cleaned up to save space.
#
# The workflow uses a configuration file ("config.yaml") for setting parameters.
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
SCATTER_COUNT = int(config.get("scatter_count", 100))
GATK_ENV = config.get("gatk_env", CONDA_ENVIRONMENT_ANNOTATION)  # Environment with GATK installed
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
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
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Determine the input samples
samples = [os.path.basename(x).replace('.vcf.gz', '') for x in get_vcf_files()]

# Generate zero-padded interval IDs
interval_ids = [str(i).zfill(4) for i in range(0, SCATTER_COUNT)]

# Set intervals list path (always genome-style scattering)
INTERVALS_LIST = os.path.join(intervals_dir, "intervals.interval_list")
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Snakemake Rules:
# The rules now handle generating intervals, scattering by intervals, annotation,
# and concatenation.

rule all:
    input:
        # The final outputs are the annotated VCFs for each sample
        expand("{annotation_dir}/{sample}.annotated.vcf.gz", 
               annotation_dir=annotation_dir, 
               sample=samples)

# Rule to generate intervals using ScatterIntervalsByNs.
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

# Rule to split intervals into equal parts using SplitIntervals.
rule split_intervals:
    input:
        intervals_list=INTERVALS_LIST,
        reference=REFERENCE_GENOME
    output:
        interval_files=expand(os.path.join(intervals_dir, "{interval_id}-scattered.interval_list"),
                              interval_id=interval_ids)
    params:
        scatter_count=SCATTER_COUNT
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
            --scatter-count {params.scatter_count} \
            -O {intervals_dir} &>> {log}
        echo "Finished SplitIntervals at: $(date)" >> {log}
        """

# Rule scatter_vcf: Scatter each input VCF into intervals defined above
rule scatter_vcf:
    input:
        vcf_file = lambda wildcards: os.path.join(VCF_INPUT_FOLDER, f"{wildcards.sample}.vcf.gz"),
        intervals = expand(os.path.join(intervals_dir, "{interval_id}-scattered.interval_list"),
                           interval_id=interval_ids)
    output:
        # One VCF per interval_id
        expand("{annotation_dir}/{sample}.{interval_id}.vcf.gz",
               annotation_dir=annotation_dir,
               sample=[lambda wildcards: wildcards.sample],
               interval_id=interval_ids)
    log:
        os.path.join(log_dir, 'scatter_vcf.{sample}.log')
    conda:
        GATK_ENV
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR
    run:
        with open(log[0], 'a') as lf:
            lf.write("Starting scatter_vcf at: {}\n".format(os.popen("date").read().strip()))
        for interval_id in interval_ids:
            interval_list = os.path.join(intervals_dir, f"{interval_id}-scattered.interval_list")
            out_vcf = f"{annotation_dir}/{wildcards.sample}.{interval_id}.vcf.gz"
            shell(
                r'''
                gatk --java-options "-Xmx8g" SelectVariants \
                  -R {REFERENCE_GENOME} \
                  -V {input.vcf_file} \
                  -L {interval_list} \
                  -O {out_vcf} 2>> {log[0]}
                bcftools index --threads {threads} -t {out_vcf} 2>> {log[0]}
                '''
            )
        with open(log[0], 'a') as lf:
            lf.write("Finished scatter_vcf at: {}\n".format(os.popen("date").read().strip()))

# Rule snpeff_annotation: Annotate scattered VCFs using snpEff
rule snpeff_annotation:
    input:
        vcf_file = "{annotation_dir}/{sample}.{interval_id}.vcf.gz"
    output:
        ann_vcf = "{annotation_dir}/{sample}.{interval_id}.ann.vcf.gz"
    log:
        os.path.join(log_dir, 'snpeff_annotation.{sample}.{interval_id}.log')
    conda:
        CONDA_ENVIRONMENT_ANNOTATION
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR          
    shell:
        '''
        echo "Starting snpeff_annotation at: $(date)" >> {log}
        snpEff -Xms4000m -Xmx8g {SNPEFF_ANNOTATION_DB} {SNPEFF_ADDITIONAL_FLAGS} -stats {output.ann_vcf}.html {input.vcf_file} | bgzip -c > {output.ann_vcf} 2>> {log}
        bcftools index --threads {threads} -t {output.ann_vcf} 2>> {log}
        echo "Finished snpeff_annotation at: $(date)" >> {log}
        '''

# Rule snpsift_variant_type: Add variant type info
rule snpsift_variant_type:
    input:
        ann_vcf = "{annotation_dir}/{sample}.{interval_id}.ann.vcf.gz"
    output:
        ann_vartype_vcf = "{annotation_dir}/{sample}.{interval_id}.ann.vartype.vcf.gz"
    log:
        os.path.join(log_dir, 'snpsift_variant_type.{sample}.{interval_id}.log')
    conda:
        CONDA_ENVIRONMENT_ANNOTATION
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR        
    shell:
        '''
        echo "Starting snpsift_variant_type at: $(date)" >> {log}
        SnpSift -Xms4000m -Xmx8g varType {input.ann_vcf} | bgzip -c > {output.ann_vartype_vcf} 2>> {log}
        bcftools index --threads {threads} -t {output.ann_vartype_vcf} 2>> {log}
        echo "Finished snpsift_variant_type at: $(date)" >> {log}
        '''

# Rule snpsift_annotation_dbnsfp: Add dbNSFP annotations
rule snpsift_annotation_dbnsfp:
    input:
        ann_vartype_vcf = "{annotation_dir}/{sample}.{interval_id}.ann.vartype.vcf.gz"
    output:
        ann_dbnsfp_vcf = "{annotation_dir}/{sample}.{interval_id}.ann.dbnsfp.vcf.gz"
    log:
        os.path.join(log_dir, 'snpsift_annotation_dbnsfp.{sample}.{interval_id}.log')
    conda:
        CONDA_ENVIRONMENT_ANNOTATION
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR        
    shell:
        '''
        echo "Starting snpsift_annotation_dbnsfp at: $(date)" >> {log}
        SnpSift -Xms4000m -Xmx8g dbnsfp -f {SNPSIFT_DBNSFP_FIELDS} -db {SNPSIFT_DB_LOCATION} {input.ann_vartype_vcf} | bgzip -c > {output.ann_dbnsfp_vcf} 2>> {log}
        bcftools index --threads {threads} -t {output.ann_dbnsfp_vcf} 2>> {log}
        echo "Finished snpsift_annotation_dbnsfp at: $(date)" >> {log}
        '''

# Rule extra_annotations: Add arbitrary VCF-based annotations
rule extra_annotations:
    input:
        vcf = "{annotation_dir}/{sample}.{interval_id}.ann.dbnsfp.vcf.gz"
    output:
        vcf = "{annotation_dir}/{sample}.{interval_id}.annotated.vcf.gz"
    log:
        os.path.join(log_dir, 'extra_annotations.{sample}.{interval_id}.log')
    conda:
        CONDA_ENVIRONMENT_ANNOTATION
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR
    params:
        annotations = formatted_extra_annotations
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

# Rule concatenate_annotated_vcfs: Concatenate all annotated chunks per sample
rule concatenate_annotated_vcfs:
    input:
        lambda wildcards: expand("{annotation_dir}/{sample}.{interval_id}.annotated.vcf.gz",
                                 annotation_dir=annotation_dir,
                                 sample=[wildcards.sample],
                                 interval_id=interval_ids)
    output:
        concatenated_vcf = "{annotation_dir}/{sample}.annotated.vcf.gz"
    log:
        os.path.join(log_dir, 'concatenate_annotated_vcfs.{sample}.log')
    conda:
        CONDA_ENVIRONMENT_ANNOTATION
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR
    shell:
        '''
        echo "Starting concatenate_annotated_vcfs at: $(date)" >> {log}
        bcftools concat -Oz --threads {threads} {input} -o {output.concatenated_vcf} 2>> {log}
        bcftools index --threads {threads} -t {output.concatenated_vcf} 2>> {log}
        echo "Finished concatenate_annotated_vcfs at: $(date)" >> {log}
        '''

# Rule clean_intermediate_files: Delete intermediate files if desired
rule clean_intermediate_files:
    input:
        "{annotation_dir}/{sample}.{interval_id}.ann.vcf.gz"
    shell:
        '''
        rm {input} {input}.tbi
        '''
# ----------------------------------------------------------------------------------- #
