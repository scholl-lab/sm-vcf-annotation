import os
import glob
import functools
import json

# ----------------------------------------------------------------------------------- #
# Script Description:
# This script orchestrates a workflow for annotating VCF files using Snakemake. 
# The workflow has three main annotation steps:
#   1. Annotation through snpEff.
#   2. Adding variant type information using SnpSift varType.
#   3. Detailed annotation using SnpSift with the dbNSFP database.
#   4. Optional: Adding extra annotations from arbitrary VCF files.
# Intermediate files generated during the workflow can be deleted to save space.
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Read the configuration file
configfile: "config.yaml"
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Define temporary directory using an environment variable (usually set by the cluster scheduler)
SCRATCH_DIR = os.environ.get('TMPDIR')
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
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Define result directories using functools.partial to join paths with the output folder
prefix_results = functools.partial(os.path.join, OUTPUT_FOLDER)
annotation_dir = prefix_results(ANNOTATION_SUBFOLDER)
log_dir = prefix_results(LOG_SUBFOLDER)
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Helper Functions:
# - get_vcf_files: Retrieves all VCF files from the specified input folder.
# - get_mem_from_threads: Calculates the amount of memory to allocate based on the number of threads.
def get_vcf_files():
    """Retrieve all VCF files in the specified input folder."""
    return glob.glob(os.path.join(VCF_INPUT_FOLDER, "*.vcf.gz"))

def get_mem_from_threads(wildcards, threads):
    """Calculate the amount of memory to allocate based on the number of threads."""
    return 2048 * threads

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
# Snakemake Rules:
# The rules define the workflow operations, from annotating VCF files using snpEff and SnpSift 
# to cleaning up intermediate files generated during the workflow.

# Rule "all": Specifies the final outputs of the workflow, i.e., the annotated VCF files.
rule all:
    input:
        expand("{annotation_dir}/{sample}.annotated.vcf.gz", 
               annotation_dir=annotation_dir, 
               sample=[os.path.basename(x).replace('.vcf.gz', '') for x in get_vcf_files()])

# Rule "snpeff_annotation": Annotates VCF files using snpEff.
rule snpeff_annotation:
    input:
        vcf_file = lambda wildcards: os.path.join(VCF_INPUT_FOLDER, f"{wildcards.sample}.vcf.gz")
    output:
        ann_vcf = os.path.join(annotation_dir, '{sample}.ann.vcf.gz')
    log:
        os.path.join(log_dir, 'snpeff_annotation.{sample}.log')
    conda:
        CONDA_ENVIRONMENT_ANNOTATION
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,      # Memory in MB based on the number of threads
        time = '72:00:00',                  # Time limit for the job
        tmpdir = SCRATCH_DIR,               # Temporary directory          
    shell:
        '''
        echo "Starting snpeff_annotation at: $(date)" >> {log}
        mkdir -p {annotation_dir}
        snpEff -Xms4000m -Xmx8g {SNPEFF_ANNOTATION_DB} {SNPEFF_ADDITIONAL_FLAGS} -stats {output}.html {input.vcf_file} | bgzip -c > {output} 2>> {log}
        bcftools index --threads {threads} -t {output} 2>> {log}
        echo "Finished snpeff_annotation at: $(date)" >> {log}
        '''

# Rule "snpsift_variant_type": Adds variant type information using SnpSift varType.
rule snpsift_variant_type:
    input:
        ann_vcf = os.path.join(annotation_dir, '{sample}.ann.vcf.gz')
    output:
        ann_vartype_vcf = os.path.join(annotation_dir, '{sample}.ann.vartype.vcf.gz')
    log:
        os.path.join(log_dir, 'snpsift_variant_type.{sample}.log')
    conda:
        CONDA_ENVIRONMENT_ANNOTATION
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,      # Memory in MB based on the number of threads
        time = '72:00:00',                  # Time limit for the job
        tmpdir = SCRATCH_DIR,               # Temporary directory        
    shell:
        '''
        echo "Starting snpsift_variant_type at: $(date)" >> {log}
        SnpSift -Xms4000m -Xmx8g varType {input.ann_vcf} | bgzip -c > {output} 2>> {log}
        bcftools index --threads {threads} -t {output} 2>> {log}
        echo "Finished snpsift_variant_type at: $(date)" >> {log}
        '''

# Rule "snpsift_annotation_dbnsfp": Further annotates the files using SnpSift with the dbNSFP database.
rule snpsift_annotation_dbnsfp:
    input:
        ann_vartype_vcf = os.path.join(annotation_dir, '{sample}.ann.vartype.vcf.gz')
    output:
        ann_dbnsfp_vcf = os.path.join(annotation_dir, '{sample}.ann.dbnsfp.vcf.gz')
    log:
        os.path.join(log_dir, 'snpsift_annotation_dbnsfp.{sample}.log')
    conda:
        CONDA_ENVIRONMENT_ANNOTATION
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,      # Memory in MB based on the number of threads
        time = '72:00:00',                  # Time limit for the job
        tmpdir = SCRATCH_DIR,               # Temporary directory        
    shell:
        '''
        echo "Starting snpsift_annotation_dbnsfp at: $(date)" >> {log}
        SnpSift -Xms4000m -Xmx8g dbnsfp -f {SNPSIFT_DBNSFP_FIELDS} -db {SNPSIFT_DB_LOCATION} {input.ann_vartype_vcf} | bgzip -c > {output} 2>> {log}
        bcftools index --threads {threads} -t {output} 2>> {log}
        echo "Finished snpsift_annotation_dbnsfp at: $(date)" >> {log}
        '''

# Rule "extra_annotations": Adds arbitrary VCF annotations using SnpSift annotate.
rule extra_annotations:
    input:
        vcf = os.path.join(annotation_dir, '{sample}.ann.dbnsfp.vcf.gz')
    output:
        vcf = os.path.join(annotation_dir, '{sample}.annotated.vcf.gz')
    log:
        os.path.join(log_dir, 'extra_annotations.{sample}.log')
    conda:
        CONDA_ENVIRONMENT_ANNOTATION
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,      # Memory in MB based on the number of threads
        time = '72:00:00',                  # Time limit for the job
        tmpdir = SCRATCH_DIR,               # Temporary directory
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

# Rule "clean_intermediate_files": Deletes the intermediate files generated during the annotation to free up space.
rule clean_intermediate_files:
    input:
        "{annotation_dir}/{sample}.ann.vcf.gz"
    shell:
        '''
        rm {input} {input}.tbi
        '''
# ----------------------------------------------------------------------------------- #
