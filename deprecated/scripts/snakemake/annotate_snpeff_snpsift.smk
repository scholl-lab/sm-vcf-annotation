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
    return 8192 * threads

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

# ----------------------------------------------------------------------------------- #
# STEP-BY-STEP EXTRA ANNOTATIONS
# ----------------------------------------------------------------------------------- #

N_ANNOTATIONS = len(EXTRA_VCF_ANNOTATIONS)

def annotation_step_vcf(sample, step):
    """Return the path for a given step of extra annotation."""
    if step == 0:
        # Step 0 => input is the VCF after dbNSFP
        return os.path.join(annotation_dir, f"{sample}.ann.dbnsfp.vcf.gz")
    # step 1..N => intermediate files
    return os.path.join(annotation_dir, f"{sample}.ann.step{step}.vcf.gz")


rule annotation_step:
    """
    Applies exactly ONE annotation from EXTRA_VCF_ANNOTATIONS by step index.
    step ranges from 1..N_ANNOTATIONS, referencing the i-th annotation in the config.
    """
    input:
        lambda wc: annotation_step_vcf(wc.sample, int(wc.step) - 1)
    output:
        annotation_dir + "/{sample}.ann.step{step}.vcf.gz"
    log:
        log_dir + "/extra_annotations.step{step}.{sample}.log"
    params:
        # We pass the annotation details to shell so we can run SnpSift in conda environment
        step_index = lambda wc: int(wc.step) - 1,
        vcf_file = lambda wc: EXTRA_VCF_ANNOTATIONS[int(wc.step) - 1]["vcf_file"],
        info_field = lambda wc: EXTRA_VCF_ANNOTATIONS[int(wc.step) - 1]["info_field"],
        annotation_prefix = lambda wc: EXTRA_VCF_ANNOTATIONS[int(wc.step) - 1]["annotation_prefix"],
    conda:
        CONDA_ENVIRONMENT_ANNOTATION
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,      # Memory in MB based on the number of threads
        time = '72:00:00',                  # Time limit for the job
        tmpdir = SCRATCH_DIR,               # Temporary directory
    shell:
        r'''
        echo "Starting annotation_step {wildcards.step} at: $(date)" >> {log}

        SnpSift -Xms4000m -Xmx8g annotate \
          -info {params.info_field} \
          -name {params.annotation_prefix} \
          {params.vcf_file} {input} \
        | bgzip -c > {output} 2>> {log}

        bcftools index --threads {threads} -t {output} 2>> {log}

        echo "Finished annotation_step {wildcards.step} at: $(date)" >> {log}
        '''


rule extra_annotations_final:
    """
    After the last step (step=N_ANNOTATIONS) is done, rename final file to .annotated.vcf.gz.
    """
    input:
        lambda wc: annotation_step_vcf(wc.sample, N_ANNOTATIONS)
    output:
        annotation_dir + "/{sample}.annotated.vcf.gz"
    log:
        log_dir + "/extra_annotations.final.{sample}.log"
    conda:
        CONDA_ENVIRONMENT_ANNOTATION
    threads: 4
    resources:
        mem_mb = get_mem_from_threads,      # Memory in MB based on the number of threads
        time = '72:00:00',                  # Time limit for the job
        tmpdir = SCRATCH_DIR,               # Temporary directory
    shell:
        r'''
        echo "Renaming final-step VCF to .annotated.vcf.gz at: $(date)" > {log}
        mv {input} {output}
        bcftools index --threads {threads} -t {output} 2>> {log}
        echo "Done final rename at: $(date)" >> {log}
        '''


# Handle the case when there are no extra annotations defined
if N_ANNOTATIONS == 0:
    rule no_extra_annotations:
        input:
            os.path.join(annotation_dir, "{sample}.ann.dbnsfp.vcf.gz")
        output:
            os.path.join(annotation_dir, "{sample}.annotated.vcf.gz")
        log:
            os.path.join(log_dir, "no_extra_annotations.{sample}.log")
        conda:
            CONDA_ENVIRONMENT_ANNOTATION
        threads: 4
        resources:
            mem_mb = get_mem_from_threads,
            time = '1:00:00',
            tmpdir = SCRATCH_DIR
        shell:
            r'''
            echo "No extra annotations configured, just renaming file at: $(date)" > {log}
            cp {input} {output}
            bcftools index --threads {threads} -t {output} 2>> {log}
            echo "Done at: $(date)" >> {log}
            '''

# Rule "clean_intermediate_files": Deletes all intermediate files generated during the annotation to free up space.
rule clean_intermediate_files:
    input:
        final_vcf = os.path.join(annotation_dir, "{sample}.annotated.vcf.gz")
    params:
        annotation_dir = annotation_dir,
        sample = "{sample}"
    shell:
        '''
        # Remove the initial snpEff output
        rm -f {params.annotation_dir}/{params.sample}.ann.vcf.gz {params.annotation_dir}/{params.sample}.ann.vcf.gz.tbi {params.annotation_dir}/{params.sample}.ann.vcf.gz.html
        
        # Remove the vartype annotated files
        rm -f {params.annotation_dir}/{params.sample}.ann.vartype.vcf.gz {params.annotation_dir}/{params.sample}.ann.vartype.vcf.gz.tbi
        
        # Remove the dbnsfp annotated files
        rm -f {params.annotation_dir}/{params.sample}.ann.dbnsfp.vcf.gz {params.annotation_dir}/{params.sample}.ann.dbnsfp.vcf.gz.tbi
        
        # Remove all step-wise annotation files
        rm -f {params.annotation_dir}/{params.sample}.ann.step*.vcf.gz {params.annotation_dir}/{params.sample}.ann.step*.vcf.gz.tbi
        
        echo "Removed all intermediate files for {params.sample}"
        '''
