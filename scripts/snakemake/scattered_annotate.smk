import os
import glob
import functools
import csv

# ----------------------------------------------------------------------------------- #
# Configuration file specifying user-defined settings
configfile: "config.yaml"
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Define temporary directory using an environment variable (default to 'tmp' if not set)
SCRATCH_DIR = os.environ.get('TMPDIR', 'tmp')

# Extract parameters from the configuration file
OUTPUT_FOLDER = config["output_folder"]
VCF_INPUT_FOLDER = config["vcf_input_folder"]
SNPEFF_ANNOTATION_DB = config["snpeff_annotation_db"]
SNPEFF_ADDITIONAL_FLAGS = config["snpeff_additional_flags"]
SNPSIFT_DB_LOCATION = config["snpsift_db_location"]
ANNOTATION_SUBFOLDER = config["annotation_subfolder"]
LOG_SUBFOLDER = config["log_subfolder"]
CONDA_ENVIRONMENT_ANNOTATION = config["conda_environment_annotation"]
SNPSIFT_DBNSFP_FIELDS = config["snpsift_dbnsfp_fields"]
INPUT_DIR = config["scattered_vcf_folder"]
METADATA_FILE = config["metadata_file"]
VCF_FILE_EXTENSION = config.get("vcf_file_extension", ".vcf.gz")
SCATTER_NAME_PREFIX = config.get("scatter_name_prefix", "")
SCATTER_NAME_DELIMITER = config.get("scatter_name_delimiter", ".")
SCATTERING_INTERVAL_RANGES = config.get("scattering_interval_ranges", [str(i) for i in range(1, 23)] + ['X', 'Y'])

# Define result directories
prefix_results = functools.partial(os.path.join, config['output_folder'])
CONCATENATED_VCF_DIR = prefix_results('concatenated_vcfs')
annotation_dir = prefix_results(ANNOTATION_SUBFOLDER)
log_dir = prefix_results(LOG_SUBFOLDER)

# Ensure output directories exist
os.makedirs(CONCATENATED_VCF_DIR, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Helper Functions

def get_vcf_files():
    """Retrieve all VCF files in the specified input folder."""
    return glob.glob(os.path.join(VCF_INPUT_FOLDER, "*.vcf.gz"))

def get_mem_from_threads(wildcards, threads):
    """Calculate the amount of memory to allocate based on the number of threads."""
    return 2048 * threads

def expand_interval_names(interval_ranges, prefix=""):
    """
    Expand a list of interval ranges into a list of interval names.
    
    Parameters:
    - interval_ranges (list): A list of interval ranges specified as strings.
                              For example: ["1-3", "5", "X-Y"].
    - prefix (str): A prefix to add to each interval name. For example, "chr".
    
    Returns:
    - list: A list of interval names with the prefix added to each name.
    """
    interval_names = []
    for interval_range in interval_ranges:
        if '-' in str(interval_range):
            start, end = map(int, str(interval_range).split('-'))
            interval_names.extend([f"{prefix}{i}" for i in range(start, end+1)])
        else:
            interval_names.append(f"{prefix}{interval_range}")
    return interval_names

# Generate the list of interval names based on the ranges and prefixes specified
SCATTERING_INTERVAL_NAMES = expand_interval_names(SCATTERING_INTERVAL_RANGES, SCATTER_NAME_PREFIX)

# Read the metadata file to build a dictionary with analysis details
metadata_dict = {}
with open(METADATA_FILE, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        analysis_key = f"{row['individual1']}_{row['analysis']}"
        metadata_dict[analysis_key] = row
# ----------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------- #
# Snakemake Rules

# Rule to specify the final output files
rule all:
    input:
        expand(f"{CONCATENATED_VCF_DIR}/{{individual1}}_{{analysis}}{VCF_FILE_EXTENSION}", 
               individual1=[row['individual1'] for row in metadata_dict.values()], 
               analysis=[row['analysis'] for row in metadata_dict.values()])

# Rule to scatter the VCF file into smaller chunks
rule scatter_vcf:
    input:
        vcf_file = lambda wildcards: config["vcf_folder"] + "all_merged.vcf.gz"
    output:
        expand(f"{VCF_INPUT_FOLDER}/all_merged.{interval}.vcf.gz", interval=SCATTERING_INTERVAL_NAMES)
    shell:
        """
        mkdir -p {VCF_INPUT_FOLDER}
        bcftools +scatter {input.vcf_file} -Oz --threads 4 -o {VCF_INPUT_FOLDER} -n 100000 -p all_merged.
        """

# Rule to annotate the scattered VCF files using snpEff
rule snpeff_annotation:
    input:
        vcf_file = os.path.join(VCF_INPUT_FOLDER, '{sample}.vcf.gz')
    output:
        ann_vcf = os.path.join(annotation_dir, '{sample}.ann.vcf.gz')
    log:
        os.path.join(log_dir, 'snpeff_annotation.{sample}.log')
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
        mkdir -p {annotation_dir}
        snpEff -Xms4000m -Xmx8g {SNPEFF_ANNOTATION_DB} {SNPEFF_ADDITIONAL_FLAGS} -stats {output.ann_vcf}.html {input.vcf_file} | bgzip -c > {output.ann_vcf} 2>> {log}
        bcftools index --threads {threads} -t {output.ann_vcf} 2>> {log}
        echo "Finished snpeff_annotation at: $(date)" >> {log}
        '''

# Rule to further annotate the VCF files using SnpSift with the dbNSFP database
rule snpsift_annotation_dbnsfp:
    input:
        ann_vcf = os.path.join(annotation_dir, '{sample}.ann.vcf.gz')
    output:
        ann_dbnsfp_vcf = os.path.join(annotation_dir, '{sample}.ann.dbnsfp.vcf.gz')
    log:
        os.path.join(log_dir, 'snpsift_annotation_dbnsfp.{sample}.log')
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
        SnpSift -Xms4000m -Xmx8g dbnsfp -f {SNPSIFT_DBNSFP_FIELDS} -db {SNPSIFT_DB_LOCATION} {input.ann_vcf} | bgzip -c > {output.ann_dbnsfp_vcf} 2>> {log}
        bcftools index --threads {threads} -t {output.ann_dbnsfp_vcf} 2>> {log}
        echo "Finished snpsift_annotation_dbnsfp at: $(date)" >> {log}
        '''

# Rule to clean up intermediate files generated during the annotation
rule clean_intermediate_files:
    input:
        os.path.join(annotation_dir, '{sample}.ann.vcf.gz')
    shell:
        '''
        rm {input} {input}.tbi
        '''

# Rule to concatenate the annotated VCF files into a single VCF file
rule concatenate_vcfs:
    input:
        lambda wildcards: expand(f"{annotation_dir}/{{individual1}}_{{analysis}}{SCATTER_NAME_DELIMITER}{{interval}}.ann.dbnsfp.vcf.gz", interval=SCATTERING_INTERVAL_NAMES)
    output:
        concatenated_vcf = f"{CONCATENATED_VCF_DIR}/{{individual1}}_{{analysis}}{VCF_FILE_EXTENSION}"
    threads: 2
    resources:
        mem_mb = get_mem_from_threads,
        time = '72:00:00',
        tmpdir = SCRATCH_DIR
    conda:
        "bcftools"
    log:
        f"{LOG_DIR}/{{individual1}}_{{analysis}}_concat.log"
    shell:
        """
        echo "Starting concatenate_vcfs at: $(date)" >> {log}
        bcftools concat {input} -o {output.concatenated_vcf} -Oz 2>> {log}
        bcftools index --threads {threads} -t {output.concatenated_vcf} 2>> {log}
        echo "Finished concatenate_vcfs at: $(date)" >> {log}
        """
# ----------------------------------------------------------------------------------- #
