# VCF Annotation Pipeline

This repository contains a Snakemake workflow script for annotating VCF files using the snpEff and SnpSift tools. Below, you will find the details of the script, the configuration file, and the resources required to run the pipeline.

## Script Description

The Snakemake script reads the configuration parameters from "config.yaml" and orchestrates a pipeline with the following rules to process VCF files:

1. **snpeff_annotation**: This rule annotates VCF files using snpEff with a specified database and additional flags. The output, including an HTML file containing statistics for each annotation, is compressed using bgzip.

2. **snpsift_annotation_dbnsfp**: This rule takes the output of the `snpeff_annotation` rule and further annotates it using SnpSift with the dbNSFP database specified in the configuration file. The output is then compressed using bgzip.

3. **clean_intermediate_files**: This rule deletes the intermediate files created after the snpEff annotation to conserve storage space.

4. **all**: This rule serves as a checkpoint to ensure all individual annotations are performed on all VCF files in the specified input folder.

## Configuration File

The configuration file, "config.yaml," contains the following user-defined parameters:

```yaml
output_folder: "path/to/output/folder"
vcf_input_folder: "path/to/vcf/input/folder"
snpeff_annotation_db: "path/to/snpeff/annotation/database"
snpeff_additional_flags: "additional flags for snpEff"
snpsift_db_location: "path/to/snpsift/database/location"
annotation_subfolder: "annotation"  # Subfolder name for the annotation outputs
log_subfolder: "logs"  # Subfolder name for the log files
conda_environment_annotation: "annotation"  # Name of the conda environment for annotation
snpsift_dbnsfp_fields: "fields for SnpSift annotation" # no spaces allowed
```

Ensure that the paths and parameters in "config.yaml" are correctly set up before running the script.

## Required Resources

Before running the script, ensure the following:

- snpEff and SnpSift tools are correctly installed and configured, including their respective databases.
- Bcftools is installed and accessible in your conda environment.
- The conda environment specified for annotation is set up with all necessary packages installed.

## Running the Script

To execute the workflow, navigate to the directory containing the Snakemake file and use the following command in your terminal:

```sh
snakemake -s annotate_snpeff_snpsift.smk --use-conda --cores 50
```

Replace `50` with the number of cores you wish to allocate to the workflow.

## Logging

The script logs the start and end times of each rule to facilitate performance profiling and troubleshooting. Log files are stored in the directory specified under `log_subfolder` in the configuration file.

## Contribution

Feel free to fork the repository and submit pull requests for any enhancements or bug fixes. Contributions to improve the script or documentation are welcome.

## TODOs/Open Issues for Improvement

1. **Resource Allocation**: Move resource and thread allocation settings, including runtime, to the configuration file or dynamically compute them based on the available system resources.

2. **Java Memory Limits**: Remove hardcoded Java memory limits and allow them to be specified in the configuration file or computed dynamically based on the available memory.

3. **Flexible Input**: Allow more flexible input by accepting non-gzipped VCF files and making the file extension a parameter in the configuration file.

4. **Improved Error Handling**: Implement better error handling to gracefully handle any issues that may arise during execution and provide informative error messages to facilitate debugging.

5. **Documentation**: Enhance the documentation to provide more detailed instructions for setting up and using the various databases required by snpEff and SnpSift, possibly including a guide for setting up the necessary conda environments.

6. **Testing**: Develop a suite of tests to verify that the script is functioning correctly, including tests for different input file formats and annotations.
