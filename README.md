# VCF Annotation Pipeline

This repository contains a Snakemake workflow script for annotating VCF files using snpEff and SnpSift tools. Below are the details of the script, configuration file, and the resources required to run the pipeline.

## Script Description

The Snakemake script reads the configuration parameters from "config.yaml" and orchestrates a pipeline with the following rules to process VCF files:

1. **snpeff_annotation**: This rule annotates VCF files using snpEff with a specified database and additional flags. The output is compressed using bgzip, and an HTML file containing statistics for each annotation is generated.
   
2. **snpsift_annotation**: This rule takes the output of the `snpeff_annotation` rule and further annotates it using SnpSift with a specified database file. The output is then compressed using bgzip.

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
```

Ensure that the paths and parameters in "config.yaml" are correctly set up before running the script.

## Required Resources

Before running the script, ensure the following:

- The snpEff and SnpSift tools, along with their respective databases, are correctly installed and accessible in your runtime environment.
- The necessary databases for snpEff and SnpSift are available and correctly configured in your environment.
- Bcftools is correctly installed and accessible in your conda environment.
- The conda environment for annotation is correctly set up and activated.

## Running the Script

Once everything is set up, you can run your workflow using the schedular script with the following command in the terminal (in the directory where your Snakemake file is):

```sh
sbatch run_annotate_snpeff_snpsift.sh 50
```
```

## Logging

The script logs the start and end times of each rule to facilitate performance profiling and diagnosis of potential issues. Log files are stored in the directory specified under `log_subfolder` in the configuration file.

## Contribution

Feel free to fork the repository and submit pull requests for any enhancements or bug fixes.