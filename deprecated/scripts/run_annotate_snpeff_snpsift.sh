#!/bin/bash
#
#SBATCH --job-name=sm_annotate_main_job
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --mem=2000M
#SBATCH --output=slurm_logs/%x-%j.log

################################################################################
# Usage:
#   sbatch run_annotate_snpeff_snpsift.sh [SNAKEMAKE_FILE] [CONFIG_FILE] [MAX_JOBS]
#
#   - SNAKEMAKE_FILE:   Path to the Snakemake workflow file
#                       (default: "snakemake/annotate_snpeff_snpsift.smk")
#   - CONFIG_FILE:      Path to a Snakemake config file (default: "config.yaml")
#   - MAX_JOBS:         Number of Snakemake jobs (in parallel) to use (default: 50)
#
################################################################################

# ------------------------------------------------------------------------------
# 1) Parse command-line arguments with defaults
# ------------------------------------------------------------------------------
SNAKEMAKE_FILE=${1:-"snakemake/annotate_snpeff_snpsift.smk"}
CONFIG_FILE=${2:-"config.yaml"}
MAX_JOBS=${3:-50}

# ------------------------------------------------------------------------------
# 2) HPC environment setup
# ------------------------------------------------------------------------------
# First, point TMPDIR to the scratch in your home as mktemp will use this
export TMPDIR="$HOME/scratch/tmp"
# Second, create another unique temporary directory within this directory
export TMPDIR=$(mktemp -d)
# Finally, setup the cleanup trap
trap "rm -rf $TMPDIR" EXIT

# Create the slurm_logs directory if it doesn't exist
mkdir -p slurm_logs

# Export default SBATCH outputs
export SBATCH_DEFAULTS="--output=slurm_logs/%x-%j.log"

echo "DEBUG: Running annotation pipeline with:"
echo "       Snakefile:    $SNAKEMAKE_FILE"
echo "       Config file:  $CONFIG_FILE"
echo "       Max jobs:     $MAX_JOBS"
echo "       TMPDIR:       $TMPDIR"

date

# ------------------------------------------------------------------------------
# 3) Run Snakemake (e.g. via srun), using conda environment & cluster profile
# ------------------------------------------------------------------------------
srun snakemake \
    -s "$SNAKEMAKE_FILE" \
    --use-conda \
    --profile=cubi-v1 \
    -j "$MAX_JOBS" \
    --configfile "$CONFIG_FILE"

date