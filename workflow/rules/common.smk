"""Shared configuration shortcuts, metadata accessors, and helper wrappers."""

import os

import pandas as pd

from rules.helpers import (
    get_java_opts as _get_java_opts_impl,
    get_samples as _get_samples_impl,
    get_scatter_units as _get_scatter_units_impl,
    get_vcf_path as _get_vcf_path_impl,
    annotation_step_vcf as _annotation_step_vcf_impl,
)


# -- Config shortcuts --------------------------------------------------------
REF = config.get("ref", {}).get("genome", "")
REF_DICT = config.get("ref", {}).get("dict", "")

VCF_FOLDER = config["paths"]["vcf_folder"]
OUTPUT_DIR = config["paths"]["output_folder"]
LOG_SUBDIR = config["paths"].get("log_subdir", "logs")
ANNOTATION_SUBDIR = config["paths"].get("annotation_subdir", "annotation")

ANNOTATION_DIR = os.path.join(OUTPUT_DIR, ANNOTATION_SUBDIR)
LOG_DIR = os.path.join(OUTPUT_DIR, LOG_SUBDIR)
INTERVALS_DIR = os.path.join(OUTPUT_DIR, "intervals")

SNPEFF_DB = config["snpeff"]["database"]
SNPEFF_EXTRA_FLAGS = config["snpeff"].get("extra_flags", "")

SNPSIFT_DBNSFP_DB = config["snpsift"]["dbnsfp_db"]
SNPSIFT_DBNSFP_FIELDS = config["snpsift"]["dbnsfp_fields"]

SCATTER_MODE = config["scatter"]["mode"]
SCATTER_COUNT = config["scatter"].get("count", 100)
CANONICAL_CONTIGS = config["scatter"].get("canonical_contigs", [])

EXTRA_ANNOTATIONS = config.get("extra_annotations", [])
N_ANNOTATIONS = len(EXTRA_ANNOTATIONS)

SCRATCH_DIR = os.environ.get("TMPDIR", "/tmp")


# -- Scatter units -----------------------------------------------------------
SCATTER_UNITS = _get_scatter_units_impl(SCATTER_MODE, SCATTER_COUNT)


# -- Sample accessors --------------------------------------------------------
def get_sample_list():
    return _get_samples_impl(samples_df)


# -- Input helper functions --------------------------------------------------
def get_input_vcf(wildcards):
    """Return the input VCF path for a sample."""
    return _get_vcf_path_impl(wildcards.sample, VCF_FOLDER, samples_df)


def get_java_opts(wildcards, resources):
    """Derive Java options from allocated resources."""
    return _get_java_opts_impl(resources.mem_mb, resources.tmpdir)


def get_annotation_step_input(wildcards):
    """Return input path for a step-wise extra annotation."""
    step = int(wildcards.step)
    return _annotation_step_vcf_impl(
        ANNOTATION_DIR, wildcards.sample, step - 1, wildcards.scatter_unit
    )


def get_final_annotation_input(wildcards):
    """Return input path for the finalize step (last annotation step or dbNSFP)."""
    return _annotation_step_vcf_impl(
        ANNOTATION_DIR, wildcards.sample, N_ANNOTATIONS, wildcards.scatter_unit
    )


# -- Final output dispatcher -------------------------------------------------
def get_final_outputs():
    """Return all expected final output files."""
    samples = get_sample_list()
    if SCATTER_MODE == "none":
        return expand(
            os.path.join(ANNOTATION_DIR, "{sample}.annotated.vcf.gz"),
            sample=samples,
        )
    else:
        # Scatter mode: final output is the concatenated VCF per sample
        return expand(
            os.path.join(ANNOTATION_DIR, "{sample}.annotated.vcf.gz"),
            sample=samples,
        )
