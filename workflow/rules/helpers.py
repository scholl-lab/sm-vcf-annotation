"""Pure helper functions for the sm-vcf-annotation workflow.

All functions are free of Snakemake imports so they can be unit-tested.
"""

import os

import pandas as pd


def get_java_opts(mem_mb: int, tmpdir: str) -> str:
    """Derive Java options from allocated resources (80% heap / 20% stack)."""
    xmx = int(mem_mb * 0.8)
    xms = int(mem_mb * 0.2)
    return f"-Xms{xms}m -Xmx{xmx}m -Djava.io.tmpdir={tmpdir}"


def get_samples(samples_df: pd.DataFrame) -> list[str]:
    """Return sorted list of unique sample identifiers."""
    return sorted(samples_df["sample"].unique().tolist())


def get_scatter_units(mode: str, scatter_count: int = 100) -> list[str]:
    """Return the list of scatter units based on mode.

    - mode="interval": returns ["0000-scattered", "0001-scattered", ...]
    - mode="none": returns ["all"]
    """
    if mode == "interval":
        return [f"{i:04d}-scattered" for i in range(scatter_count)]
    return ["all"]


def get_vcf_path(sample: str, vcf_folder: str, samples_df: pd.DataFrame) -> str:
    """Resolve the input VCF path for a sample."""
    row = samples_df.loc[samples_df["sample"] == sample].iloc[0]
    basename = row["vcf_basename"]
    return os.path.join(vcf_folder, f"{basename}.vcf.gz")


def annotation_step_vcf(annotation_dir: str, sample: str, step: int, scatter_unit: str) -> str:
    """Return the path for a given step of extra annotation.

    step=0 means the VCF after dbNSFP; step=1..N are intermediate files.
    """
    if scatter_unit == "all":
        unit_part = ""
    else:
        unit_part = f".{scatter_unit}"

    if step == 0:
        return os.path.join(annotation_dir, f"{sample}{unit_part}.ann.dbnsfp.vcf.gz")
    return os.path.join(annotation_dir, f"{sample}{unit_part}.ann.step{step}.vcf.gz")


def format_extra_annotations(annotations: list[dict]) -> str:
    """Format extra annotations into a space-separated string."""
    parts = []
    for ann in annotations:
        parts.append(f"{ann['vcf_file']},{ann['info_field']},{ann['annotation_prefix']}")
    return " ".join(parts)
