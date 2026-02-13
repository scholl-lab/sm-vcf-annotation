"""Pure helper functions for the sm-vcf-annotation workflow.

All functions are free of Snakemake imports so they can be unit-tested.
"""

import os
from collections.abc import Iterable

import pandas as pd


def get_java_opts(mem_mb: int, tmpdir: str) -> str:
    """Derive Java options from allocated resources (80% max heap / 20% initial heap)."""
    xmx = int(mem_mb * 0.8)
    xms = int(mem_mb * 0.2)
    return f"-Xms{xms}m -Xmx{xmx}m -Djava.io.tmpdir={tmpdir}"


def get_samples(samples_df: pd.DataFrame) -> list[str]:
    """Return sorted list of unique sample identifiers."""
    return sorted(samples_df["sample"].unique().tolist())


def get_scatter_units(
    mode: str,
    chromosomes: list[str] | None = None,
    scatter_count: int = 100,
) -> list[str]:
    """Return the list of scatter units based on mode.

    - mode="chromosome": returns the chromosome list (e.g. ["chr1", ..., "chrY"])
    - mode="interval": returns ["0000-scattered", "0001-scattered", ...]
    - mode="none": returns ["all"]
    """
    if mode == "chromosome":
        if not chromosomes:
            raise ValueError(
                "get_scatter_units(mode='chromosome') requires a non-empty 'chromosomes' list."
            )
        return chromosomes
    if mode == "interval":
        return [f"{i:04d}-scattered" for i in range(scatter_count)]
    return ["all"]


def get_vcf_path(sample: str, vcf_folder: str, samples_df: pd.DataFrame) -> str:
    """Resolve the input VCF path for a sample."""
    row = samples_df.loc[samples_df["sample"] == sample].iloc[0]
    basename = row["vcf_basename"]
    return os.path.join(vcf_folder, f"{basename}.vcf.gz")


def parse_annotation_completeness(
    lines: Iterable[str],
    field_names: list[str],
) -> dict[str, dict[str, int | float]]:
    """Parse bcftools query output and compute annotation completeness per field.

    Parameters
    ----------
    lines:
        Iterable of TSV lines from ``bcftools query``.  Each line has
        CHROM, POS, then one column per field in *field_names*.  Missing
        values are represented by ``"."``.  Accepts a file-like object
        (for streaming from ``Popen.stdout``) or a list of strings.
    field_names:
        Column names corresponding to positions 2..N in the TSV.

    Returns
    -------
    dict mapping each field name to
        ``{"total": int, "annotated": int, "rate": float}``.
    """
    stats: dict[str, dict[str, int | float]] = {
        f: {"total": 0, "annotated": 0, "rate": 0.0} for f in field_names
    }
    for line in lines:
        line = line.rstrip("\n")
        if not line:
            continue
        cols = line.split("\t")
        # First two columns are CHROM, POS
        values = cols[2:]
        for i, field in enumerate(field_names):
            if i < len(values):
                stats[field]["total"] = int(stats[field]["total"]) + 1
                val = values[i]
                if val not in (".", ""):
                    stats[field]["annotated"] = int(stats[field]["annotated"]) + 1
    for field in field_names:
        total = stats[field]["total"]
        if total > 0:
            stats[field]["rate"] = round(int(stats[field]["annotated"]) / int(total), 4)
    return stats


def parse_snpsift_tstv(tstv_output: str) -> dict[str, str]:
    """Parse SnpSift tstv output into a dict of metric → value.

    SnpSift tstv produces rows like::

        Sample        : 1       Total
        Transitions   : 150488  150488
        Transversions : 70878   70878
        Ts/Tv         : 2.123   2.123

    Returns the last column value for each metric row, keyed by row label
    (e.g. ``{"Transitions": "150488", "Transversions": "70878", "Ts/Tv": "2.123"}``).
    When multiple samples are present the last column is *Total*; for a
    single-sample VCF it is the only value column.
    """
    result: dict[str, str] = {}
    for line in tstv_output.strip().splitlines():
        if ":" not in line:
            continue
        label, rest = line.split(":", 1)
        label = label.strip()
        if label == "Sample":
            continue
        values = rest.strip().split()
        # Last column is Total when multiple samples; only column otherwise
        result[label] = values[-1] if values else ""
    return result


def annotation_step_vcf(annotation_dir: str, sample: str, step: int, scatter_unit: str) -> str:
    """Return the path for a given step of extra annotation.

    step=0 means the VCF after dbNSFP; step=1..N are intermediate files.
    The scatter_unit is always included in the filename to match rule output patterns.
    """
    if step == 0:
        return os.path.join(annotation_dir, f"{sample}.{scatter_unit}.ann.dbnsfp.vcf.gz")
    return os.path.join(annotation_dir, f"{sample}.{scatter_unit}.ann.step{step}.vcf.gz")
