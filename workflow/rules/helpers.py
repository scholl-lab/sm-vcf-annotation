"""Pure helper functions for the sm-vcf-annotation workflow.

All functions are free of Snakemake imports so they can be unit-tested.
"""

import os
import re
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

    - mode="chromosome": returns the chromosome list plus "other" (e.g. ["chr1", ..., "chrM", "other"])
    - mode="interval": returns ["0000-scattered", "0001-scattered", ...]
    - mode="none": returns ["all"]
    """
    if mode == "chromosome":
        if not chromosomes:
            raise ValueError(
                "get_scatter_units(mode='chromosome') requires a non-empty 'chromosomes' list."
            )
        return [*chromosomes, "other"]
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

    SnpSift 5.2 produces comma-separated rows like::

        Sample ,SampleA,Total
        Transitions ,150488,150488
        Transversions ,70878,70878
        Ts/Tv ,2.123,2.123

    Colon-separated output is also supported::

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
        if "," in line:
            label, rest = line.split(",", 1)
            values = [value.strip() for value in rest.split(",")]
        elif ":" in line:
            label, rest = line.split(":", 1)
            values = rest.strip().split()
        else:
            continue
        label = label.strip()
        if label == "Sample":
            continue
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


AGGREGATORS = {
    "max": lambda v: max(v, key=float),
    "min": lambda v: min(v, key=float),
    "first": lambda v: v[0],
}


def normalize_info_header(line: str, new_ids: set[str], collapse: dict[str, str]) -> str:
    """Fix the Number of an INFO header line written by SnpSift annotate.

    SnpSift copies the database's Number (often 1) but joins the values of every matching
    database record ("1,1,1"). Collapsed fields become Number=1; other new fields Number=.,
    unless they are per-allele (A, R) or per-genotype (G).
    """
    if not line.startswith("##INFO=<ID="):
        return line
    fid = line[len("##INFO=<ID=") :].split(",", 1)[0]
    if fid in collapse:
        return re.sub(r",Number=[^,]+,", ",Number=1,", line, count=1)
    if fid in new_ids:
        return re.sub(r",Number=(?![ARG.],)[^,]+,", ",Number=.,", line, count=1)
    return line


def collapse_info(info: str, collapse: dict[str, str]) -> str:
    """Collapse joined multi-record values ("1,1,3") of the given INFO keys to one value."""
    if not collapse or info == ".":
        return info
    out = []
    for kv in info.split(";"):
        key, sep, val = kv.partition("=")
        if sep and key in collapse and "," in val:
            vals = [x for x in val.split(",") if x not in ("", ".")]
            val = AGGREGATORS[collapse[key]](vals) if vals else "."
        out.append(f"{key}{sep}{val}")
    return ";".join(out)


def normalize_annotation_vcf(lines, new_ids: list[str], collapse: dict[str, str]):
    """Stream a VCF (text lines) and apply normalize_info_header / collapse_info."""
    ids = set(new_ids)
    for line in lines:
        if line.startswith("#"):
            yield normalize_info_header(line, ids, collapse)
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) > 7:
            f[7] = collapse_info(f[7], collapse)
        yield "\t".join(f) + "\n"
