#!/usr/bin/env python3
"""Generate config/samples.tsv and config/config.yaml for the sm-vcf-annotation pipeline.

Scans a VCF directory, discovers samples, and auto-detects annotation databases
(reference genome, dbNSFP, extra annotation VCFs) to generate pipeline config files.

Two modes:
  Interactive:  python scripts/generate_config.py          (guided wizard)
  Flags:        python scripts/generate_config.py --vcf-folder /path/to/vcfs

Usage:
    python scripts/generate_config.py
    python scripts/generate_config.py --vcf-folder /data/vcfs
    python scripts/generate_config.py --vcf-folder /data/vcfs --dry-run
    python scripts/generate_config.py --vcf-folder /data/vcfs --ref-dir /ref --genome-build GRCh38
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Any

try:
    import pandas as pd
except ImportError:
    sys.exit(
        "Error: pandas is required but not installed.\n"
        "Install it with: pip install pandas\n"
        "(pandas is already available in Snakemake conda environments.)"
    )


# ---------------------------------------------------------------------------
# Constants & Search Path Configuration
# ---------------------------------------------------------------------------

# Genome build mapping: build -> snpEff database, dbNSFP suffix, hg alias
GENOME_BUILDS: dict[str, dict[str, str]] = {
    "GRCh38": {"snpeff_db": "GRCh38.p14", "dbnsfp_suffix": "grch38", "hg": "hg38"},
    "GRCh37": {"snpeff_db": "GRCh37.p13", "dbnsfp_suffix": "grch37", "hg": "hg19"},
}

# Reference genome file extensions (reused from sm-calling)
GENOME_EXTENSIONS = ("*.fna", "*.fa", "*.fasta")

# Search dirs for reference genome (relative to project root)
REF_SEARCH_DIRS = [
    "resources/ref/GRCh38",
    "resources/ref/GRCh37",
    "resources/ref",
    "analysis/ref/GRCh38",
    "analysis/ref/GRCh37",
    "analysis/ref",
    "../resources/ref/GRCh38",
    "../resources/ref/GRCh37",
    "../resources/ref",
]

# Search dirs for dbNSFP
DBNSFP_SEARCH_DIRS = [
    "dbnsfp",
    "resources/dbnsfp",
    "annotation/dbnsfp",
    "../dbnsfp",
    "../resources/dbnsfp",
]

# Search dirs for extra annotation VCFs
ANNOTATION_SEARCH_DIRS = [
    "annotation",
    "resources/annotation",
    "../annotation",
    "../resources/annotation",
]

# Well-known shared locations on BIH HPC
SHARED_REF_DIRS = [
    "/data/cephfs-1/work/groups/scholl/shared/ref/GRCh38",
    "/data/cephfs-1/work/groups/scholl/shared/ref/GRCh37",
    "/data/cephfs-1/work/groups/scholl/shared/ref",
]
SHARED_DBNSFP_DIRS = [
    "/data/cephfs-1/work/groups/scholl/shared/annotation/hg38/dbnsfp",
    "/data/cephfs-1/work/groups/scholl/shared/annotation/hg19/dbnsfp",
    "/data/cephfs-1/work/groups/scholl/shared/annotation/dbnsfp",
]
SHARED_ANNOTATION_DIRS = [
    "/data/cephfs-1/work/groups/scholl/shared/annotation/hg38",
    "/data/cephfs-1/work/groups/scholl/shared/annotation/hg19",
    "/data/cephfs-1/work/groups/scholl/shared/annotation",
]

# Known extra annotation databases with their default config
KNOWN_ANNOTATIONS: dict[str, dict[str, str]] = {
    "dbscSNV": {
        "glob": "**/dbscSNV*vcf.gz",
        "info_field": "dbscSNV_ada_score,dbscSNV_rf_score",
        "annotation_prefix": "splice_",
        "description": "dbscSNV splice-site predictions",
    },
    "spidex": {
        "glob": "**/spidex*vcf.gz",
        "info_field": "spidex_dpsi_max_tissue,spidex_dpsi_zscore",
        "annotation_prefix": "splice_",
        "description": "SPIDEX splice-site predictions",
    },
    "clinvar": {
        "glob": "**/clinvar*vcf.gz",
        "info_field": "CLNSIG",
        "annotation_prefix": "ClinVar_",
        "description": "ClinVar clinical significance",
    },
    "hgmd": {
        "glob": "**/HGMD*vcf.gz",
        "info_field": "CLASS",
        "annotation_prefix": "hgmd_",
        "description": "HGMD variant classification",
    },
    "cosmic": {
        "glob": "**/Cosmic*vcf.gz",
        "info_field": "GENOME_SCREEN_SAMPLE_COUNT",
        "annotation_prefix": "COSMIC_",
        "description": "COSMIC genome screen counts",
    },
}

# Default canonical contigs
DEFAULT_CANONICAL_CONTIGS_CHR = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
DEFAULT_CANONICAL_CONTIGS_NOCHR = [str(i) for i in range(1, 23)] + ["X", "Y"]

# dbNSFP version-specific field lists
# v4.x fields (dbNSFP 4.x series)
DBNSFP_FIELDS_V4 = ",".join(
    [
        "aaref",
        "aaalt",
        "aapos",
        "SIFT_pred",
        "SIFT4G_score",
        "Polyphen2_HDIV_pred",
        "Polyphen2_HVAR_score",
        "LRT_pred",
        "MutationTaster_score",
        "MutationTaster_converted_rankscore",
        "MutationTaster_pred",
        "MutationAssessor_score",
        "MutationAssessor_rankscore",
        "MutationAssessor_pred",
        "FATHMM_score",
        "FATHMM_converted_rankscore",
        "FATHMM_pred",
        "PROVEAN_score",
        "PROVEAN_converted_rankscore",
        "PROVEAN_pred",
        "MetaSVM_score",
        "MetaSVM_rankscore",
        "MetaSVM_pred",
        "M-CAP_score",
        "M-CAP_rankscore",
        "M-CAP_pred",
        "REVEL_score",
        "REVEL_rankscore",
        "MVP_score",
        "MVP_rankscore",
        "MPC_score",
        "MPC_rankscore",
        "CADD_raw_rankscore",
        "CADD_phred",
        "GERP++_NR",
        "GERP++_RS",
        "GERP++_RS_rankscore",
        "phastCons100way_vertebrate",
        "1000Gp3_AC",
        "1000Gp3_AF",
        "1000Gp3_AFR_AC",
        "1000Gp3_AFR_AF",
        "1000Gp3_EUR_AC",
        "1000Gp3_EUR_AF",
        "1000Gp3_AMR_AC",
        "1000Gp3_AMR_AF",
        "1000Gp3_EAS_AC",
        "1000Gp3_EAS_AF",
        "1000Gp3_SAS_AC",
        "1000Gp3_SAS_AF",
        "gnomAD_exomes_flag",
        "gnomAD_exomes_AC",
        "gnomAD_exomes_AN",
        "gnomAD_exomes_AF",
        "gnomAD_exomes_nhomalt",
        "gnomAD_genomes_flag",
        "gnomAD_genomes_AC",
        "gnomAD_genomes_AN",
        "gnomAD_genomes_AF",
        "gnomAD_genomes_nhomalt",
        "ALFA_Asian_AF",
        "ALFA_Total_AC",
        "ALFA_Total_AN",
        "ALFA_Total_AF",
        "clinvar_id",
        "clinvar_clnsig",
    ]
)

# v5.x fields (dbNSFP 5.x series)
# Changes from v4: LRT retired; FATHMM retired -> fathmm-XF; MutationTaster
# rankscore renamed; gnomAD exomes/genomes merged into gnomAD4.1 joint
DBNSFP_FIELDS_V5 = ",".join(
    [
        "aaref",
        "aaalt",
        "aapos",
        "SIFT_pred",
        "SIFT4G_score",
        "Polyphen2_HDIV_pred",
        "Polyphen2_HVAR_score",
        "MutationTaster_score",
        "MutationTaster_rankscore",
        "MutationTaster_pred",
        "MutationAssessor_score",
        "MutationAssessor_rankscore",
        "MutationAssessor_pred",
        "fathmm-XF_coding_score",
        "fathmm-XF_coding_rankscore",
        "fathmm-XF_coding_pred",
        "PROVEAN_score",
        "PROVEAN_converted_rankscore",
        "PROVEAN_pred",
        "MetaSVM_score",
        "MetaSVM_rankscore",
        "MetaSVM_pred",
        "M-CAP_score",
        "M-CAP_rankscore",
        "M-CAP_pred",
        "REVEL_score",
        "REVEL_rankscore",
        "MVP_score",
        "MVP_rankscore",
        "MPC_score",
        "MPC_rankscore",
        "CADD_raw_rankscore",
        "CADD_phred",
        "GERP++_NR",
        "GERP++_RS",
        "GERP++_RS_rankscore",
        "phastCons100way_vertebrate",
        "1000Gp3_AC",
        "1000Gp3_AF",
        "1000Gp3_AFR_AC",
        "1000Gp3_AFR_AF",
        "1000Gp3_EUR_AC",
        "1000Gp3_EUR_AF",
        "1000Gp3_AMR_AC",
        "1000Gp3_AMR_AF",
        "1000Gp3_EAS_AC",
        "1000Gp3_EAS_AF",
        "1000Gp3_SAS_AC",
        "1000Gp3_SAS_AF",
        "gnomAD4.1_joint_flag",
        "gnomAD4.1_joint_AC",
        "gnomAD4.1_joint_AN",
        "gnomAD4.1_joint_AF",
        "gnomAD4.1_joint_nhomalt",
        "ALFA_Asian_AF",
        "ALFA_Total_AC",
        "ALFA_Total_AN",
        "ALFA_Total_AF",
        "clinvar_id",
        "clinvar_clnsig",
    ]
)


# ---------------------------------------------------------------------------
# Section A: VCF Discovery
# ---------------------------------------------------------------------------


def discover_vcf_files(vcf_folder: str | Path) -> list[dict]:
    """Scan a directory for VCF files and return metadata.

    Returns a list of dicts with keys: path, basename, size_bytes.
    """
    folder = Path(vcf_folder)
    if not folder.is_dir():
        return []

    entries = []
    for vcf in sorted(folder.glob("*.vcf.gz")):
        if vcf.is_file():
            entries.append(
                {
                    "path": str(vcf),
                    "basename": vcf.name.replace(".vcf.gz", ""),
                    "size_bytes": vcf.stat().st_size,
                }
            )
    return entries


def detect_dbnsfp_version(dbnsfp_path: str | None) -> int:
    """Detect the major dbNSFP version from a file path.

    Returns 4 or 5 based on the filename. Defaults to 4 if undetectable.
    """
    if dbnsfp_path is None:
        return 4
    match = re.search(r"dbNSFP(\d+)", Path(dbnsfp_path).name, re.IGNORECASE)
    if match:
        return int(match.group(1))
    return 4


def get_dbnsfp_fields(version: int) -> str:
    """Return the appropriate dbNSFP field list for the given major version."""
    if version >= 5:
        return DBNSFP_FIELDS_V5
    return DBNSFP_FIELDS_V4


def _infer_output_folder(vcf_folder: str) -> str:
    """Infer the annotation output folder from the VCF input folder.

    Walks up the path to find a 'variant_calls' component and replaces
    it (and everything below) with 'annotation'. Falls back to replacing
    the last path component.

    Examples:
        ../results/A5297/variant_calls/mutect2/filtered_vcfs
            -> ../results/A5297/annotation
        ../results/exomes/variant_calls/filtered_vcfs
            -> ../results/exomes/annotation
        results/final
            -> results/annotation
    """
    parts = Path(vcf_folder).parts
    # Find the 'variant_calls' component
    for i, part in enumerate(parts):
        if part == "variant_calls":
            return str(Path(*parts[:i]) / "annotation").replace("\\", "/")
    # Fallback: replace last component
    return str(Path(vcf_folder).parent / "annotation").replace("\\", "/")


def format_size(size_bytes: int) -> str:
    """Convert bytes to human-readable format."""
    value = float(size_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024:
            return f"{value:.1f} {unit}"
        value /= 1024
    return f"{value:.1f} PB"


# ---------------------------------------------------------------------------
# Section B: Reference & Annotation Discovery
# ---------------------------------------------------------------------------


def detect_genome_build(genome_name: str) -> str:
    """Detect genome build from a filename or path string.

    Returns "GRCh38" or "GRCh37" based on heuristics. Defaults to "GRCh38".
    """
    lower = genome_name.lower()
    if any(tag in lower for tag in ("grch37", "hg19", "hs37", "b37")):
        return "GRCh37"
    return "GRCh38"


def _find_genome_fastas(search_dir: Path) -> list[dict[str, Any]]:
    """Find FASTA genome files in a directory.

    Returns list of dicts with keys: path, has_fai, name.
    """
    results: list[dict[str, Any]] = []
    if not search_dir.is_dir():
        return results

    for ext in GENOME_EXTENSIONS:
        for fasta in sorted(search_dir.glob(ext)):
            if fasta.is_file():
                fai = Path(str(fasta) + ".fai")
                results.append(
                    {
                        "path": str(fasta),
                        "has_fai": fai.is_file(),
                        "name": fasta.name,
                    }
                )
    return results


def _find_dict_file(genome_path: str) -> str | None:
    """Find a .dict file for a given genome FASTA path.

    Tries: genome.dict (replacing extension), genome.fa.dict (appended), and
    sibling .dict files in the same directory.
    """
    gp = Path(genome_path)

    # Try replacing extension: genome.fa -> genome.dict
    dict_path = gp.with_suffix(".dict")
    if dict_path.is_file():
        return str(dict_path)

    # Try appending: genome.fa.dict
    appended = Path(str(gp) + ".dict")
    if appended.is_file():
        return str(appended)

    # Try any .dict in the same directory
    parent = gp.parent
    if parent.is_dir():
        dicts = sorted(parent.glob("*.dict"))
        if len(dicts) == 1:
            return str(dicts[0])

    return None


def discover_reference_data(
    ref_dir: Path | None = None,
    project_root: Path | None = None,
) -> dict[str, Any]:
    """Discover reference genome FASTA, dict, and build.

    Searches ref_dir (if given), then candidate directories relative to
    project_root, then shared HPC paths.

    Returns dict with keys: genome, dict, build, search_log.
    """
    search_log: list[str] = []
    if project_root is None:
        project_root = Path.cwd()

    # Build search list
    search_dirs: list[Path] = []
    if ref_dir is not None:
        search_dirs.append(Path(ref_dir))
    for rel in REF_SEARCH_DIRS:
        search_dirs.append(project_root / rel)
    for shared in SHARED_REF_DIRS:
        search_dirs.append(Path(shared))

    # Search each directory
    for sd in search_dirs:
        fastas = _find_genome_fastas(sd)
        if fastas:
            # Prefer indexed genomes
            indexed = [f for f in fastas if f["has_fai"]]
            best = indexed[0] if indexed else fastas[0]
            genome_path = best["path"]
            build = detect_genome_build(genome_path)
            dict_path = _find_dict_file(genome_path)
            search_log.append(f"Found genome: {genome_path} (build={build})")
            return {
                "genome": genome_path,
                "dict": dict_path,
                "build": build,
                "search_log": search_log,
            }
        if sd.is_dir():
            search_log.append(f"No FASTA in: {sd}")
        else:
            search_log.append(f"Dir not found: {sd}")

    search_log.append("No reference genome found in any search directory.")
    return {"genome": None, "dict": None, "build": "GRCh38", "search_log": search_log}


def discover_dbnsfp(
    dbnsfp_dir: Path | None = None,
    project_root: Path | None = None,
    build: str = "GRCh38",
) -> dict[str, Any]:
    """Discover dbNSFP database file matching the genome build.

    Returns dict with keys: path, search_log.
    """
    search_log: list[str] = []
    if project_root is None:
        project_root = Path.cwd()

    build_info = GENOME_BUILDS.get(build, GENOME_BUILDS["GRCh38"])
    suffix = build_info["dbnsfp_suffix"]

    search_dirs: list[Path] = []
    if dbnsfp_dir is not None:
        search_dirs.append(Path(dbnsfp_dir))
    for rel in DBNSFP_SEARCH_DIRS:
        search_dirs.append(project_root / rel)
    for shared in SHARED_DBNSFP_DIRS:
        search_dirs.append(Path(shared))

    for sd in search_dirs:
        if not sd.is_dir():
            search_log.append(f"Dir not found: {sd}")
            continue

        # Look for build-matching dbNSFP files
        matches = sorted(sd.glob(f"dbNSFP*{suffix}.gz"))
        if not matches:
            # Try case-insensitive via broader glob
            matches = [
                f
                for f in sorted(sd.glob("dbNSFP*.gz"))
                if suffix in f.name.lower() and not f.name.endswith(".tbi.gz")
            ]

        if matches:
            best = str(matches[0])
            search_log.append(f"Found dbNSFP: {best}")
            return {"path": best, "search_log": search_log}

        search_log.append(f"No dbNSFP for {build} in: {sd}")

    search_log.append("No dbNSFP database found.")
    return {"path": None, "search_log": search_log}


def discover_extra_annotations(
    annotation_dir: Path | None = None,
    project_root: Path | None = None,
    build: str = "GRCh38",
) -> dict[str, Any]:
    """Discover known extra annotation VCFs for the given genome build.

    Returns dict with keys: annotations (list of dicts), search_log.
    """
    search_log: list[str] = []
    if project_root is None:
        project_root = Path.cwd()

    build_info = GENOME_BUILDS.get(build, GENOME_BUILDS["GRCh38"])
    hg = build_info["hg"]  # e.g. "hg38"

    search_dirs: list[Path] = []
    if annotation_dir is not None:
        search_dirs.append(Path(annotation_dir))
    for rel in ANNOTATION_SEARCH_DIRS:
        search_dirs.append(project_root / rel)
    for shared in SHARED_ANNOTATION_DIRS:
        search_dirs.append(Path(shared))

    found: dict[str, dict[str, str]] = {}

    for sd in search_dirs:
        if not sd.is_dir():
            continue

        for name, info in KNOWN_ANNOTATIONS.items():
            if name in found:
                continue

            matches = sorted(sd.glob(info["glob"]))

            # Filter to build-matching paths if possible
            build_matches = [
                m
                for m in matches
                if hg in str(m).lower()
                or build.lower() in str(m).lower()
                or (hg not in str(sd).lower() and build.lower() not in str(sd).lower())
            ]
            candidates = build_matches if build_matches else matches

            if candidates:
                vcf_path = str(candidates[0])
                found[name] = {
                    "name": name,
                    "vcf_file": vcf_path,
                    "info_field": info["info_field"],
                    "annotation_prefix": info["annotation_prefix"],
                }
                search_log.append(f"Found {name}: {vcf_path}")

    if not found:
        search_log.append("No extra annotation VCFs found.")

    return {"annotations": list(found.values()), "search_log": search_log}


def parse_canonical_contigs(dict_path: str | None, build: str = "GRCh38") -> list[str]:
    """Parse canonical contigs from a .dict file.

    Extracts @SQ SN: lines and filters to chr1-22 + X + Y (or 1-22 + X + Y).
    Falls back to defaults if dict_path is None or unreadable.
    """
    if dict_path is None:
        if build == "GRCh37":
            return list(DEFAULT_CANONICAL_CONTIGS_NOCHR)
        return list(DEFAULT_CANONICAL_CONTIGS_CHR)

    try:
        contigs: list[str] = []
        with open(dict_path) as fh:
            for line in fh:
                if line.startswith("@SQ"):
                    match = re.search(r"SN:(\S+)", line)
                    if match:
                        contigs.append(match.group(1))
    except OSError:
        if build == "GRCh37":
            return list(DEFAULT_CANONICAL_CONTIGS_NOCHR)
        return list(DEFAULT_CANONICAL_CONTIGS_CHR)

    if not contigs:
        if build == "GRCh37":
            return list(DEFAULT_CANONICAL_CONTIGS_NOCHR)
        return list(DEFAULT_CANONICAL_CONTIGS_CHR)

    # Detect chr prefix from first contig
    has_chr = contigs[0].startswith("chr")

    # Filter to canonical autosomes + sex chromosomes
    if has_chr:
        canonical = set(DEFAULT_CANONICAL_CONTIGS_CHR)
    else:
        canonical = set(DEFAULT_CANONICAL_CONTIGS_NOCHR)

    return [c for c in contigs if c in canonical]


def discover_sibling_config(
    project_root: Path | None = None,
) -> dict[str, Any]:
    """Look for a sibling sm-calling project config and extract useful paths.

    Returns dict with keys: vcf_folder, ref_genome, found, search_log.
    """
    search_log: list[str] = []
    if project_root is None:
        project_root = Path.cwd()

    sibling = project_root / ".." / "sm-calling" / "config" / "config.yaml"
    sibling = sibling.resolve()

    result: dict[str, Any] = {
        "vcf_folder": None,
        "ref_genome": None,
        "found": False,
        "search_log": search_log,
    }

    if not sibling.is_file():
        search_log.append(f"No sibling config at: {sibling}")
        return result

    search_log.append(f"Found sibling config: {sibling}")
    result["found"] = True

    try:
        text = sibling.read_text()
    except OSError:
        search_log.append("Could not read sibling config.")
        return result

    # Extract ref genome path (simple regex, no YAML dependency)
    genome_match = re.search(r'^\s*genome:\s*["\']?([^"\'#\n]+)', text, re.MULTILINE)
    if genome_match:
        result["ref_genome"] = genome_match.group(1).strip()
        search_log.append(f"  ref.genome: {result['ref_genome']}")

    # Extract output_folder to infer VCF location
    output_match = re.search(r'^\s*output_folder:\s*["\']?([^"\'#\n]+)', text, re.MULTILINE)
    if output_match:
        output_folder = output_match.group(1).strip()
        # Common pattern: sm-calling output has filtered VCFs
        sibling_root = sibling.parent.parent
        vcf_candidates = [
            sibling_root / output_folder / "mutect2" / "filtered_vcfs",
            sibling_root / output_folder / "freebayes" / "filtered_vcfs",
            sibling_root / output_folder / "filtered_vcfs",
        ]
        for candidate in vcf_candidates:
            if candidate.is_dir():
                result["vcf_folder"] = str(candidate)
                search_log.append(f"  vcf_folder: {result['vcf_folder']}")
                break
        if result["vcf_folder"] is None:
            search_log.append(f"  output_folder: {output_folder} (no filtered VCFs found)")

    return result


# ---------------------------------------------------------------------------
# Section C: Sample & Config Generation
# ---------------------------------------------------------------------------


def build_samples_dataframe(vcf_entries: list[dict]) -> pd.DataFrame:
    """Create a samples DataFrame from discovered VCF entries."""
    rows = []
    for entry in vcf_entries:
        rows.append(
            {
                "sample": entry["basename"],
                "vcf_basename": entry["basename"],
            }
        )
    return pd.DataFrame(rows)


def write_samples_tsv(
    df: pd.DataFrame,
    output_path: str | Path = "config/samples.tsv",
    *,
    dry_run: bool = False,
    overwrite: bool = False,
) -> str:
    """Write samples DataFrame to TSV.

    Returns the content that was (or would be) written.
    """
    content: str = df.to_csv(sep="\t", index=False)
    path = Path(output_path)

    if dry_run:
        return content

    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists. Use --overwrite or delete it first.")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return content


def generate_config_template(
    vcf_folder: str = "results/final",
    output_folder: str | None = None,
    output_path: str | Path = "config/config.yaml",
    ref_data: dict[str, Any] | None = None,
    dbnsfp_data: dict[str, Any] | None = None,
    extra_annotations: list[dict[str, str]] | None = None,
    scatter_config: dict[str, Any] | None = None,
    *,
    dry_run: bool = False,
    overwrite: bool = False,
) -> str:
    """Generate a config.yaml template with discovered paths.

    Returns the YAML content that was (or would be) written.
    """
    # Resolve reference data
    if ref_data and ref_data.get("genome"):
        genome_path = ref_data["genome"]
        dict_path = ref_data.get("dict") or "EDIT_ME: /path/to/reference/genome.dict"
        build = ref_data.get("build", "GRCh38")
    else:
        genome_path = "EDIT_ME: /path/to/reference/genome.fa"
        dict_path = "EDIT_ME: /path/to/reference/genome.dict"
        build = "GRCh38"

    build_info = GENOME_BUILDS.get(build, GENOME_BUILDS["GRCh38"])

    # Resolve dbNSFP
    if dbnsfp_data and dbnsfp_data.get("path"):
        dbnsfp_path = dbnsfp_data["path"]
    else:
        dbnsfp_path = f"EDIT_ME: /path/to/dbNSFP4.9a_{build_info['dbnsfp_suffix']}.gz"

    # Resolve dbNSFP fields based on detected version
    dbnsfp_version = detect_dbnsfp_version(dbnsfp_path)
    dbnsfp_fields = get_dbnsfp_fields(dbnsfp_version)

    # Resolve output folder
    resolved_output_folder = output_folder or _infer_output_folder(vcf_folder)

    # Resolve scatter config
    if scatter_config:
        scatter_mode = scatter_config.get("mode", "none")
        scatter_count = scatter_config.get("count", 100)
        contigs = scatter_config.get("canonical_contigs", [])
    else:
        scatter_mode = "none"
        scatter_count = 100
        contigs = []

    # Format contigs as YAML list
    if contigs:
        contigs_yaml = "\n" + "\n".join(f'    - "{c}"' for c in contigs)
    else:
        contigs_yaml = " []"

    # Format extra annotations as YAML
    if extra_annotations:
        ann_lines = []
        for ann in extra_annotations:
            ann_lines.append(f'  - vcf_file: "{ann["vcf_file"]}"')
            ann_lines.append(f'    info_field: "{ann["info_field"]}"')
            ann_lines.append(f'    annotation_prefix: "{ann["annotation_prefix"]}"')
        extra_ann_yaml = "\n" + "\n".join(ann_lines)
    else:
        extra_ann_yaml = " []"

    content = f"""\
# =============================================================================
# sm-vcf-annotation configuration
# =============================================================================
# Generated by: python scripts/generate_config.py

ref:
  genome: "{genome_path}"
  dict: "{dict_path}"

paths:
  samples: "config/samples.tsv"
  vcf_folder: "{vcf_folder}"
  output_folder: "{resolved_output_folder}"
  log_subdir: "logs"
  annotation_subdir: "annotation"

snpeff:
  database: "{build_info["snpeff_db"]}"
  extra_flags: "-lof -noInteraction -noMotif -noNextProt -spliceRegionIntronMax 12"

snpsift:
  dbnsfp_db: "{dbnsfp_path}"
  dbnsfp_fields: "{dbnsfp_fields}"

scatter:
  mode: "{scatter_mode}"
  count: {scatter_count}
  canonical_contigs:{contigs_yaml}

extra_annotations:{extra_ann_yaml}
"""

    path = Path(output_path)

    if dry_run:
        return content

    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists. Use --overwrite or delete it first.")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return content


# ---------------------------------------------------------------------------
# Section D: Prompt Helpers
# ---------------------------------------------------------------------------


def _prompt(prompt: str, default: str = "") -> str:
    """Prompt the user for input with an optional default value."""
    if default:
        result = input(f"{prompt} [{default}]: ").strip()
        return result if result else default
    return input(f"{prompt}: ").strip()


def _prompt_yn(prompt: str, default: bool = True) -> bool:
    """Prompt for yes/no with a default."""
    suffix = "[Y/n]" if default else "[y/N]"
    result = input(f"{prompt} {suffix}: ").strip().lower()
    if not result:
        return default
    return result in ("y", "yes")


def _prompt_path(prompt: str, default: str = "", *, must_exist: bool = True) -> str:
    """Prompt for a filesystem path, re-prompting on invalid input."""
    while True:
        raw = _prompt(prompt, default)
        if not raw:
            if not must_exist:
                return ""
            print("  Path cannot be empty. Please try again.")
            continue
        p = Path(raw).expanduser()
        if must_exist and not p.exists():
            print(f"  Path does not exist: {p}")
            retry = _prompt_yn("  Try again?", default=True)
            if not retry:
                return str(p)
            continue
        return str(p)


def _prompt_choice(prompt: str, choices: list[str], default: str = "") -> str:
    """Prompt with numbered choices."""
    for i, choice in enumerate(choices, 1):
        marker = " <-- current" if choice == default else ""
        print(f"    [{i}] {choice}{marker}")
    while True:
        raw = _prompt(f"  {prompt}", default=default)
        if raw.isdigit():
            idx = int(raw) - 1
            if 0 <= idx < len(choices):
                return choices[idx]
        if raw in choices:
            return raw
        print(f"  Invalid choice. Enter 1-{len(choices)} or a value.")


# ---------------------------------------------------------------------------
# Section E: Interactive Mode
# ---------------------------------------------------------------------------


def _interactive_wizard() -> None:
    """Run the interactive config generation wizard."""
    print("=" * 60)
    print("  sm-vcf-annotation -- Config Generator")
    print("=" * 60)
    print()

    project_root = Path.cwd()

    # --- Step 1: Check sibling project ---
    sibling = discover_sibling_config(project_root)
    default_vcf_folder = "results/final"
    if sibling["found"]:
        print("Found sibling sm-calling project config.")
        if sibling.get("vcf_folder"):
            default_vcf_folder = sibling["vcf_folder"]
            print(f"  Suggested VCF folder: {default_vcf_folder}")
        print()

    # --- Step 2: VCF folder ---
    vcf_folder = _prompt_path("Path to VCF folder", default=default_vcf_folder, must_exist=False)

    # --- Step 3: Discover VCFs ---
    entries = discover_vcf_files(vcf_folder)
    if not entries:
        print(f"\nNo VCF files found in: {vcf_folder}")
        print("Continuing anyway -- you can add samples manually to samples.tsv.")
        print()
    else:
        print(f"\nFound {len(entries)} VCF file(s):")
        for e in entries:
            print(f"  {e['basename']}.vcf.gz  ({format_size(e['size_bytes'])})")
        print()

    # --- Step 3b: Output folder ---
    default_output = _infer_output_folder(vcf_folder)
    output_folder = _prompt("Output folder", default=default_output)

    # --- Step 4: Reference genome ---
    print("Scanning for reference genome...")
    ref_data = discover_reference_data(project_root=project_root)
    if ref_data.get("genome"):
        print(f"  Found: {ref_data['genome']}  (build: {ref_data['build']})")
        if not _prompt_yn("  Use this reference?", default=True):
            ref_path = _prompt_path("  Path to reference genome FASTA", must_exist=False)
            ref_data = {
                "genome": ref_path,
                "dict": _find_dict_file(ref_path),
                "build": detect_genome_build(ref_path),
                "search_log": [],
            }
    else:
        print("  No reference genome found automatically.")
        ref_path = _prompt_path(
            "  Path to reference genome FASTA (or press Enter to skip)",
            must_exist=False,
        )
        if ref_path:
            ref_data = {
                "genome": ref_path,
                "dict": _find_dict_file(ref_path),
                "build": detect_genome_build(ref_path),
                "search_log": [],
            }
    print()

    # --- Step 5: Genome build ---
    build = ref_data.get("build", "GRCh38")
    build = _prompt_choice("Genome build", choices=["GRCh38", "GRCh37"], default=build)
    ref_data["build"] = build
    print()

    # --- Step 6: dbNSFP ---
    print("Scanning for dbNSFP database...")
    dbnsfp_data = discover_dbnsfp(project_root=project_root, build=build)
    if dbnsfp_data.get("path"):
        print(f"  Found: {dbnsfp_data['path']}")
        if not _prompt_yn("  Use this dbNSFP?", default=True):
            db_path = _prompt_path("  Path to dbNSFP file", must_exist=False)
            dbnsfp_data = {"path": db_path, "search_log": []}
    else:
        print("  No dbNSFP database found automatically.")
        db_path = _prompt_path("  Path to dbNSFP file (or press Enter to skip)", must_exist=False)
        if db_path:
            dbnsfp_data = {"path": db_path, "search_log": []}
    print()

    # --- Step 7: Extra annotations ---
    print("Scanning for extra annotation VCFs...")
    ann_data = discover_extra_annotations(project_root=project_root, build=build)
    extra_annotations: list[dict[str, str]] = []
    if ann_data["annotations"]:
        print(f"  Found {len(ann_data['annotations'])} annotation database(s):")
        for ann in ann_data["annotations"]:
            desc = KNOWN_ANNOTATIONS.get(ann["name"], {}).get("description", ann["name"])
            print(f"    - {ann['name']}: {desc}")
            print(f"      {ann['vcf_file']}")
        if _prompt_yn("  Include all found annotations?", default=True):
            extra_annotations = ann_data["annotations"]
        else:
            for ann in ann_data["annotations"]:
                if _prompt_yn(f"  Include {ann['name']}?", default=True):
                    extra_annotations.append(ann)
    else:
        print("  No extra annotation VCFs found.")
    print()

    # --- Step 8: Scatter mode ---
    scatter_mode = _prompt_choice(
        "Scatter mode", choices=["none", "chromosome", "interval"], default="none"
    )
    scatter_config: dict[str, Any] = {"mode": scatter_mode, "count": 100, "canonical_contigs": []}

    if scatter_mode in ("chromosome", "interval"):
        if scatter_mode == "interval":
            count_str = _prompt("  Number of scatter intervals", default="100")
            try:
                scatter_config["count"] = int(count_str)
            except ValueError:
                print("  Invalid number, using default 100.")

        # Parse canonical contigs from .dict or use defaults
        contigs = parse_canonical_contigs(ref_data.get("dict"), build)
        print(
            f"  Using {len(contigs)} canonical contigs (from {'dict file' if ref_data.get('dict') else 'defaults'})."
        )
        scatter_config["canonical_contigs"] = contigs
    print()

    # --- Step 9: Confirm & write ---
    print("=" * 60)
    print("  Summary")
    print("=" * 60)
    print(f"  VCF folder:     {vcf_folder}")
    print(f"  Output folder:  {output_folder}")
    print(f"  Samples:        {len(entries)}")
    print(f"  Reference:      {ref_data.get('genome', 'not set')}")
    print(f"  Build:          {build}")
    print(f"  dbNSFP:         {dbnsfp_data.get('path', 'not set')}")
    print(f"  Annotations:    {len(extra_annotations)}")
    print(f"  Scatter mode:   {scatter_mode}")
    print()

    if not _prompt_yn("Write config files?", default=True):
        print("Aborted.")
        return

    # Write samples.tsv
    if entries:
        df = build_samples_dataframe(entries)
        try:
            write_samples_tsv(df, overwrite=True)
            print("  Wrote: config/samples.tsv")
        except Exception as e:
            print(f"  Error writing samples.tsv: {e}")

    # Write config.yaml
    try:
        generate_config_template(
            vcf_folder=vcf_folder,
            output_folder=output_folder,
            ref_data=ref_data,
            dbnsfp_data=dbnsfp_data,
            extra_annotations=extra_annotations if extra_annotations else None,
            scatter_config=scatter_config,
            overwrite=True,
        )
        print("  Wrote: config/config.yaml")
    except Exception as e:
        print(f"  Error writing config.yaml: {e}")

    # Final message
    has_edit_me = not ref_data.get("genome") or not dbnsfp_data.get("path")
    if has_edit_me:
        print("\nDone! Search for EDIT_ME in config/config.yaml to fill in missing paths.")
    else:
        print("\nDone! Review config/config.yaml before running the pipeline.")


# ---------------------------------------------------------------------------
# Section F: CLI Entry Point
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate config files for sm-vcf-annotation.",
    )
    parser.add_argument(
        "--vcf-folder",
        type=str,
        default=None,
        help="Path to directory containing input VCF files.",
    )
    parser.add_argument(
        "--ref-dir",
        type=str,
        default=None,
        help="Reference genome directory to search.",
    )
    parser.add_argument(
        "--dbnsfp-dir",
        type=str,
        default=None,
        help="dbNSFP database directory to search.",
    )
    parser.add_argument(
        "--annotation-dir",
        type=str,
        default=None,
        help="Extra annotations directory to search.",
    )
    parser.add_argument(
        "--genome-build",
        type=str,
        default=None,
        choices=["GRCh37", "GRCh38"],
        help="Genome build (default: auto-detect from reference).",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Pipeline output directory (default: inferred from --vcf-folder).",
    )
    parser.add_argument(
        "--scatter-mode",
        type=str,
        default="none",
        choices=["none", "chromosome", "interval"],
        help="Scatter mode (default: none).",
    )
    parser.add_argument(
        "--scatter-count",
        type=int,
        default=100,
        help="Number of scatter intervals (default: 100).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print output without writing files.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing config files.",
    )

    args = parser.parse_args()

    # Interactive mode if no VCF folder specified
    if args.vcf_folder is None:
        _interactive_wizard()
        return

    # --- CLI mode ---
    project_root = Path.cwd()

    # Discover VCFs
    entries = discover_vcf_files(args.vcf_folder)
    if not entries:
        print(f"No VCF files found in: {args.vcf_folder}", file=sys.stderr)
        sys.exit(1)

    # Resolve output folder
    output_folder = args.output_dir or _infer_output_folder(args.vcf_folder)

    print(f"Found {len(entries)} VCF file(s) in {args.vcf_folder}")
    print(f"Output folder: {output_folder}")
    for e in entries:
        print(f"  {e['basename']}.vcf.gz  ({format_size(e['size_bytes'])})")

    # Discover reference data
    ref_dir = Path(args.ref_dir) if args.ref_dir else None
    ref_data = discover_reference_data(ref_dir=ref_dir, project_root=project_root)
    build = args.genome_build or ref_data.get("build", "GRCh38")
    ref_data["build"] = build

    if ref_data.get("genome"):
        print(f"Reference: {ref_data['genome']}  (build: {build})")
    else:
        print(f"No reference genome found (build: {build})")

    # Discover dbNSFP
    dbnsfp_dir = Path(args.dbnsfp_dir) if args.dbnsfp_dir else None
    dbnsfp_data = discover_dbnsfp(dbnsfp_dir=dbnsfp_dir, project_root=project_root, build=build)
    if dbnsfp_data.get("path"):
        print(f"dbNSFP: {dbnsfp_data['path']}")

    # Discover extra annotations
    annotation_dir = Path(args.annotation_dir) if args.annotation_dir else None
    ann_data = discover_extra_annotations(
        annotation_dir=annotation_dir, project_root=project_root, build=build
    )
    extra_annotations = ann_data["annotations"] if ann_data["annotations"] else None
    if extra_annotations:
        print(f"Extra annotations: {len(extra_annotations)} found")
        for ann in extra_annotations:
            print(f"  - {ann['name']}: {ann['vcf_file']}")

    # Build scatter config
    scatter_config: dict[str, Any] = {
        "mode": args.scatter_mode,
        "count": args.scatter_count,
        "canonical_contigs": [],
    }
    if args.scatter_mode in ("chromosome", "interval"):
        scatter_config["canonical_contigs"] = parse_canonical_contigs(ref_data.get("dict"), build)

    # Build samples
    df = build_samples_dataframe(entries)

    if args.dry_run:
        print("\n--- samples.tsv ---")
        print(write_samples_tsv(df, dry_run=True))
        print("--- config.yaml ---")
        print(
            generate_config_template(
                vcf_folder=args.vcf_folder,
                output_folder=output_folder,
                ref_data=ref_data,
                dbnsfp_data=dbnsfp_data,
                extra_annotations=extra_annotations,
                scatter_config=scatter_config,
                dry_run=True,
            )
        )
        return

    write_samples_tsv(df, overwrite=args.overwrite)
    print("Wrote: config/samples.tsv")

    generate_config_template(
        vcf_folder=args.vcf_folder,
        output_folder=output_folder,
        ref_data=ref_data,
        dbnsfp_data=dbnsfp_data,
        extra_annotations=extra_annotations,
        scatter_config=scatter_config,
        overwrite=args.overwrite,
    )
    print("Wrote: config/config.yaml")

    has_edit_me = not ref_data.get("genome") or not dbnsfp_data.get("path")
    if has_edit_me:
        print("\nDone! Search for EDIT_ME in config/config.yaml to fill in missing paths.")
    else:
        print("\nDone! Review config/config.yaml before running the pipeline.")


if __name__ == "__main__":
    main()
