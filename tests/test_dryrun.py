"""Snakemake dry-run integration test."""

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _write_dryrun_fixture(tmp_path: Path, mode: str) -> Path:
    """Create temporary config and dummy input files for a Snakemake dry-run."""
    vcf_dir = tmp_path / "vcfs"
    vcf_dir.mkdir()
    src_vcf_dir = Path(REPO_ROOT) / "tests" / "data" / "vcfs"
    for name in ("SampleA.vcf.gz", "SampleB.vcf.gz"):
        dst = vcf_dir / name
        shutil.copy(src_vcf_dir / name, dst)
        (vcf_dir / f"{name}.tbi").write_bytes(b"")

    ref_fa = tmp_path / "genome.fa"
    ref_fa.write_text(">chr1\nACGT\n")
    dbnsfp = tmp_path / "dbNSFP4.9a_grch38.gz"
    dbnsfp.write_bytes(b"")
    (tmp_path / "dbNSFP4.9a_grch38.gz.tbi").write_bytes(b"")

    with open(Path(REPO_ROOT) / "tests" / "data" / "config" / "config.yaml") as fh:
        cfg = yaml.safe_load(fh)

    cfg["ref"] = {"genome": str(ref_fa), "dict": str(tmp_path / "nonexistent.dict")}
    cfg["paths"]["vcf_folder"] = str(vcf_dir)
    cfg["paths"]["output_folder"] = str(tmp_path / "out")
    cfg["snpsift"]["dbnsfp_db"] = str(dbnsfp)
    cfg["scatter"] = {
        "mode": mode,
        "count": 2,
        "canonical_contigs": ["chr1", "chrM"] if mode == "chromosome" else [],
    }

    cfg_path = tmp_path / f"config.{mode}.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))
    return cfg_path


@pytest.mark.dryrun
@pytest.mark.skipif(
    shutil.which("snakemake") is None,
    reason="snakemake not found on PATH",
)
@pytest.mark.parametrize("mode", ["none", "chromosome", "interval"])
def test_dryrun_modes(tmp_path, mode):
    """Verify the workflow DAG resolves across scatter modes."""
    cfg_path = _write_dryrun_fixture(tmp_path, mode)
    result = subprocess.run(
        [
            "snakemake",
            "-s",
            "workflow/Snakefile",
            "--configfile",
            str(cfg_path),
            "-n",
            "--cores",
            "1",
        ],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    assert result.returncode == 0, (
        f"Dry-run ({mode}) failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )


@pytest.mark.dryrun
@pytest.mark.skipif(
    shutil.which("snakemake") is None,
    reason="snakemake not found on PATH",
)
def test_dryrun_missing_tbi_fails(tmp_path):
    """Verify DAG construction fails when a configured extra annotation .vcf.gz.tbi is missing."""
    cfg_path = _write_dryrun_fixture(tmp_path, "none")
    cfg = yaml.safe_load(cfg_path.read_text())
    extra_vcf = tmp_path / "clinvar.vcf.gz"
    extra_vcf.write_bytes(b"")
    cfg["extra_annotations"] = [
        {"vcf_file": str(extra_vcf), "info_field": "CLNSIG", "annotation_prefix": "ClinVar_"}
    ]
    cfg_path.write_text(yaml.safe_dump(cfg))

    result = subprocess.run(
        [
            "snakemake",
            "-s",
            "workflow/Snakefile",
            "--configfile",
            str(cfg_path),
            "-n",
            "--cores",
            "1",
        ],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    assert result.returncode != 0


def test_scatter_rules_and_wildcard_constraints():
    """Verify scatter_unit constraints and scatter rule definitions without snakemake binary."""
    snakefile = (Path(REPO_ROOT) / "workflow" / "Snakefile").read_text()
    match = re.search(r'scatter_unit="([^"]+)"', snakefile)
    assert match is not None
    pattern = re.compile(f"^{match.group(1)}$")
    for unit in ("all", "0000-scattered", "chr1", "chrX", "chrM", "1", "MT", "other"):
        assert pattern.match(unit), f"{unit} did not match scatter_unit constraint"

    scatter_smk = (Path(REPO_ROOT) / "workflow" / "rules" / "scatter.smk").read_text()
    assert "_get_noncanonical_contigs" not in scatter_smk
    assert "-XL" not in scatter_smk
    assert "--ignore-non-starting-in-interval" in scatter_smk
    assert "-t ^" in scatter_smk
    assert "chr-prefix mismatches VCF contig" in scatter_smk
    assert scatter_smk.count("rule concatenate_annotated_vcfs:") == 1
    assert "bcftools index -n" in scatter_smk


def test_snpsift_rules_and_tracked_inputs():
    """Verify annotation_step uses -noId -tabix and tracks db_vcf/db_tbi and dbnsfp inputs."""
    snpsift_smk = (Path(REPO_ROOT) / "workflow" / "rules" / "snpsift.smk").read_text()
    assert "-noId -tabix" in snpsift_smk
    assert "dbnsfp_db=SNPSIFT_DBNSFP_DB," in snpsift_smk
    assert 'dbnsfp_tbi=f"{SNPSIFT_DBNSFP_DB}.tbi",' in snpsift_smk
    assert "-db {input.dbnsfp_db}" in snpsift_smk
    assert 'db_vcf=lambda wc: EXTRA_ANNOTATIONS[int(wc.step) - 1]["vcf_file"],' in snpsift_smk
    assert (
        "db_tbi=lambda wc: f\"{EXTRA_ANNOTATIONS[int(wc.step) - 1]['vcf_file']}.tbi\","
        in snpsift_smk
    )
    assert "{input.db_vcf} {input.vcf}" in snpsift_smk
