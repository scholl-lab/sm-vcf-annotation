"""Snakemake dry-run integration test."""

import os
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
        idx_ext = "csi" if mode == "chromosome" and name == "SampleB.vcf.gz" else "tbi"
        (vcf_dir / f"{name}.{idx_ext}").write_bytes(b"")

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
    """Verify the workflow DAG resolves across scatter modes and renders expected commands."""
    cfg_path = _write_dryrun_fixture(tmp_path, mode)
    result = subprocess.run(
        [
            "snakemake",
            "-s",
            "workflow/Snakefile",
            "--configfile",
            str(cfg_path),
            "-n",
            "-p",
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
    if mode == "interval":
        assert "--variant-output-filtering STARTS_IN" in result.stdout
        assert "--ignore-non-starting-in-interval" not in result.stdout
        assert "bcftools index -n" in result.stdout
        assert "record count mismatch" in result.stdout
    elif mode == "chromosome":
        assert "-t ^chr1,chrM" in result.stdout
        assert "bcftools index -n" in result.stdout
        assert "record count mismatch" in result.stdout


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
    combined = result.stdout + result.stderr
    assert "MissingInputException" in combined
    assert "clinvar.vcf.gz.tbi" in combined


@pytest.mark.dryrun
@pytest.mark.skipif(
    shutil.which("snakemake") is None,
    reason="snakemake not found on PATH",
)
def test_dryrun_extra_annotations_and_canon_list(tmp_path):
    """Verify positive dry-run with extra_annotations (.tbi present) and snpeff.canon_list."""
    cfg_path = _write_dryrun_fixture(tmp_path, "none")
    cfg = yaml.safe_load(cfg_path.read_text())
    extra_vcf = tmp_path / "clinvar.vcf.gz"
    extra_vcf.write_bytes(b"")
    (tmp_path / "clinvar.vcf.gz.tbi").write_bytes(b"")
    canon_list = tmp_path / "mane_canon.txt"
    canon_list.write_text("NM_000720.4\nNM_001128840.3\n")
    cfg["snpeff"]["canon_list"] = str(canon_list)
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
            "-p",
            "--cores",
            "1",
        ],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    assert result.returncode == 0, (
        f"Dry-run failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )
    assert "-noId -tabix" in result.stdout
    assert f"-onlyTr {canon_list}" in result.stdout


@pytest.mark.skipif(
    shutil.which("bcftools") is None,
    reason="bcftools not found on PATH",
)
def test_chromosome_scatter_other_shard_and_count_check(tmp_path):
    """Verify -r and -t ^ split captures chrUn/pCMV in 'other' and record count check works."""
    raw_vcf = tmp_path / "raw.vcf"
    raw_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        "##contig=<ID=chrM,length=1000>\n"
        "##contig=<ID=chrUn_KI270442v1,length=1000>\n"
        "##contig=<ID=pCMV,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\t.\tA\tG\t.\tPASS\t.\n"
        "chrM\t20\t.\tC\tT\t.\tPASS\t.\n"
        "chrUn_KI270442v1\t30\t.\tG\tA\t.\tPASS\t.\n"
        "pCMV\t40\t.\tT\tC\t.\tPASS\t.\n"
    )
    raw_gz = tmp_path / "raw.vcf.gz"
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(raw_gz), str(raw_vcf)], check=True)
    subprocess.run(["bcftools", "index", "-t", str(raw_gz)], check=True)

    chr1_gz = tmp_path / "chr1.vcf.gz"
    chrm_gz = tmp_path / "chrM.vcf.gz"
    other_gz = tmp_path / "other.vcf.gz"
    subprocess.run(
        ["bcftools", "view", "-r", "chr1", str(raw_gz), "-Oz", "-o", str(chr1_gz)], check=True
    )
    subprocess.run(
        ["bcftools", "view", "-r", "chrM", str(raw_gz), "-Oz", "-o", str(chrm_gz)], check=True
    )
    subprocess.run(
        ["bcftools", "view", "-t", "^chr1,chrM", str(raw_gz), "-Oz", "-o", str(other_gz)],
        check=True,
    )

    other_records = subprocess.run(
        ["bcftools", "view", "-H", str(other_gz)], capture_output=True, text=True, check=True
    ).stdout.splitlines()
    assert [line.split("\t")[0] for line in other_records] == ["chrUn_KI270442v1", "pCMV"]

    merged_gz = tmp_path / "merged.vcf.gz"
    subprocess.run(
        [
            "bcftools",
            "concat",
            "-Oz",
            str(chr1_gz),
            str(chrm_gz),
            str(other_gz),
            "-o",
            str(merged_gz),
        ],
        check=True,
    )
    subprocess.run(["bcftools", "index", "-t", str(merged_gz)], check=True)

    check_cmd = (
        f'[ "$(bcftools view -H {raw_gz} | wc -l)" -eq "$(bcftools index -n {merged_gz})" ] '
        '|| { echo "record count mismatch"; exit 1; }'
    )
    assert subprocess.run(["bash", "-c", check_cmd]).returncode == 0

    incomplete_gz = tmp_path / "incomplete.vcf.gz"
    subprocess.run(
        ["bcftools", "concat", "-Oz", str(chr1_gz), str(chrm_gz), "-o", str(incomplete_gz)],
        check=True,
    )
    subprocess.run(["bcftools", "index", "-t", str(incomplete_gz)], check=True)
    bad_cmd = (
        f'[ "$(bcftools view -H {raw_gz} | wc -l)" -eq "$(bcftools index -n {incomplete_gz})" ] '
        '|| { echo "record count mismatch"; exit 1; }'
    )
    res = subprocess.run(["bash", "-c", bad_cmd], capture_output=True, text=True)
    assert res.returncode != 0
    assert "record count mismatch" in res.stdout
