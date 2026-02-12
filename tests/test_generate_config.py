"""Tests for scripts/generate_config.py."""

import os
import sys

import pytest

# Add scripts to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from generate_config import (
    build_samples_dataframe,
    discover_vcf_files,
    format_size,
    generate_config_template,
    write_samples_tsv,
)


class TestDiscoverVcfFiles:
    def test_finds_vcfs(self, tmp_path):
        (tmp_path / "sample1.vcf.gz").write_bytes(b"\x1f\x8b" + b"\x00" * 10)
        (tmp_path / "sample2.vcf.gz").write_bytes(b"\x1f\x8b" + b"\x00" * 20)
        (tmp_path / "not_a_vcf.txt").write_text("hello")

        entries = discover_vcf_files(tmp_path)
        assert len(entries) == 2
        basenames = [e["basename"] for e in entries]
        assert "sample1" in basenames
        assert "sample2" in basenames

    def test_empty_directory(self, tmp_path):
        entries = discover_vcf_files(tmp_path)
        assert entries == []

    def test_nonexistent_directory(self):
        entries = discover_vcf_files("/nonexistent/path")
        assert entries == []

    def test_size_captured(self, tmp_path):
        data = b"\x1f\x8b" + b"\x00" * 100
        (tmp_path / "test.vcf.gz").write_bytes(data)
        entries = discover_vcf_files(tmp_path)
        assert entries[0]["size_bytes"] == len(data)


class TestFormatSize:
    def test_bytes(self):
        assert format_size(500) == "500.0 B"

    def test_kilobytes(self):
        assert format_size(1024) == "1.0 KB"

    def test_megabytes(self):
        assert format_size(1024 * 1024) == "1.0 MB"

    def test_gigabytes(self):
        assert format_size(1024**3) == "1.0 GB"


class TestBuildSamplesDataframe:
    def test_basic(self):
        entries = [
            {"path": "/a/s1.vcf.gz", "basename": "s1", "size_bytes": 100},
            {"path": "/a/s2.vcf.gz", "basename": "s2", "size_bytes": 200},
        ]
        df = build_samples_dataframe(entries)
        assert list(df.columns) == ["sample", "vcf_basename"]
        assert len(df) == 2
        assert df.iloc[0]["sample"] == "s1"
        assert df.iloc[1]["vcf_basename"] == "s2"

    def test_empty(self):
        df = build_samples_dataframe([])
        assert len(df) == 0


class TestWriteSamplesTsv:
    def test_dry_run(self):
        import pandas as pd

        df = pd.DataFrame({"sample": ["A"], "vcf_basename": ["A"]})
        content = write_samples_tsv(df, dry_run=True)
        assert "sample" in content
        assert "A" in content

    def test_writes_file(self, tmp_path):
        import pandas as pd

        df = pd.DataFrame({"sample": ["A"], "vcf_basename": ["A"]})
        out = tmp_path / "samples.tsv"
        write_samples_tsv(df, output_path=out)
        assert out.exists()
        assert "sample" in out.read_text()

    def test_no_overwrite(self, tmp_path):
        import pandas as pd

        df = pd.DataFrame({"sample": ["A"], "vcf_basename": ["A"]})
        out = tmp_path / "samples.tsv"
        out.write_text("existing")
        with pytest.raises(FileExistsError):
            write_samples_tsv(df, output_path=out)


class TestGenerateConfigTemplate:
    def test_dry_run(self):
        content = generate_config_template(dry_run=True)
        assert "snpeff:" in content
        assert "snpsift:" in content
        assert "scatter:" in content

    def test_vcf_folder_in_output(self):
        content = generate_config_template(vcf_folder="/my/vcfs", dry_run=True)
        assert "/my/vcfs" in content

    def test_writes_file(self, tmp_path):
        out = tmp_path / "config.yaml"
        generate_config_template(output_path=out)
        assert out.exists()
        assert "snpeff:" in out.read_text()
