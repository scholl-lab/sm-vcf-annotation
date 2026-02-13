"""Tests for scripts/generate_config.py."""

import os
import sys
import textwrap

import pytest

# Add scripts to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from generate_config import (
    DEFAULT_CANONICAL_CONTIGS_CHR,
    DEFAULT_CANONICAL_CONTIGS_NOCHR,
    _infer_output_folder,
    build_samples_dataframe,
    detect_genome_build,
    discover_dbnsfp,
    discover_extra_annotations,
    discover_reference_data,
    discover_sibling_config,
    discover_vcf_files,
    format_size,
    generate_config_template,
    parse_canonical_contigs,
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


class TestInferOutputFolder:
    def test_variant_calls_mutect2(self):
        result = _infer_output_folder("../results/A5297/variant_calls/mutect2/filtered_vcfs")
        assert result == "../results/A5297/annotation"

    def test_variant_calls_freebayes(self):
        result = _infer_output_folder("../results/exomes/variant_calls/freebayes/filtered_vcfs")
        assert result == "../results/exomes/annotation"

    def test_variant_calls_direct(self):
        result = _infer_output_folder("../results/cohort/variant_calls/filtered_vcfs")
        assert result == "../results/cohort/annotation"

    def test_variant_calls_top_level(self):
        result = _infer_output_folder("../results/cohort/variant_calls")
        assert result == "../results/cohort/annotation"

    def test_fallback_no_variant_calls(self):
        result = _infer_output_folder("results/final")
        assert result == "results/annotation"

    def test_fallback_simple_path(self):
        result = _infer_output_folder("/data/vcfs")
        assert result == "/data/annotation"

    def test_forward_slashes_on_windows(self):
        # Ensure no backslashes in output
        result = _infer_output_folder("../results/A5297/variant_calls/mutect2/filtered")
        assert "\\" not in result


class TestDetectGenomeBuild:
    def test_grch38_default(self):
        assert detect_genome_build("genome.fa") == "GRCh38"

    def test_grch38_explicit(self):
        assert detect_genome_build("GRCh38.p14.fna") == "GRCh38"

    def test_grch37_from_name(self):
        assert detect_genome_build("GRCh37.p13.fna") == "GRCh37"

    def test_hg19(self):
        assert detect_genome_build("hg19.fa") == "GRCh37"

    def test_hs37(self):
        assert detect_genome_build("hs37d5.fa") == "GRCh37"

    def test_b37(self):
        assert detect_genome_build("b37.fa") == "GRCh37"

    def test_grch37_in_path(self):
        assert detect_genome_build("/ref/GRCh37/genome.fa") == "GRCh37"

    def test_case_insensitive(self):
        assert detect_genome_build("GRCH37.fa") == "GRCh37"
        assert detect_genome_build("HG19.fa") == "GRCh37"


class TestDiscoverReferenceData:
    def test_finds_genome_in_explicit_dir(self, tmp_path):
        ref_dir = tmp_path / "ref"
        ref_dir.mkdir()
        genome = ref_dir / "GRCh38.p14.fna"
        genome.write_text(">chr1\nACGT\n")
        fai = ref_dir / "GRCh38.p14.fna.fai"
        fai.write_text("chr1\t4\t6\t4\t5\n")

        result = discover_reference_data(ref_dir=ref_dir, project_root=tmp_path)
        assert result["genome"] == str(genome)
        assert result["build"] == "GRCh38"

    def test_finds_genome_in_search_dir(self, tmp_path):
        # Create resources/ref/GRCh38 relative to project root
        ref_dir = tmp_path / "resources" / "ref" / "GRCh38"
        ref_dir.mkdir(parents=True)
        genome = ref_dir / "genome.fa"
        genome.write_text(">chr1\nACGT\n")

        result = discover_reference_data(project_root=tmp_path)
        assert result["genome"] == str(genome)

    def test_prefers_indexed_genome(self, tmp_path):
        ref_dir = tmp_path / "ref"
        ref_dir.mkdir()
        # Unindexed
        (ref_dir / "aaa.fa").write_text(">chr1\nACGT\n")
        # Indexed
        indexed = ref_dir / "zzz.fa"
        indexed.write_text(">chr1\nACGT\n")
        (ref_dir / "zzz.fa.fai").write_text("chr1\t4\n")

        result = discover_reference_data(ref_dir=ref_dir, project_root=tmp_path)
        assert result["genome"] == str(indexed)

    def test_finds_dict_file(self, tmp_path):
        ref_dir = tmp_path / "ref"
        ref_dir.mkdir()
        genome = ref_dir / "genome.fa"
        genome.write_text(">chr1\nACGT\n")
        dict_file = ref_dir / "genome.dict"
        dict_file.write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n")

        result = discover_reference_data(ref_dir=ref_dir, project_root=tmp_path)
        assert result["dict"] == str(dict_file)

    def test_no_genome_found(self, tmp_path):
        result = discover_reference_data(project_root=tmp_path)
        assert result["genome"] is None
        assert result["dict"] is None
        assert result["build"] == "GRCh38"

    def test_search_log_populated(self, tmp_path):
        result = discover_reference_data(project_root=tmp_path)
        assert len(result["search_log"]) > 0

    def test_detects_grch37_from_path(self, tmp_path):
        ref_dir = tmp_path / "resources" / "ref" / "GRCh37"
        ref_dir.mkdir(parents=True)
        genome = ref_dir / "hs37d5.fa"
        genome.write_text(">1\nACGT\n")

        result = discover_reference_data(project_root=tmp_path)
        assert result["build"] == "GRCh37"


class TestDiscoverDbnsfp:
    def test_finds_dbnsfp_matching_build(self, tmp_path):
        db_dir = tmp_path / "dbnsfp"
        db_dir.mkdir()
        db_file = db_dir / "dbNSFP4.9a_grch38.gz"
        db_file.write_bytes(b"\x1f\x8b" + b"\x00" * 10)

        result = discover_dbnsfp(dbnsfp_dir=db_dir, project_root=tmp_path, build="GRCh38")
        assert result["path"] == str(db_file)

    def test_finds_dbnsfp_in_search_dirs(self, tmp_path):
        db_dir = tmp_path / "resources" / "dbnsfp"
        db_dir.mkdir(parents=True)
        db_file = db_dir / "dbNSFP4.9a_grch37.gz"
        db_file.write_bytes(b"\x1f\x8b" + b"\x00" * 10)

        result = discover_dbnsfp(project_root=tmp_path, build="GRCh37")
        assert result["path"] == str(db_file)

    def test_no_match_wrong_build(self, tmp_path):
        db_dir = tmp_path / "dbnsfp"
        db_dir.mkdir()
        # grch37 file but searching for grch38
        (db_dir / "dbNSFP4.9a_grch37.gz").write_bytes(b"\x1f\x8b" + b"\x00" * 10)

        result = discover_dbnsfp(dbnsfp_dir=db_dir, project_root=tmp_path, build="GRCh38")
        assert result["path"] is None

    def test_no_dbnsfp_found(self, tmp_path):
        result = discover_dbnsfp(project_root=tmp_path, build="GRCh38")
        assert result["path"] is None
        assert len(result["search_log"]) > 0

    def test_search_log_tracks_dirs(self, tmp_path):
        result = discover_dbnsfp(project_root=tmp_path, build="GRCh38")
        assert any("not found" in msg or "No dbNSFP" in msg for msg in result["search_log"])


class TestDiscoverExtraAnnotations:
    def test_finds_known_annotations(self, tmp_path):
        ann_dir = tmp_path / "annotation"
        ann_dir.mkdir()
        clinvar = ann_dir / "clinvar_20230101.vcf.gz"
        clinvar.write_bytes(b"\x1f\x8b" + b"\x00" * 10)

        result = discover_extra_annotations(
            annotation_dir=ann_dir, project_root=tmp_path, build="GRCh38"
        )
        found_names = [a["name"] for a in result["annotations"]]
        assert "clinvar" in found_names

    def test_finds_multiple_annotations(self, tmp_path):
        ann_dir = tmp_path / "annotation"
        ann_dir.mkdir()
        (ann_dir / "clinvar.vcf.gz").write_bytes(b"\x1f\x8b" + b"\x00" * 10)
        (ann_dir / "HGMD_Pro.vcf.gz").write_bytes(b"\x1f\x8b" + b"\x00" * 10)

        result = discover_extra_annotations(
            annotation_dir=ann_dir, project_root=tmp_path, build="GRCh38"
        )
        found_names = {a["name"] for a in result["annotations"]}
        assert "clinvar" in found_names
        assert "hgmd" in found_names

    def test_no_annotations_found(self, tmp_path):
        result = discover_extra_annotations(project_root=tmp_path, build="GRCh38")
        assert result["annotations"] == []

    def test_annotations_have_required_keys(self, tmp_path):
        ann_dir = tmp_path / "annotation"
        ann_dir.mkdir()
        (ann_dir / "clinvar.vcf.gz").write_bytes(b"\x1f\x8b" + b"\x00" * 10)

        result = discover_extra_annotations(
            annotation_dir=ann_dir, project_root=tmp_path, build="GRCh38"
        )
        for ann in result["annotations"]:
            assert "name" in ann
            assert "vcf_file" in ann
            assert "info_field" in ann
            assert "annotation_prefix" in ann

    def test_finds_in_subdirectories(self, tmp_path):
        ann_dir = tmp_path / "annotation" / "hg38"
        ann_dir.mkdir(parents=True)
        (ann_dir / "dbscSNV1.1_hg38.vcf.gz").write_bytes(b"\x1f\x8b" + b"\x00" * 10)

        result = discover_extra_annotations(
            annotation_dir=tmp_path / "annotation",
            project_root=tmp_path,
            build="GRCh38",
        )
        found_names = [a["name"] for a in result["annotations"]]
        assert "dbscSNV" in found_names


class TestParseCanonicalContigs:
    def test_parse_dict_with_chr_prefix(self, tmp_path):
        dict_file = tmp_path / "genome.dict"
        lines = ["@HD\tVN:1.6\n"]
        for c in DEFAULT_CANONICAL_CONTIGS_CHR:
            lines.append(f"@SQ\tSN:{c}\tLN:1000\n")
        # Add non-canonical
        lines.append("@SQ\tSN:chrM\tLN:16569\n")
        lines.append("@SQ\tSN:chr1_random\tLN:500\n")
        dict_file.write_text("".join(lines))

        contigs = parse_canonical_contigs(str(dict_file))
        assert contigs == DEFAULT_CANONICAL_CONTIGS_CHR
        assert "chrM" not in contigs

    def test_parse_dict_without_chr_prefix(self, tmp_path):
        dict_file = tmp_path / "genome.dict"
        lines = ["@HD\tVN:1.6\n"]
        for c in DEFAULT_CANONICAL_CONTIGS_NOCHR:
            lines.append(f"@SQ\tSN:{c}\tLN:1000\n")
        lines.append("@SQ\tSN:MT\tLN:16569\n")
        dict_file.write_text("".join(lines))

        contigs = parse_canonical_contigs(str(dict_file), build="GRCh37")
        assert contigs == DEFAULT_CANONICAL_CONTIGS_NOCHR
        assert "MT" not in contigs

    def test_default_contigs_grch38(self):
        contigs = parse_canonical_contigs(None, build="GRCh38")
        assert contigs == DEFAULT_CANONICAL_CONTIGS_CHR

    def test_default_contigs_grch37(self):
        contigs = parse_canonical_contigs(None, build="GRCh37")
        assert contigs == DEFAULT_CANONICAL_CONTIGS_NOCHR

    def test_missing_file_returns_defaults(self):
        contigs = parse_canonical_contigs("/nonexistent/genome.dict", build="GRCh38")
        assert contigs == DEFAULT_CANONICAL_CONTIGS_CHR

    def test_empty_dict_returns_defaults(self, tmp_path):
        dict_file = tmp_path / "empty.dict"
        dict_file.write_text("")
        contigs = parse_canonical_contigs(str(dict_file), build="GRCh38")
        assert contigs == DEFAULT_CANONICAL_CONTIGS_CHR

    def test_preserves_order_from_dict(self, tmp_path):
        dict_file = tmp_path / "genome.dict"
        # Write in reverse order
        lines = ["@HD\tVN:1.6\n"]
        for c in reversed(DEFAULT_CANONICAL_CONTIGS_CHR):
            lines.append(f"@SQ\tSN:{c}\tLN:1000\n")
        dict_file.write_text("".join(lines))

        contigs = parse_canonical_contigs(str(dict_file))
        # Order should match dict file (reversed)
        assert contigs == list(reversed(DEFAULT_CANONICAL_CONTIGS_CHR))


class TestDiscoverSiblingConfig:
    def test_no_sibling(self, tmp_path):
        result = discover_sibling_config(project_root=tmp_path)
        assert result["found"] is False
        assert result["vcf_folder"] is None

    def test_finds_sibling_config(self, tmp_path):
        # Create sibling sm-calling config
        sibling_dir = tmp_path / "sm-calling" / "config"
        sibling_dir.mkdir(parents=True)
        config = sibling_dir / "config.yaml"
        config.write_text(
            textwrap.dedent("""\
            ref:
              genome: "/ref/GRCh38.fa"
            paths:
              output_folder: "results/variant_calls"
            """)
        )

        # Project root is a sibling directory
        project_root = tmp_path / "sm-vcf-annotation"
        project_root.mkdir()

        result = discover_sibling_config(project_root=project_root)
        assert result["found"] is True
        assert result["ref_genome"] == "/ref/GRCh38.fa"

    def test_finds_vcf_folder_from_sibling(self, tmp_path):
        # Create sibling sm-calling with output directories
        sibling_root = tmp_path / "sm-calling"
        sibling_config_dir = sibling_root / "config"
        sibling_config_dir.mkdir(parents=True)
        (sibling_config_dir / "config.yaml").write_text(
            textwrap.dedent("""\
            paths:
              output_folder: "results/variant_calls"
            """)
        )
        # Create the filtered VCFs directory
        vcf_dir = sibling_root / "results" / "variant_calls" / "mutect2" / "filtered_vcfs"
        vcf_dir.mkdir(parents=True)

        project_root = tmp_path / "sm-vcf-annotation"
        project_root.mkdir()

        result = discover_sibling_config(project_root=project_root)
        assert result["found"] is True
        assert result["vcf_folder"] is not None
        assert "filtered_vcfs" in result["vcf_folder"]


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

    def test_no_overwrite(self, tmp_path):
        out = tmp_path / "config.yaml"
        out.write_text("existing")
        with pytest.raises(FileExistsError):
            generate_config_template(output_path=out)

    def test_ref_data_populates_genome(self):
        ref_data = {
            "genome": "/ref/GRCh38.fa",
            "dict": "/ref/GRCh38.dict",
            "build": "GRCh38",
        }
        content = generate_config_template(ref_data=ref_data, dry_run=True)
        assert "/ref/GRCh38.fa" in content
        assert "/ref/GRCh38.dict" in content
        assert "GRCh38.p14" in content

    def test_grch37_build_sets_snpeff_db(self):
        ref_data = {
            "genome": "/ref/GRCh37.fa",
            "dict": "/ref/GRCh37.dict",
            "build": "GRCh37",
        }
        content = generate_config_template(ref_data=ref_data, dry_run=True)
        assert "GRCh37.p13" in content

    def test_dbnsfp_data_populates_path(self):
        dbnsfp_data = {"path": "/db/dbNSFP4.9a_grch38.gz"}
        content = generate_config_template(dbnsfp_data=dbnsfp_data, dry_run=True)
        assert "/db/dbNSFP4.9a_grch38.gz" in content

    def test_extra_annotations_in_output(self):
        extra = [
            {
                "name": "clinvar",
                "vcf_file": "/ann/clinvar.vcf.gz",
                "info_field": "CLNSIG",
                "annotation_prefix": "ClinVar_",
            }
        ]
        content = generate_config_template(extra_annotations=extra, dry_run=True)
        assert "/ann/clinvar.vcf.gz" in content
        assert "CLNSIG" in content
        assert "ClinVar_" in content

    def test_scatter_config_in_output(self):
        scatter = {
            "mode": "interval",
            "count": 50,
            "canonical_contigs": ["chr1", "chr2", "chr3"],
        }
        content = generate_config_template(scatter_config=scatter, dry_run=True)
        assert 'mode: "interval"' in content
        assert "count: 50" in content
        assert '"chr1"' in content
        assert '"chr2"' in content
        assert '"chr3"' in content

    def test_edit_me_placeholders_when_no_data(self):
        content = generate_config_template(dry_run=True)
        assert "EDIT_ME:" in content

    def test_no_edit_me_when_data_provided(self):
        ref_data = {
            "genome": "/ref/genome.fa",
            "dict": "/ref/genome.dict",
            "build": "GRCh38",
        }
        dbnsfp_data = {"path": "/db/dbNSFP.gz"}
        content = generate_config_template(
            ref_data=ref_data, dbnsfp_data=dbnsfp_data, dry_run=True
        )
        assert "EDIT_ME:" not in content

    def test_output_folder_inferred_from_vcf_folder(self):
        content = generate_config_template(
            vcf_folder="../results/A5297/variant_calls/mutect2/filtered",
            dry_run=True,
        )
        assert "../results/A5297/annotation" in content

    def test_output_folder_explicit(self):
        content = generate_config_template(
            vcf_folder="/data/vcfs",
            output_folder="../results/myproject/annotation",
            dry_run=True,
        )
        assert "../results/myproject/annotation" in content

    def test_output_folder_default_fallback(self):
        content = generate_config_template(vcf_folder="results/final", dry_run=True)
        assert "results/annotation" in content
