"""Unit tests for workflow/rules/helpers.py."""

import os
import sys

# Add workflow/rules to the path so helpers.py can be imported directly
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "rules"))

from helpers import (
    annotation_step_vcf,
    get_java_opts,
    get_samples,
    get_scatter_units,
    get_vcf_path,
    parse_annotation_completeness,
    parse_snpsift_tstv,
)


class TestGetJavaOpts:
    def test_standard_allocation(self):
        result = get_java_opts(32768, "/tmp")
        assert "-Xms6553m" in result
        assert "-Xmx26214m" in result
        assert "-Djava.io.tmpdir=/tmp" in result

    def test_small_allocation(self):
        result = get_java_opts(1000, "/scratch")
        assert "-Xms200m" in result
        assert "-Xmx800m" in result


class TestGetSamples:
    def test_returns_sorted(self, samples_df):
        result = get_samples(samples_df)
        assert result == ["SampleA", "SampleB"]

    def test_deduplicates(self, samples_df):
        result = get_samples(samples_df)
        assert len(result) == 2


class TestGetScatterUnits:
    def test_interval_mode(self):
        result = get_scatter_units("interval", scatter_count=5)
        assert result == [
            "0000-scattered",
            "0001-scattered",
            "0002-scattered",
            "0003-scattered",
            "0004-scattered",
        ]

    def test_none_mode(self):
        result = get_scatter_units("none")
        assert result == ["all"]

    def test_interval_default_count(self):
        result = get_scatter_units("interval")
        assert len(result) == 100
        assert result[0] == "0000-scattered"
        assert result[99] == "0099-scattered"

    def test_chromosome_mode(self):
        chroms = ["chr1", "chr2", "chr3", "chrX", "chrY"]
        result = get_scatter_units("chromosome", chromosomes=chroms)
        assert result == chroms

    def test_chromosome_mode_no_chromosomes_raises(self):
        import pytest

        with pytest.raises(ValueError, match="non-empty"):
            get_scatter_units("chromosome")

    def test_chromosome_mode_empty_list_raises(self):
        import pytest

        with pytest.raises(ValueError, match="non-empty"):
            get_scatter_units("chromosome", chromosomes=[])

    def test_chromosome_mode_grch37_contigs(self):
        contigs = ["1", "2", "3", "X", "Y"]
        result = get_scatter_units("chromosome", chromosomes=contigs)
        assert result == contigs

    def test_interval_mode_unaffected_by_chromosomes(self):
        result = get_scatter_units("interval", chromosomes=["chr1"], scatter_count=3)
        assert result == ["0000-scattered", "0001-scattered", "0002-scattered"]

    def test_none_mode_unaffected_by_chromosomes(self):
        result = get_scatter_units("none", chromosomes=["chr1"])
        assert result == ["all"]


class TestGetVcfPath:
    def test_resolves_path(self, samples_df):
        result = get_vcf_path("SampleA", "/data/vcfs", samples_df)
        assert result == os.path.join("/data/vcfs", "SampleA.vcf.gz")

    def test_different_basename(self):
        import pandas as pd

        df = pd.DataFrame({"sample": ["S1"], "vcf_basename": ["different_name"]}).set_index(
            "sample", drop=False
        )
        result = get_vcf_path("S1", "/vcfs", df)
        assert result == os.path.join("/vcfs", "different_name.vcf.gz")


class TestAnnotationStepVcf:
    def test_step_zero_no_scatter(self):
        result = annotation_step_vcf("/out/ann", "SampleA", 0, "all")
        assert result == os.path.join("/out/ann", "SampleA.all.ann.dbnsfp.vcf.gz")

    def test_step_one_no_scatter(self):
        result = annotation_step_vcf("/out/ann", "SampleA", 1, "all")
        assert result == os.path.join("/out/ann", "SampleA.all.ann.step1.vcf.gz")

    def test_step_zero_with_scatter(self):
        result = annotation_step_vcf("/out/ann", "SampleA", 0, "0005-scattered")
        assert result == os.path.join("/out/ann", "SampleA.0005-scattered.ann.dbnsfp.vcf.gz")

    def test_step_two_with_scatter(self):
        result = annotation_step_vcf("/out/ann", "SampleA", 2, "0005-scattered")
        assert result == os.path.join("/out/ann", "SampleA.0005-scattered.ann.step2.vcf.gz")

    def test_step_zero_with_chromosome(self):
        result = annotation_step_vcf("/out/ann", "SampleA", 0, "chr1")
        assert result == os.path.join("/out/ann", "SampleA.chr1.ann.dbnsfp.vcf.gz")

    def test_step_one_with_chromosome(self):
        result = annotation_step_vcf("/out/ann", "SampleA", 1, "chr1")
        assert result == os.path.join("/out/ann", "SampleA.chr1.ann.step1.vcf.gz")


class TestParseAnnotationCompleteness:
    def test_all_annotated(self):
        query = "chr1\t100\tA|missense\tD\t0.95\nchr1\t200\tA|synonymous\tT\t0.5\n"
        fields = ["ANN", "SIFT_pred", "REVEL_score"]
        result = parse_annotation_completeness(query, fields)
        assert result["ANN"]["total"] == 2
        assert result["ANN"]["annotated"] == 2
        assert result["ANN"]["rate"] == 1.0
        assert result["SIFT_pred"]["annotated"] == 2
        assert result["REVEL_score"]["annotated"] == 2

    def test_partial_annotation(self):
        query = "chr1\t100\tA|missense\t.\t0.95\nchr1\t200\t.\tT\t.\n"
        fields = ["ANN", "SIFT_pred", "REVEL_score"]
        result = parse_annotation_completeness(query, fields)
        assert result["ANN"]["total"] == 2
        assert result["ANN"]["annotated"] == 1
        assert result["ANN"]["rate"] == 0.5
        assert result["SIFT_pred"]["annotated"] == 1
        assert result["REVEL_score"]["annotated"] == 1

    def test_empty_input(self):
        result = parse_annotation_completeness("", ["ANN", "SIFT_pred"])
        assert result["ANN"]["total"] == 0
        assert result["ANN"]["annotated"] == 0
        assert result["ANN"]["rate"] == 0.0

    def test_all_missing(self):
        query = "chr1\t100\t.\t.\nchr1\t200\t.\t.\n"
        fields = ["ANN", "SIFT_pred"]
        result = parse_annotation_completeness(query, fields)
        assert result["ANN"]["annotated"] == 0
        assert result["SIFT_pred"]["annotated"] == 0
        assert result["ANN"]["rate"] == 0.0

    def test_single_field(self):
        query = "chr1\t100\tA|missense\nchr1\t200\t.\nchr1\t300\tB|stop\n"
        fields = ["ANN"]
        result = parse_annotation_completeness(query, fields)
        assert result["ANN"]["total"] == 3
        assert result["ANN"]["annotated"] == 2
        assert result["ANN"]["rate"] == round(2 / 3, 4)


class TestParseSnpsiftTstv:
    def test_multi_sample(self):
        output = (
            "Sample        : 1       2       Total\n"
            "Transitions   : 150488  160000  310488\n"
            "Transversions : 70878   80000   150878\n"
            "Ts/Tv         : 2.123   2.000   2.058\n"
        )
        result = parse_snpsift_tstv(output)
        assert result["Transitions"] == "310488"
        assert result["Transversions"] == "150878"
        assert result["Ts/Tv"] == "2.058"

    def test_single_sample(self):
        output = (
            "Sample        : 1       Total\n"
            "Transitions   : 150488  150488\n"
            "Transversions : 70878   70878\n"
            "Ts/Tv         : 2.123   2.123\n"
        )
        result = parse_snpsift_tstv(output)
        assert result["Transitions"] == "150488"
        assert result["Transversions"] == "70878"
        assert result["Ts/Tv"] == "2.123"

    def test_empty_input(self):
        result = parse_snpsift_tstv("")
        assert result == {}

    def test_skips_sample_row(self):
        output = (
            "Sample        : 1       Total\n"
            "Transitions   : 100     100\n"
        )
        result = parse_snpsift_tstv(output)
        assert "Sample" not in result
        assert result["Transitions"] == "100"
