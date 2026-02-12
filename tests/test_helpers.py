"""Unit tests for workflow/rules/helpers.py."""

import os
import sys

# Add workflow/rules to the path so helpers.py can be imported directly
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "rules"))

from helpers import (
    annotation_step_vcf,
    format_extra_annotations,
    get_java_opts,
    get_samples,
    get_scatter_units,
    get_vcf_path,
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


class TestGetVcfPath:
    def test_resolves_path(self, samples_df):
        result = get_vcf_path("SampleA", "/data/vcfs", samples_df)
        assert result == os.path.join("/data/vcfs", "SampleA.vcf.gz")

    def test_different_basename(self):
        import pandas as pd

        df = pd.DataFrame(
            {"sample": ["S1"], "vcf_basename": ["different_name"]}
        ).set_index("sample", drop=False)
        result = get_vcf_path("S1", "/vcfs", df)
        assert result == os.path.join("/vcfs", "different_name.vcf.gz")


class TestAnnotationStepVcf:
    def test_step_zero_no_scatter(self):
        result = annotation_step_vcf("/out/ann", "SampleA", 0, "all")
        assert result == os.path.join("/out/ann", "SampleA.ann.dbnsfp.vcf.gz")

    def test_step_one_no_scatter(self):
        result = annotation_step_vcf("/out/ann", "SampleA", 1, "all")
        assert result == os.path.join("/out/ann", "SampleA.ann.step1.vcf.gz")

    def test_step_zero_with_scatter(self):
        result = annotation_step_vcf("/out/ann", "SampleA", 0, "0005-scattered")
        assert result == os.path.join(
            "/out/ann", "SampleA.0005-scattered.ann.dbnsfp.vcf.gz"
        )

    def test_step_two_with_scatter(self):
        result = annotation_step_vcf("/out/ann", "SampleA", 2, "0005-scattered")
        assert result == os.path.join(
            "/out/ann", "SampleA.0005-scattered.ann.step2.vcf.gz"
        )


class TestFormatExtraAnnotations:
    def test_empty(self):
        assert format_extra_annotations([]) == ""

    def test_single(self):
        annotations = [
            {
                "vcf_file": "/db/gnomad.vcf.gz",
                "info_field": "AF",
                "annotation_prefix": "GNOMAD_",
            }
        ]
        result = format_extra_annotations(annotations)
        assert result == "/db/gnomad.vcf.gz,AF,GNOMAD_"

    def test_multiple(self):
        annotations = [
            {"vcf_file": "/db/a.vcf.gz", "info_field": "AF", "annotation_prefix": "A_"},
            {"vcf_file": "/db/b.vcf.gz", "info_field": "AC", "annotation_prefix": "B_"},
        ]
        result = format_extra_annotations(annotations)
        assert result == "/db/a.vcf.gz,AF,A_ /db/b.vcf.gz,AC,B_"
