"""Shared pytest fixtures for sm-vcf-annotation tests."""

import pandas as pd
import pytest


@pytest.fixture
def samples_df():
    """Sample DataFrame matching config/samples.tsv schema."""
    return pd.DataFrame(
        {
            "sample": ["SampleA", "SampleB"],
            "vcf_basename": ["SampleA", "SampleB"],
        }
    ).set_index("sample", drop=False)


@pytest.fixture
def config_dict():
    """Minimal valid configuration dictionary."""
    return {
        "ref": {
            "genome": "/ref/GRCh37.p13.genome.fa",
            "dict": "/ref/GRCh37.p13.genome.dict",
        },
        "paths": {
            "samples": "config/samples.tsv",
            "vcf_folder": "results/final",
            "output_folder": "results/annotation",
            "log_subdir": "logs",
            "annotation_subdir": "annotation",
        },
        "snpeff": {
            "database": "GRCh37.p13",
            "extra_flags": "-lof -noInteraction",
        },
        "snpsift": {
            "dbnsfp_db": "dbnsfp/dbNSFP4.9a_grch37.gz",
            "dbnsfp_fields": "aaref,aaalt,SIFT_pred",
        },
        "scatter": {
            "mode": "none",
            "count": 100,
            "canonical_contigs": [],
        },
        "extra_annotations": [],
    }
