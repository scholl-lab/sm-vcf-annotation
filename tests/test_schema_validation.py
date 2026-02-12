"""Tests for config and samples schema validation."""

import copy

import pytest
import yaml

try:
    import jsonschema
except ImportError:
    jsonschema = None

SCHEMA_DIR = "workflow/schemas"


def _load_schema(name):
    with open(f"{SCHEMA_DIR}/{name}") as f:
        return yaml.safe_load(f)


@pytest.fixture
def config_schema():
    return _load_schema("config.schema.yaml")


@pytest.fixture
def samples_schema():
    return _load_schema("samples.schema.yaml")


@pytest.mark.skipif(jsonschema is None, reason="jsonschema not installed")
class TestConfigSchema:
    def test_valid_config(self, config_dict, config_schema):
        jsonschema.validate(config_dict, config_schema)

    def test_missing_paths(self, config_dict, config_schema):
        broken = copy.deepcopy(config_dict)
        del broken["paths"]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, config_schema)

    def test_missing_snpeff(self, config_dict, config_schema):
        broken = copy.deepcopy(config_dict)
        del broken["snpeff"]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, config_schema)

    def test_missing_snpsift(self, config_dict, config_schema):
        broken = copy.deepcopy(config_dict)
        del broken["snpsift"]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, config_schema)

    def test_missing_scatter(self, config_dict, config_schema):
        broken = copy.deepcopy(config_dict)
        del broken["scatter"]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, config_schema)

    def test_invalid_scatter_mode(self, config_dict, config_schema):
        broken = copy.deepcopy(config_dict)
        broken["scatter"]["mode"] = "invalid"
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, config_schema)

    def test_ref_optional(self, config_dict, config_schema):
        valid = copy.deepcopy(config_dict)
        del valid["ref"]
        jsonschema.validate(valid, config_schema)

    def test_extra_annotations_valid(self, config_dict, config_schema):
        valid = copy.deepcopy(config_dict)
        valid["extra_annotations"] = [
            {
                "vcf_file": "/db/ann.vcf.gz",
                "info_field": "AF",
                "annotation_prefix": "DB_",
            }
        ]
        jsonschema.validate(valid, config_schema)

    def test_extra_annotations_missing_field(self, config_dict, config_schema):
        broken = copy.deepcopy(config_dict)
        broken["extra_annotations"] = [
            {"vcf_file": "/db/ann.vcf.gz", "info_field": "AF"}
        ]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, config_schema)


@pytest.mark.skipif(jsonschema is None, reason="jsonschema not installed")
class TestSamplesSchema:
    def test_valid_sample(self, samples_schema):
        sample = {"sample": "SampleA", "vcf_basename": "SampleA"}
        jsonschema.validate(sample, samples_schema)

    def test_missing_sample(self, samples_schema):
        broken = {"vcf_basename": "SampleA"}
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, samples_schema)

    def test_missing_vcf_basename(self, samples_schema):
        broken = {"sample": "SampleA"}
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, samples_schema)
