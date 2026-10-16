"""Tests for config and samples schema validation."""

import copy
import re
from pathlib import Path

import pytest
import yaml

try:
    import jsonschema
except ImportError:
    jsonschema = None  # type: ignore[assignment]

REPO_ROOT = Path(__file__).resolve().parent.parent
SCHEMA_DIR = REPO_ROOT / "workflow" / "schemas"


def _load_schema(name):
    with open(SCHEMA_DIR / name) as f:
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

    def test_ref_optional_when_no_scatter(self, config_dict, config_schema):
        valid = copy.deepcopy(config_dict)
        valid["scatter"]["mode"] = "none"
        del valid["ref"]
        jsonschema.validate(valid, config_schema)

    def test_ref_required_when_interval_scatter(self, config_dict, config_schema):
        broken = copy.deepcopy(config_dict)
        broken["scatter"]["mode"] = "interval"
        del broken["ref"]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, config_schema)

    def test_valid_chromosome_mode(self, config_dict_chromosome, config_schema):
        jsonschema.validate(config_dict_chromosome, config_schema)

    def test_ref_required_when_chromosome_scatter(self, config_dict, config_schema):
        broken = copy.deepcopy(config_dict)
        broken["scatter"]["mode"] = "chromosome"
        del broken["ref"]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, config_schema)

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
        broken["extra_annotations"] = [{"vcf_file": "/db/ann.vcf.gz", "info_field": "AF"}]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, config_schema)

    def test_snpeff_canon_list_valid(self, config_dict, config_schema):
        valid = copy.deepcopy(config_dict)
        valid["snpeff"]["canon_list"] = "/resources/mane_canon.txt"
        jsonschema.validate(valid, config_schema)


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

    def test_sample_with_dot_rejected(self, samples_schema):
        broken = {"sample": "APA01.filtered", "vcf_basename": "APA01.filtered"}
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(broken, samples_schema)


class TestProfiles:
    def test_profiles_have_no_mem_mb_per_cpu_and_matching_resources(self):
        profile_files = sorted((REPO_ROOT / "profiles").glob("*/config.v8+.yaml"))
        assert profile_files
        for pf in profile_files:
            data = yaml.safe_load(pf.read_text())
            default_res = data.get("default-resources", {})
            assert "mem_mb_per_cpu" not in default_res, f"{pf} still sets mem_mb_per_cpu"
            set_threads = data.get("set-threads", {})
            set_resources = data.get("set-resources", {})
            for rule, threads in set_threads.items():
                assert rule in set_resources, f"{rule} missing from set-resources in {pf}"
                assert set_resources[rule]["cpus_per_task"] == threads
                assert set_resources[rule]["mem_mb"] >= 4000

    def test_profile_runtime_values_have_explicit_units(self):
        profile_files = sorted((REPO_ROOT / "profiles").glob("*/config*.yaml"))
        assert profile_files
        for pf in profile_files:
            data = yaml.safe_load(pf.read_text())
            default_res = data.get("default-resources", {})
            if "runtime" in default_res:
                rt = default_res["runtime"]
                assert isinstance(rt, str) and not rt.isdigit(), (
                    f"{pf} default-resources runtime is a bare integer: {rt!r}"
                )
            for rule, res in data.get("set-resources", {}).items():
                if "runtime" in res:
                    rt = res["runtime"]
                    assert isinstance(rt, str) and not rt.isdigit(), (
                        f"{pf} set-resources.{rule}.runtime is a bare integer: {rt!r}"
                    )

    def test_conda_envs_are_pinned_and_metadata_consistent(self):
        env_files = sorted((REPO_ROOT / "workflow" / "envs").glob("*.yaml"))
        assert env_files
        for ef in env_files:
            data = yaml.safe_load(ef.read_text())
            deps = data.get("dependencies", [])
            assert deps, f"{ef} has no dependencies"
            for dep in deps:
                assert re.match(r"^[A-Za-z0-9_.-]+=[A-Za-z0-9_.-]+$", str(dep)), (
                    f"Unpinned dependency {dep!r} in {ef}"
                )
            if ef.name == "snpeff.yaml":
                assert all(not str(d).startswith("numpy") for d in deps)

        pyproject_text = (REPO_ROOT / "pyproject.toml").read_text()
        m = re.search(r'^version\s*=\s*"([^"]+)"', pyproject_text, re.MULTILINE)
        assert m is not None, "pyproject.toml is missing [project] version"
        citation = yaml.safe_load((REPO_ROOT / "CITATION.cff").read_text())
        assert m.group(1) == str(citation["version"])

        ci_workflow = REPO_ROOT / ".github" / "workflows" / "ci.yml"
        assert ci_workflow.is_file()
        ci_data = yaml.safe_load(ci_workflow.read_text())
        steps = ci_data["jobs"]["ci"]["steps"]
        run_cmds = "\n".join(s.get("run", "") for s in steps)
        assert "make lint test-unit" in run_cmds
        for pkg in ("types-PyYAML", "types-jsonschema", "pandas==2.2.*", "numpy<2.5"):
            assert pkg in run_cmds, f"CI install missing {pkg!r}"
