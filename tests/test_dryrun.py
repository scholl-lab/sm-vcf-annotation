"""Snakemake dry-run integration test."""

import os
import subprocess

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.mark.dryrun
def test_dryrun_no_scatter():
    """Verify the workflow DAG resolves with scatter.mode=none."""
    result = subprocess.run(
        [
            "snakemake",
            "-s",
            "workflow/Snakefile",
            "--configfile",
            "tests/data/config/config.yaml",
            "-n",
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
