#!/usr/bin/env python3
"""Run SnpSift tstv and write a MultiQC-compatible TSV.

Called via Snakemake's ``script:`` directive.  Receives the standard
``snakemake`` object with input/output/params/wildcards/log attributes.

SnpSift tstv writes its statistics table to **stdout**.  Progress messages
go to stderr.  We capture both streams explicitly via ``subprocess`` to
avoid the known limitations of Snakemake's ``shell(read=True)`` in
``run:`` blocks (see GitHub issue #23).
"""

import os
import subprocess
import sys

# Make workflow/rules/helpers.py importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "rules"))
from helpers import parse_snpsift_tstv


def main():
    vcf = snakemake.input.vcf
    output_tsv = snakemake.output.tstv
    log_file = snakemake.log[0]
    java_opts = snakemake.params.java_opts
    sample = snakemake.wildcards.sample

    with open(log_file, "w") as lf:
        lf.write("Starting snpsift_tstv\n")

    cmd = ["SnpSift", *java_opts.split(), "tstv", vcf]
    with open(log_file, "a") as lf:
        lf.write(f"Command: {' '.join(cmd)}\n")

    result = subprocess.run(cmd, capture_output=True, text=True)

    with open(log_file, "a") as lf:
        if result.stderr:
            lf.write(result.stderr)

    if result.returncode != 0:
        with open(log_file, "a") as lf:
            lf.write(f"SnpSift tstv failed with exit code {result.returncode}\n")
        sys.exit(result.returncode)

    # SnpSift tstv writes stats to stdout
    stdout = result.stdout
    if not stdout.strip():
        # Fallback: some versions may write to stderr
        stdout = result.stderr

    tstv = parse_snpsift_tstv(stdout)

    with open(output_tsv, "w") as fh:
        fh.write("# id: 'snpsift_tstv'\n")
        fh.write("# section_name: 'SnpSift Ts/Tv'\n")
        fh.write("# description: 'Transition / transversion ratios from SnpSift'\n")
        fh.write("# format: 'tsv'\n")
        fh.write("# plot_type: 'generalstats'\n")
        fh.write("# pconfig:\n")
        fh.write("#     - Transitions:\n")
        fh.write("#         title: 'Transitions'\n")
        fh.write("#         format: '{:,.0f}'\n")
        fh.write("#     - Transversions:\n")
        fh.write("#         title: 'Transversions'\n")
        fh.write("#         format: '{:,.0f}'\n")
        fh.write("#     - Ts/Tv:\n")
        fh.write("#         title: 'Ts/Tv'\n")
        fh.write("#         min: 0\n")
        fh.write("#         format: '{:.3f}'\n")
        cols = ["Transitions", "Transversions", "Ts/Tv"]
        fh.write("Sample\t" + "\t".join(cols) + "\n")
        vals = [tstv.get(c, "") for c in cols]
        fh.write(f"{sample}\t" + "\t".join(vals) + "\n")

    with open(log_file, "a") as lf:
        lf.write("Finished snpsift_tstv\n")


main()
