#!/usr/bin/env python3
"""Compute annotation completeness and write a MultiQC-compatible TSV.

Called via Snakemake's ``script:`` directive.  Uses ``subprocess.run``
instead of Snakemake's ``shell(read=True)`` for reliable stdout capture
(same rationale as snpsift_tstv.py — see GitHub issue #23).
"""

import os
import subprocess
import sys

# Make workflow/rules/helpers.py importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "rules"))
from helpers import parse_annotation_completeness


def main():
    vcf = snakemake.input.vcf
    output_tsv = snakemake.output.tsv
    log_file = snakemake.log[0]
    fields = list(snakemake.params.fields)
    sample = snakemake.wildcards.sample

    with open(log_file, "w") as lf:
        lf.write("Starting annotation_completeness\n")

    field_fmt = "\\t".join([f"%{f}" for f in fields])
    cmd = f"bcftools query -f '%CHROM\\t%POS\\t{field_fmt}\\n' {vcf}"

    with open(log_file, "a") as lf:
        lf.write(f"Command: {cmd}\n")

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    with open(log_file, "a") as lf:
        if result.stderr:
            lf.write(result.stderr)

    if result.returncode != 0:
        with open(log_file, "a") as lf:
            lf.write(f"bcftools query failed with exit code {result.returncode}\n")
        sys.exit(result.returncode)

    stats = parse_annotation_completeness(result.stdout.splitlines(), fields)

    with open(output_tsv, "w") as fh:
        fh.write("# id: 'annotation_completeness'\n")
        fh.write("# section_name: 'Annotation Completeness'\n")
        fh.write("# description: 'Fraction of variants with non-missing annotation values'\n")
        fh.write("# format: 'tsv'\n")
        fh.write("# plot_type: 'generalstats'\n")
        fh.write("# pconfig:\n")
        for field in fields:
            fh.write(f"#     - {field}_rate:\n")
            fh.write(f"#         title: '{field} rate'\n")
            fh.write("#         min: 0\n")
            fh.write("#         max: 1\n")
            fh.write("#         format: '{:.2%}'\n")
        fh.write("Sample\t" + "\t".join([f"{f}_rate" for f in fields]) + "\n")
        rates = [str(stats[f]["rate"]) for f in fields]
        fh.write(f"{sample}\t" + "\t".join(rates) + "\n")

    with open(log_file, "a") as lf:
        lf.write("Finished annotation_completeness\n")


main()
