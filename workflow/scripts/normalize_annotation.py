"""Normalize INFO fields added by SnpSift annotate (stdin VCF text -> stdout).

Usage: normalize_annotation.py <comma-separated new INFO IDs> [ID:max|min|first,...]
See helpers.normalize_annotation_vcf.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "rules"))
from helpers import AGGREGATORS, normalize_annotation_vcf

new_ids = [i for i in sys.argv[1].split(",") if i]
collapse = dict(
    c.split(":", 1) for c in (sys.argv[2].split(",") if len(sys.argv) > 2 and sys.argv[2] else [])
)
bad = {k: v for k, v in collapse.items() if v not in AGGREGATORS or k not in new_ids}
if bad:
    sys.exit(f"normalize_annotation.py: invalid collapse spec {bad}")
sys.stdout.writelines(normalize_annotation_vcf(sys.stdin, new_ids, collapse))
