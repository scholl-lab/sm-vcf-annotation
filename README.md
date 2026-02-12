# sm-vcf-annotation

Snakemake 8+ pipeline for annotating VCF files using snpEff, SnpSift, and dbNSFP.

## Features

- **Unified workflow** — single Snakefile handles both single-pass and scattered (interval-based) annotation
- **Profile-driven resources** — no hardcoded threads/memory; all tuning via `profiles/`
- **Schema validation** — config and samples validated at startup
- **YAML conda environments** — pinned, reproducible tool versions
- **Cluster auto-detection** — launcher script detects BIH, Charité, or local execution

## Quick Start

### 1. Generate config files

```bash
# Interactive wizard
python scripts/generate_config.py

# Or with flags
python scripts/generate_config.py --vcf-folder /path/to/vcfs
```

### 2. Edit configuration

Edit `config/config.yaml` to set database paths:

```yaml
snpeff:
  database: "GRCh37.p13"
snpsift:
  dbnsfp_db: "/path/to/dbNSFP4.9a_grch37.gz"
```

For scatter mode (large multisample VCFs), also set:

```yaml
ref:
  genome: "/path/to/reference.fa"
  dict: "/path/to/reference.dict"
scatter:
  mode: "interval"
  count: 100
  canonical_contigs: ["chr1", "chr2", ...]
```

### 3. Run the pipeline

```bash
# On HPC (auto-detects cluster)
sbatch scripts/run_snakemake.sh

# With custom config
sbatch scripts/run_snakemake.sh config/config.yaml

# Local execution
snakemake -s workflow/Snakefile --workflow-profile profiles/default --profile profiles/local
```

## Configuration Reference

See [`config/README.md`](config/README.md) for full field reference.

### Config Key Migration (from legacy layout)

| Old Key | New Key |
|---------|---------|
| `snpeff_annotation_db` | `snpeff.database` |
| `snpeff_additional_flags` | `snpeff.extra_flags` |
| `snpsift_db_location` | `snpsift.dbnsfp_db` |
| `snpsift_dbnsfp_fields` | `snpsift.dbnsfp_fields` |
| `vcf_input_folder` | `paths.vcf_folder` |
| `output_folder` | `paths.output_folder` |
| `annotation_subfolder` | `paths.annotation_subdir` |
| `log_subfolder` | `paths.log_subdir` |
| `conda_environment_annotation` | *(removed — uses workflow/envs/*.yaml)* |
| `extra_vcf_annotations` | `extra_annotations` |
| `reference_genome` | `ref.genome` |
| `reference_dict` | `ref.dict` |
| `scatter_count` | `scatter.count` |
| `canonical_contigs` | `scatter.canonical_contigs` |

## Profiles

| Profile | Description |
|---------|-------------|
| `profiles/default/` | Workflow profile: resource allocation, conda, shared-fs |
| `profiles/bih/` | BIH HPC: SLURM executor, medium partition, scholl-lab account |
| `profiles/charite/` | Charité HPC: SLURM executor, compute partition |
| `profiles/local/` | Local execution: 4 jobs, conda only |

## Repository Structure

```
workflow/
  Snakefile              # Entry point (Snakemake 8+)
  rules/
    common.smk           # Config shortcuts, input functions
    snpeff.smk           # snpEff annotation rule
    snpsift.smk          # SnpSift varType, dbNSFP, extra annotations
    scatter.smk          # Conditional scatter/gather (GATK)
    helpers.py           # Pure Python helpers (unit-testable)
  envs/
    snpeff.yaml          # snpEff + SnpSift + bcftools
    gatk.yaml            # GATK4 + samtools
  schemas/
    config.schema.yaml   # Config validation
    samples.schema.yaml  # Samples validation
config/
  config.yaml            # Pipeline configuration
  samples.tsv            # Sample sheet
  README.md              # Config field reference
profiles/
  default/               # Workflow profile (resource allocation)
  bih/                   # BIH HPC cluster profile
  charite/               # Charité HPC cluster profile
  local/                 # Local execution profile
scripts/
  run_snakemake.sh       # Unified launcher with cluster detection
  generate_config.py     # Config file generator
tests/
  test_helpers.py        # Unit tests for helpers.py
  test_schema_validation.py  # Schema validation tests
  test_generate_config.py    # Config generator tests
  test_dryrun.py         # Snakemake dry-run integration test
deprecated/              # Legacy files (for reference)
```

## Development

```bash
# Install dev tools
make install-dev

# Run all checks
make lint

# Auto-format
make format

# Run unit tests
make test-unit

# Run all tests (including dry-run)
make test
```

## Annotation Pipeline

The workflow performs these steps per sample (per scatter interval if `scatter.mode: "interval"`):

1. **snpEff** — variant annotation (gene, transcript, functional impact)
2. **SnpSift varType** — variant type classification (SNP, INS, DEL, etc.)
3. **SnpSift dbNSFP** — functional prediction scores (SIFT, PolyPhen, CADD, etc.)
4. **Extra annotations** (optional) — step-wise SnpSift annotate from custom VCFs
5. **Finalize** — copy/rename to `.annotated.vcf.gz`
6. **Concatenate** (scatter mode only) — merge intervals back into one VCF per sample

Intermediate files use Snakemake `temp()` for automatic cleanup.
