# Deprecated Files

These files are from the legacy pipeline layout and have been replaced by the
modernized `workflow/` structure. They are kept here for reference.

## Migration Table

| Old Path | New Path |
|----------|----------|
| `scripts/snakemake/annotate_snpeff_snpsift.smk` | `workflow/Snakefile` + `workflow/rules/*.smk` |
| `scripts/snakemake/scattered_annotate.smk` | `workflow/Snakefile` + `workflow/rules/scatter.smk` |
| `scripts/snakemake/config_example.yaml` | `config/config.yaml` |
| `scripts/run_annotate_snpeff_snpsift.sh` | `scripts/run_snakemake.sh` |
| `scripts/run_scattered_annotate.sh` | `scripts/run_snakemake.sh` |
| `config.yaml` (flat keys) | `config/config.yaml` (nested keys) |

## Key Changes

- Two separate workflow files merged into one unified workflow with `scatter.mode` config key
- Named conda environments replaced with YAML env files in `workflow/envs/`
- Hardcoded resources moved to `profiles/` directories
- Schema validation added for config and samples
- `samples.tsv` replaces directory globbing for sample discovery
