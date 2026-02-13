# Configuration

## `config.yaml`

Unified pipeline configuration. All paths are relative to the repository root.

### Required Fields

| Field | Type | Description |
|-------|------|-------------|
| `paths.samples` | string | Path to samples TSV |
| `paths.vcf_folder` | string | Directory with input VCF files |
| `paths.output_folder` | string | Root output directory |
| `snpeff.database` | string | snpEff database name (e.g., `GRCh37.p13`) |
| `snpsift.dbnsfp_db` | string | Path to dbNSFP database file (v4.x or v5.x) |
| `snpsift.dbnsfp_fields` | string | Comma-separated dbNSFP fields (auto-selected for v4/v5 by config generator) |
| `scatter.mode` | string | `"none"`, `"chromosome"`, or `"interval"` |

### Required for Scatter Modes (`scatter.mode: "chromosome"` or `"interval"`)

| Field | Description |
|-------|-------------|
| `ref.genome` | Path to reference genome FASTA |
| `ref.dict` | Path to reference sequence dictionary |
| `scatter.count` | Number of scatter intervals (default: 100, interval mode only) |
| `scatter.canonical_contigs` | List of canonical contigs to include (used by both chromosome and interval modes) |

### Optional Fields

| Field | Default | Description |
|-------|---------|-------------|
| `paths.log_subdir` | `"logs"` | Log subdirectory name |
| `paths.annotation_subdir` | `"annotation"` | Annotation subdirectory name |
| `snpeff.extra_flags` | `""` | Additional snpEff flags |
| `extra_annotations` | `[]` | List of extra VCF annotation steps |

### dbNSFP Field Lists

The config generator auto-detects the dbNSFP major version from the database filename and selects the appropriate field list:

| dbNSFP version | Key differences |
|----------------|-----------------|
| **v4.x** | LRT, FATHMM, gnomAD exomes/genomes (separate), MutationTaster_converted_rankscore |
| **v5.x** | LRT retired, FATHMM replaced by fathmm-XF, gnomAD merged into gnomAD4.1 joint, MutationTaster_rankscore |

If you change the dbNSFP database version after generating config, regenerate with `scripts/generate_config.py` or manually update `snpsift.dbnsfp_fields`.

### Extra Annotations

Each entry in `extra_annotations` applies one `SnpSift annotate` step:

```yaml
extra_annotations:
  - vcf_file: "/path/to/annotation.vcf.gz"
    info_field: "AF"
    annotation_prefix: "MYDB_"
```

## `samples.tsv`

Tab-separated sample metadata. Required columns:

| Column | Description |
|--------|-------------|
| `sample` | Unique sample identifier |
| `vcf_basename` | VCF file basename (without `.vcf.gz` extension) |
