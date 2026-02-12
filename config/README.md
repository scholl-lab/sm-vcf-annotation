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
| `snpsift.dbnsfp_db` | string | Path to dbNSFP database file |
| `snpsift.dbnsfp_fields` | string | Comma-separated dbNSFP fields |
| `scatter.mode` | string | `"interval"` or `"none"` |

### Required for Scatter Mode (`scatter.mode: "interval"`)

| Field | Description |
|-------|-------------|
| `ref.genome` | Path to reference genome FASTA |
| `ref.dict` | Path to reference sequence dictionary |
| `scatter.count` | Number of scatter intervals (default: 100) |
| `scatter.canonical_contigs` | List of canonical contigs to include |

### Optional Fields

| Field | Default | Description |
|-------|---------|-------------|
| `paths.log_subdir` | `"logs"` | Log subdirectory name |
| `paths.annotation_subdir` | `"annotation"` | Annotation subdirectory name |
| `snpeff.extra_flags` | `""` | Additional snpEff flags |
| `extra_annotations` | `[]` | List of extra VCF annotation steps |

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
