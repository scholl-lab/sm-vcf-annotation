"""snpEff annotation rules."""


rule snpeff_download_db:
    """Install the snpEff database once; parallel annotation jobs otherwise race on the same download."""
    output:
        # ponytail: marker only; if the snpeff conda env is rebuilt, delete it to re-download
        touch(os.path.join(LOG_DIR, f"snpeff_db.{SNPEFF_DB}.done")),
    log:
        os.path.join(LOG_DIR, f"snpeff_download_db.{SNPEFF_DB}.log"),
    params:
        java_opts=get_java_opts,
        db=SNPEFF_DB,
    conda:
        "../envs/snpeff.yaml"
    shell:
        "snpEff download {params.java_opts} -v {params.db} > {log} 2>&1"


rule snpeff_annotation:
    input:
        vcf_file=lambda wc: (
            get_input_vcf(wc)
            if wc.scatter_unit == "all"
            else os.path.join(ANNOTATION_DIR, f"{wc.sample}.{wc.scatter_unit}.vcf.gz")
        ),
        canon_list=SNPEFF_CANON_LIST if SNPEFF_CANON_LIST else [],
        snpeff_db=rules.snpeff_download_db.output,
    output:
        ann_vcf=temp(
            os.path.join(
                ANNOTATION_DIR,
                "{sample}.{scatter_unit}.ann.vcf.gz",
            )
        ),
        stats_csv=os.path.join(
            QC_DIR,
            "snpeff_stats.{sample}.{scatter_unit}.csv",
        ),
        genes_txt=os.path.join(
            QC_DIR,
            "snpeff_stats.{sample}.{scatter_unit}.genes.txt",
        ),
    log:
        os.path.join(LOG_DIR, "snpeff_annotation.{sample}.{scatter_unit}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/snpeff_annotation.{sample}.{scatter_unit}.tsv")
    params:
        java_opts=get_java_opts,
        db=SNPEFF_DB,
        extra_flags=SNPEFF_EXTRA_FLAGS,
        canon_arg=lambda wc, input: f"-onlyTr {input.canon_list}" if input.canon_list else "",
    conda:
        "../envs/snpeff.yaml"
    shell:
        r"""
        echo "Starting snpeff_annotation at: $(date)" >> {log}
        snpEff {params.java_opts} {params.db} {params.extra_flags} {params.canon_arg} \
            -csvStats {output.stats_csv} {input.vcf_file} \
        | bgzip -c > {output.ann_vcf} 2>> {log}
        echo "Finished snpeff_annotation at: $(date)" >> {log}
        """
