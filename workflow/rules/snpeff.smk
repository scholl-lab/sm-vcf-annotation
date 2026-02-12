"""snpEff annotation rule."""


rule snpeff_annotation:
    input:
        vcf_file=lambda wc: (
            get_input_vcf(wc)
            if wc.scatter_unit == "all"
            else os.path.join(ANNOTATION_DIR, f"{wc.sample}.{wc.scatter_unit}.vcf.gz")
        ),
    output:
        ann_vcf=temp(
            os.path.join(
                ANNOTATION_DIR,
                "{sample}.{scatter_unit}.ann.vcf.gz",
            )
        ),
    log:
        os.path.join(LOG_DIR, "snpeff_annotation.{sample}.{scatter_unit}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/snpeff_annotation.{sample}.{scatter_unit}.tsv")
    params:
        java_opts=get_java_opts,
        db=SNPEFF_DB,
        extra_flags=SNPEFF_EXTRA_FLAGS,
        stats_html=os.path.join(LOG_DIR, "snpeff_stats.{sample}.{scatter_unit}.html"),
    conda:
        "../envs/snpeff.yaml"
    shell:
        r"""
        echo "Starting snpeff_annotation at: $(date)" >> {log}
        snpEff {params.java_opts} {params.db} {params.extra_flags} \
            -stats {params.stats_html} {input.vcf_file} \
        | bgzip -c > {output.ann_vcf} 2>> {log}
        echo "Finished snpeff_annotation at: $(date)" >> {log}
        """
