"""QC rules: bcftools stats, SnpSift tstv, annotation completeness, and MultiQC."""


rule bcftools_stats:
    input:
        vcf=os.path.join(ANNOTATION_DIR, "{sample}.annotated.vcf.gz"),
    output:
        stats=os.path.join(QC_DIR, "{sample}.bcftools_stats.txt"),
    log:
        os.path.join(LOG_DIR, "bcftools_stats.{sample}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/bcftools_stats.{sample}.tsv")
    conda:
        "../envs/snpeff.yaml"
    shell:
        r"""
        echo "Starting bcftools_stats at: $(date)" > {log}
        bcftools stats --threads {threads} {input.vcf} > {output.stats} 2>> {log}
        echo "Finished bcftools_stats at: $(date)" >> {log}
        """


rule snpsift_tstv:
    input:
        vcf=os.path.join(ANNOTATION_DIR, "{sample}.annotated.vcf.gz"),
    output:
        tstv=os.path.join(QC_DIR, "{sample}.snpsift_tstv_mqc.tsv"),
    log:
        os.path.join(LOG_DIR, "snpsift_tstv.{sample}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/snpsift_tstv.{sample}.tsv")
    params:
        java_opts=get_java_opts,
    conda:
        "../envs/snpeff.yaml"
    script:
        "../scripts/snpsift_tstv.py"


rule annotation_completeness:
    input:
        vcf=os.path.join(ANNOTATION_DIR, "{sample}.annotated.vcf.gz"),
    output:
        tsv=os.path.join(QC_DIR, "{sample}.annotation_completeness_mqc.tsv"),
    log:
        os.path.join(LOG_DIR, "annotation_completeness.{sample}.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/annotation_completeness.{sample}.tsv")
    params:
        fields=["ANN", "dbNSFP_SIFT_pred", "dbNSFP_REVEL_score"],
    conda:
        "../envs/snpeff.yaml"
    script:
        "../scripts/annotation_completeness.py"


rule multiqc:
    input:
        bcftools=expand(
            os.path.join(QC_DIR, "{sample}.bcftools_stats.txt"),
            sample=get_sample_list(),
        ),
        tstv=expand(
            os.path.join(QC_DIR, "{sample}.snpsift_tstv_mqc.tsv"),
            sample=get_sample_list(),
        ),
        completeness=expand(
            os.path.join(QC_DIR, "{sample}.annotation_completeness_mqc.tsv"),
            sample=get_sample_list(),
        ),
        snpeff_csv=expand(
            os.path.join(QC_DIR, "snpeff_stats.{sample}.{scatter_unit}.csv"),
            sample=get_sample_list(),
            scatter_unit=SCATTER_UNITS,
        ),
    output:
        html=os.path.join(QC_DIR, "multiqc_report.html"),
        data=directory(os.path.join(QC_DIR, "multiqc_data")),
    log:
        os.path.join(LOG_DIR, "multiqc.log"),
    benchmark:
        os.path.join(LOG_DIR, "benchmarks/multiqc.tsv")
    params:
        qc_dir=QC_DIR,
        mqc_config=workflow.source_path("../report/multiqc_config.yaml"),
    conda:
        "../envs/multiqc.yaml"
    shell:
        r"""
        echo "Starting multiqc at: $(date)" > {log}
        multiqc {params.qc_dir} \
            --config {params.mqc_config} \
            --outdir {params.qc_dir} \
            --force 2>> {log}
        echo "Finished multiqc at: $(date)" >> {log}
        """
