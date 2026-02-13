"""QC rules: bcftools stats, SnpSift tstv, annotation completeness, and MultiQC."""

from rules.helpers import parse_annotation_completeness, parse_snpsift_tstv


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
    run:
        import subprocess

        shell("echo 'Starting snpsift_tstv at: $(date)' > {log}")
        result = subprocess.run(
            f"SnpSift {params.java_opts} tstv {input.vcf}",
            shell=True,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            with open(str(log), "a") as lf:
                lf.write(f"SnpSift tstv stderr: {result.stderr}\n")
        tstv = parse_snpsift_tstv(result.stdout)
        with open(str(output.tstv), "w") as fh:
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
            fh.write(f"{wildcards.sample}\t" + "\t".join(vals) + "\n")
        shell("echo 'Finished snpsift_tstv at: $(date)' >> {log}")


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
    run:
        import subprocess

        shell("echo 'Starting annotation_completeness at: $(date)' > {log}")
        field_fmt = "\t".join([f"%{f}" for f in params.fields])
        query_cmd = f"bcftools query -f '%CHROM\\t%POS\\t{field_fmt}\\n' {input.vcf}"
        result = subprocess.run(
            query_cmd,
            shell=True,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            with open(str(log), "a") as lf:
                lf.write(f"bcftools query stderr: {result.stderr}\n")
        query_output = result.stdout
        stats = parse_annotation_completeness(query_output, list(params.fields))

        with open(str(output.tsv), "w") as fh:
            fh.write("# id: 'annotation_completeness'\n")
            fh.write("# section_name: 'Annotation Completeness'\n")
            fh.write("# description: 'Fraction of variants with non-missing annotation values'\n")
            fh.write("# format: 'tsv'\n")
            fh.write("# plot_type: 'generalstats'\n")
            fh.write("# pconfig:\n")
            for field in params.fields:
                fh.write(f"#     - {field}_rate:\n")
                fh.write(f"#         title: '{field} rate'\n")
                fh.write("#         min: 0\n")
                fh.write("#         max: 1\n")
                fh.write("#         format: '{:.2%}'\n")
            fh.write("Sample\t" + "\t".join([f"{f}_rate" for f in params.fields]) + "\n")
            rates = [str(stats[f]["rate"]) for f in params.fields]
            fh.write(f"{wildcards.sample}\t" + "\t".join(rates) + "\n")
        shell("echo 'Finished annotation_completeness at: $(date)' >> {log}")


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
