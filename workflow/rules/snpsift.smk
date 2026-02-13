"""SnpSift annotation rules: varType, dbNSFP, extra annotations, and finalize."""


rule snpsift_variant_type:
    input:
        ann_vcf=os.path.join(ANNOTATION_DIR, "{sample}.{scatter_unit}.ann.vcf.gz"),
    output:
        ann_vartype_vcf=temp(
            os.path.join(
                ANNOTATION_DIR,
                "{sample}.{scatter_unit}.ann.vartype.vcf.gz",
            )
        ),
    log:
        os.path.join(LOG_DIR, "snpsift_variant_type.{sample}.{scatter_unit}.log"),
    benchmark:
        os.path.join(
            LOG_DIR,
            "benchmarks/snpsift_variant_type.{sample}.{scatter_unit}.tsv",
        )
    params:
        java_opts=get_java_opts,
    conda:
        "../envs/snpeff.yaml"
    shell:
        r"""
        echo "Starting snpsift_variant_type at: $(date)" >> {log}
        SnpSift {params.java_opts} varType {input.ann_vcf} \
        | bgzip -c > {output.ann_vartype_vcf} 2>> {log}
        echo "Finished snpsift_variant_type at: $(date)" >> {log}
        """


rule snpsift_annotation_dbnsfp:
    input:
        ann_vartype_vcf=os.path.join(
            ANNOTATION_DIR,
            "{sample}.{scatter_unit}.ann.vartype.vcf.gz",
        ),
    output:
        ann_dbnsfp_vcf=temp(
            os.path.join(
                ANNOTATION_DIR,
                "{sample}.{scatter_unit}.ann.dbnsfp.vcf.gz",
            )
        ),
    log:
        os.path.join(
            LOG_DIR,
            "snpsift_annotation_dbnsfp.{sample}.{scatter_unit}.log",
        ),
    benchmark:
        os.path.join(
            LOG_DIR,
            "benchmarks/snpsift_annotation_dbnsfp.{sample}.{scatter_unit}.tsv",
        )
    params:
        java_opts=get_java_opts,
        dbnsfp_db=SNPSIFT_DBNSFP_DB,
        dbnsfp_fields=SNPSIFT_DBNSFP_FIELDS,
    conda:
        "../envs/snpeff.yaml"
    shell:
        r"""
        echo "Starting snpsift_annotation_dbnsfp at: $(date)" >> {log}
        SnpSift {params.java_opts} dbnsfp \
            -f {params.dbnsfp_fields} \
            -db {params.dbnsfp_db} \
            {input.ann_vartype_vcf} \
        | bgzip -c > {output.ann_dbnsfp_vcf} 2>> {log}
        echo "Finished snpsift_annotation_dbnsfp at: $(date)" >> {log}
        """


# ---------------------------------------------------------------------------
# Step-wise extra annotations
# ---------------------------------------------------------------------------
if N_ANNOTATIONS > 0:

    rule annotation_step:
        input:
            get_annotation_step_input,
        output:
            temp(
                os.path.join(
                    ANNOTATION_DIR,
                    "{sample}.{scatter_unit}.ann.step{step}.vcf.gz",
                )
            ),
        log:
            os.path.join(
                LOG_DIR,
                "extra_annotations.step{step}.{sample}.{scatter_unit}.log",
            ),
        benchmark:
            os.path.join(
                LOG_DIR,
                "benchmarks/extra_annotations.step{step}.{sample}.{scatter_unit}.tsv",
            )
        params:
            java_opts=get_java_opts,
            vcf_file=lambda wc: EXTRA_ANNOTATIONS[int(wc.step) - 1]["vcf_file"],
            info_field=lambda wc: EXTRA_ANNOTATIONS[int(wc.step) - 1]["info_field"],
            annotation_prefix=lambda wc: EXTRA_ANNOTATIONS[int(wc.step) - 1]["annotation_prefix"],
        conda:
            "../envs/snpeff.yaml"
        shell:
            r"""
            echo "Starting annotation_step {wildcards.step} at: $(date)" >> {log}
            SnpSift {params.java_opts} annotate \
              -info {params.info_field} \
              -name {params.annotation_prefix} \
              {params.vcf_file} {input} \
            | bgzip -c > {output} 2>> {log}
            echo "Finished annotation_step {wildcards.step} at: $(date)" >> {log}
            """


rule finalize_annotation:
    input:
        get_final_annotation_input,
    output:
        temp(
            os.path.join(
                ANNOTATION_DIR,
                "{sample}.{scatter_unit}.annotated.vcf.gz",
            )
        ),
    log:
        os.path.join(LOG_DIR, "finalize_annotation.{sample}.{scatter_unit}.log"),
    conda:
        "../envs/snpeff.yaml"
    shell:
        r"""
        echo "Finalizing annotation at: $(date)" > {log}
        cp {input} {output}
        echo "Done at: $(date)" >> {log}
        """
