"""Conditional scatter/gather rules for interval-based parallelization."""

if SCATTER_MODE == "none":

    rule rename_no_scatter:
        input:
            os.path.join(
                ANNOTATION_DIR,
                "{sample}.all.annotated.vcf.gz",
            ),
        output:
            os.path.join(ANNOTATION_DIR, "{sample}.annotated.vcf.gz"),
        log:
            os.path.join(LOG_DIR, "rename_no_scatter.{sample}.log"),
        conda:
            "../envs/snpeff.yaml"
        shell:
            r"""
            echo "Renaming final annotation at: $(date)" > {log}
            cp {input} {output}
            bcftools index --threads {threads} -t {output} 2>> {log}
            echo "Done at: $(date)" >> {log}
            """

elif SCATTER_MODE == "chromosome":

    rule scatter_vcf_chromosome:
        input:
            vcf_file=get_input_vcf,
            vcf_tbi=lambda wc: f"{get_input_vcf(wc)}.tbi",
        output:
            temp(os.path.join(ANNOTATION_DIR, "{sample}.{scatter_unit}.vcf.gz")),
        params:
            region_arg=lambda wc: (
                f"-t ^{','.join(CANONICAL_CONTIGS)}"
                if wc.scatter_unit == "other"
                else f"-r {wc.scatter_unit}"
            ),
        log:
            os.path.join(LOG_DIR, "scatter_vcf_chromosome.{sample}.{scatter_unit}.log"),
        conda:
            "../envs/snpeff.yaml"
        shell:
            r"""
            echo "Starting scatter_vcf_chromosome at: $(date)" >> {log}
            vcf_contig=$(bcftools view -h {input.vcf_file} | grep -m1 '^##contig=<ID=' | sed 's/^##contig=<ID=\([^,>]*\).*/\1/' || true)
            if [ -n "$vcf_contig" ] && [ "{wildcards.scatter_unit}" != "other" ]; then
                case "{wildcards.scatter_unit}:$vcf_contig" in
                    chr*:chr*) ;;
                    chr*:*|*:chr*)
                        echo "Error: scatter_unit {wildcards.scatter_unit} chr-prefix mismatches VCF contig $vcf_contig" >> {log}
                        exit 1
                        ;;
                esac
            fi
            bcftools view {params.region_arg} {input.vcf_file} \
                -Oz --threads {threads} -o {output} 2>> {log}
            echo "Finished scatter_vcf_chromosome at: $(date)" >> {log}
            """

elif SCATTER_MODE == "interval":

    INTERVALS_LIST = os.path.join(INTERVALS_DIR, "intervals.interval_list")

    _interval_ids = [f"{i:04d}" for i in range(SCATTER_COUNT)]

    rule scatter_intervals_by_ns:
        input:
            reference=REF,
        output:
            intervals_list=INTERVALS_LIST,
        log:
            os.path.join(LOG_DIR, "scatter_intervals_by_ns.log"),
        conda:
            "../envs/gatk.yaml"
        shell:
            r"""
            set -e
            echo "Starting ScatterIntervalsByNs at: $(date)" > {log}
            gatk ScatterIntervalsByNs \
                -R {input.reference} \
                -O {output.intervals_list} &>> {log}
            echo "Finished ScatterIntervalsByNs at: $(date)" >> {log}
            """

    rule split_intervals:
        input:
            intervals_list=INTERVALS_LIST,
            reference=REF,
        output:
            interval_files=expand(
                os.path.join(
                    INTERVALS_DIR,
                    "{interval_id}-scattered.interval_list",
                ),
                interval_id=_interval_ids,
            ),
        params:
            scatter_count=SCATTER_COUNT,
            intervals_dir=INTERVALS_DIR,
        log:
            os.path.join(LOG_DIR, "split_intervals.log"),
        conda:
            "../envs/gatk.yaml"
        shell:
            r"""
            set -e
            echo "Starting SplitIntervals at: $(date)" > {log}
            gatk SplitIntervals \
                -R {input.reference} \
                -L {input.intervals_list} \
                --scatter-count {params.scatter_count} \
                -O {params.intervals_dir} &>> {log}
            echo "Finished SplitIntervals at: $(date)" >> {log}
            """

    rule scatter_vcf:
        input:
            vcf_file=get_input_vcf,
            interval_file=lambda wc: os.path.join(
                INTERVALS_DIR,
                f"{wc.scatter_unit}.interval_list",
            ),
        output:
            vcf=temp(
                os.path.join(
                    ANNOTATION_DIR,
                    "{sample}.{scatter_unit}.vcf.gz",
                )
            ),
            tbi=temp(
                os.path.join(
                    ANNOTATION_DIR,
                    "{sample}.{scatter_unit}.vcf.gz.tbi",
                )
            ),
        log:
            os.path.join(LOG_DIR, "scatter_vcf.{sample}.{scatter_unit}.log"),
        params:
            java_opts=get_java_opts,
            ref=REF,
        conda:
            "../envs/gatk.yaml"
        shell:
            r"""
            echo "Starting scatter_vcf at: $(date)" >> {log}
            gatk --java-options "{params.java_opts}" SelectVariants \
              -R {params.ref} \
              -V {input.vcf_file} \
              -L {input.interval_file} \
              --ignore-non-starting-in-interval \
              -O {output.vcf} 2>> {log}
            echo "Finished scatter_vcf at: $(date)" >> {log}
            """


if SCATTER_MODE in ("chromosome", "interval"):

    rule concatenate_annotated_vcfs:
        input:
            vcfs=lambda wc: expand(
                os.path.join(
                    ANNOTATION_DIR,
                    "{sample}.{scatter_unit}.annotated.vcf.gz",
                ),
                sample=[wc.sample],
                scatter_unit=SCATTER_UNITS,
            ),
            raw_vcf=get_input_vcf,
        output:
            os.path.join(ANNOTATION_DIR, "{sample}.annotated.vcf.gz"),
        log:
            os.path.join(LOG_DIR, "concatenate_annotated_vcfs.{sample}.log"),
        conda:
            "../envs/snpeff.yaml"
        shell:
            r"""
            echo "Starting concatenate_annotated_vcfs at: $(date)" >> {log}
            bcftools concat -Oz --threads {threads} {input.vcfs} \
                -o {output} 2>> {log}
            bcftools index --threads {threads} -t {output} 2>> {log}
            [ "$(bcftools view -H {input.raw_vcf} | wc -l)" -eq "$(bcftools index -n {output})" ]
            echo "Finished concatenate_annotated_vcfs at: $(date)" >> {log}
            """
