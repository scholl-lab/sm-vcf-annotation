"""Conditional scatter/gather rules for interval-based parallelization."""

if SCATTER_MODE == "interval":

    INTERVALS_LIST = os.path.join(INTERVALS_DIR, "intervals.interval_list")

    def _get_noncanonical_contigs(dict_file, canonical):
        """Parse reference dictionary to find non-canonical contigs."""
        noncanonical = []
        with open(dict_file) as f:
            for line in f:
                if line.startswith("@SQ"):
                    for field in line.strip().split("\t"):
                        if field.startswith("SN:"):
                            contig = field.replace("SN:", "")
                            if contig not in canonical:
                                noncanonical.append(contig)
        return noncanonical

    _noncanonical = _get_noncanonical_contigs(REF_DICT, CANONICAL_CONTIGS)
    _exclude_intervals = " ".join(["-XL " + c for c in _noncanonical])

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
            exclude_intervals=_exclude_intervals,
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
                {params.exclude_intervals} \
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
            temp(
                os.path.join(
                    ANNOTATION_DIR,
                    "{sample}.{scatter_unit}.vcf.gz",
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
              -O {output} 2>> {log}
            echo "Finished scatter_vcf at: $(date)" >> {log}
            """

    rule concatenate_annotated_vcfs:
        input:
            lambda wc: expand(
                os.path.join(
                    ANNOTATION_DIR,
                    "{sample}.{scatter_unit}.annotated.vcf.gz",
                ),
                sample=[wc.sample],
                scatter_unit=SCATTER_UNITS,
            ),
        output:
            os.path.join(ANNOTATION_DIR, "{sample}.annotated.vcf.gz"),
        log:
            os.path.join(LOG_DIR, "concatenate_annotated_vcfs.{sample}.log"),
        conda:
            "../envs/snpeff.yaml"
        shell:
            r"""
            echo "Starting concatenate_annotated_vcfs at: $(date)" >> {log}
            bcftools concat -Oz --threads {threads} {input} \
                -o {output} 2>> {log}
            bcftools index --threads {threads} -t {output} 2>> {log}
            echo "Finished concatenate_annotated_vcfs at: $(date)" >> {log}
            """
