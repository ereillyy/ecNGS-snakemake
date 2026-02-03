rule filter_read_cons:
    input:
        bam="{pipeline}/tmp/2_bam_filtered/b_blacklist/{sample}.bam"
    output:
        filtered=temp("{pipeline}/tmp/2_bam_filtered/c_params/{sample}.bam")
    threads: 1
    resources:
        mem_mb=5 * 1024,
        time="01:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/logs/filter_read_cons/{sample}.log"
    params:
        subs_threshold=2,
        as_xs_diff_threshold=10,
        n_count_threshold=15,
        indel_count_threshold=0,
        indel_length_threshold=10000000000
    shell:
        r"""
        echo "[filter_read_cons] Starting for {wildcards.sample}" >> {log}
        python ../../src/filter_bam_pysam.py \
            --input-bam {input.bam} \
            --output-bam {output.filtered} \
            --subs-threshold {params.subs_threshold} \
            --as-xs-diff-threshold {params.as_xs_diff_threshold} \
            --n-count-threshold {params.n_count_threshold} \
            --indel-count-threshold {params.indel_count_threshold} \
            --indel-length-threshold {params.indel_length_threshold} \
            > {output.filtered}.out \
            2>> {log}
        samtools stats {output.filtered} > {output.filtered}.stats
        samtools index {output.filtered}
        echo "[filter_read_cons] Finished for {wildcards.sample}" >> {log}
        """
