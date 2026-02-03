rule mn_remove_dups:
    input:
        bam="{pipeline}/tmp/2_mn/b_aligned/{sample}.bam",
        bam_bai="{pipeline}/tmp/2_mn/b_aligned/{sample}.bam.bai"
    output:
        bam="{pipeline}/tmp/2_mn/c_dedup/{sample}.bam",
        bam_bai="{pipeline}/tmp/2_mn/c_dedup/{sample}.bam.bai",
        metrics="{pipeline}/qc/mn_rmdups/{sample}_metrics.txt"
    wildcard_constraints:
        sample="|".join(MATCHED_NORMALS)
    threads: 4
    resources:
        mem_mb=16 * 1024,
        time="04:00:00"
    conda:
        "../../../../../envs/main.yaml"
    log:
        "{pipeline}/qc/mn_rmdups/{sample}.log"
    shell:
        r"""
        echo "[$(date)] Starting remove_dups for {wildcards.sample}" > {log}
        picard MarkDuplicates \
        I={input.bam} \
        O={output.bam} \
        M={output.metrics} \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        REMOVE_DUPLICATES=true \
        DUPLICATE_SCORING_STRATEGY=TOTAL_MAPPED_REFERENCE_LENGTH \
        MINIMUM_DISTANCE=300 \
        >> {log} 2>&1
        samtools index {output.bam} >> {log} 2>&1
        echo "[$(date)] Finished remove_dups for {wildcards.sample}" >> {log}
        """
